import os
import glob
import argparse
from pandas import read_csv
from math import floor
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command")
    
    parser_a = subparsers.add_parser("get-blocks", help="Extracts genome blocks of a specified length (bp) from each interval of each sequence in a FASTA file, \
                                                    e.g. first 50,000 bp as one block for every 500,000 bp in each sequence.")
    parser_a.add_argument("-p", "--prefix", metavar="", help="(Optional) Prefix for accessions.")
    parser_a.add_argument("-n", "--name", metavar="", required=True, type=str, help="Accession name.")
    parser_a.add_argument("-i", "--indir", metavar="", required=True, type=str, help="Absolute path where directories for each sample containing the corresponding bam file is located.")
    parser_a.add_argument("-S", "--blksize", metavar="", required=True, type=int, help="Size of genome block (bp).")
    parser_a.add_argument("-I", "--blkinv", metavar="", required=True, type=int, help="Size of interval where the first x bp (specified in -S/--blksize) will be extracted.")
    parser_a.set_defaults(func=extract_for_sample)
    
    parser_b = subparsers.add_parser("block-txt", help="Write a txt file containing the names of all possible blocks extracted from a genome index file (.fai) to the specified output path.")
    parser_b.add_argument("-f", "--fai", metavar="", required=True, type=str, help="Absolute path of input fai file.")
    parser_b.add_argument("-o", "--outpath", metavar="", required=True, type=str, help="Output path to generate a txt file in which the genome block names will be listed.")
    parser_b.add_argument("-S", "--blksize", metavar="", required=True, type=int, help="Size of genome block (bp).")
    parser_b.add_argument("-I", "--blkinv", metavar="", required=True, type=int, help="Size of interval where the first x bp (specified in -S/--blksize) will be extracted.")
    parser_b.set_defaults(func=generate_txt_block_names)
    
    parser_c = subparsers.add_parser("full-seq-aln", help="Merge (aligned) genome blocks column-to-column across all samples to create a FASTA file containing the final multiple sequence alignment.")
    parser_c.add_argument("-s", "--sampdir", metavar="", required=True, help="Path to a directory listing names of all samples.")
    parser_c.add_argument("-b", "--blktxt", metavar="", required=True, help="Absolute path of block name txt file.")
    parser_c.add_argument("-i", "--indir", metavar="", required=True, help="Input path of directory containing all (aligned) genome blocks fasta files.")
    parser_c.add_argument("-o", "--outdir", metavar="", required=True, help="Output path to where the full sequence alignment will be generated.")
    
    
    args = parser.parse_args()
    if args.command == "get-blocks":
        extract_for_sample(args.name, args.indir, args.blksize, args.blkinv, args.prefix)
    if args.command == "block-txt":
        generate_txt_block_names(args.fai, args.outpath, args.blksize, args.blkinv)
    if args.command == "full-seq-aln":
        merge_seq_per_samp(args.sampdir, args.blktxt, args.indir, args.outdir)
    

def chop_fasta_into_seq(fasta_file, prefix="", only_bp_more_than=0, outpath=""):
    """
    Chops a fasta file up into multiple files, one for each sequence in the original fasta file.
    Creates a directory called "seq_fasta" within the same directory as the input fasta to store all output files.
    Output file is named after the header of each sequence.
    Prefix option is provided for prefixing the output file name, e.g. sample identifiers.
    Sequence length sequence can be set to ignore sequences shorter than the threshold.
    """
    if outpath == "":
        subdir = os.path.join(os.path.dirname(fasta_file), "seq_fasta")
    else:
        subdir = outpath
    if not os.path.exists(subdir):
            os.mkdir(subdir)
    for seq in SeqIO.parse(fasta_file, "fasta"):
        if len(seq.seq) < only_bp_more_than:
            pass
        else:
            if prefix != "":
                seq.id = prefix + "." + seq.id
            seq.description = seq.id
            filename = seq.id + ".fasta"
            filepath = os.path.join(subdir,filename)
            SeqIO.write(seq, filepath, "fasta")


def extract_genom_blks(seq_fasta_file, blk_size, blk_interval, outpath=""):
    """
    Extracts genome blocks of a specified size (bp) from every interval from the fasta file. 
    x = size of block; y = size of interval.
    Only first x bp of every interval of y bp will be extracted.
    Output files will be stored in directory "seq_blocks_fasta", organised in subdirectories named after each scaffold.
    """
    if outpath == "":
        subdir = os.path.join(os.path.dirname(seq_fasta_file), "..", "seq_blocks_fasta")
    else:
        subdir = outpath
    if not os.path.exists(subdir):
        os.mkdir(subdir)
    for seq in SeqIO.parse(seq_fasta_file, "fasta"):
        nblock = floor(len(seq.seq)/blk_interval)
        for i in range(nblock+1):
            seq_blk = seq.seq[i*blk_interval:(i*blk_interval)+blk_size]
            if len(seq_blk) < 50000:
                break
            seq_blk_id = seq.id + f".{i*blk_interval}-{(i*blk_interval)+blk_size}"
            filename = seq_blk_id + ".fasta"
            filepath = os.path.join(subdir, filename)
            record = SeqRecord(Seq(str(seq_blk)), id=seq_blk_id, description="")
            SeqIO.write(record, filepath, "fasta")


def extract_for_sample(identifier, sample_bam_dir, blk_size, blk_interval, prefix=None):
    if prefix != None:
        iden = prefix + "." + identifier
    else:
        iden = identifier
    iden_fasta_path = os.path.join(sample_bam_dir, identifier, "fasta", iden+".fa")
    chop_fasta_into_seq(iden_fasta_path, iden, blk_size)
    
    seq_fasta = glob.glob(os.path.join(sample_bam_dir, identifier, "fasta", "seq_fasta", "*"))
    for seq in seq_fasta:
        extract_genom_blks(seq, blk_size, blk_interval)


def generate_txt_block_names(genom_fai, txt_outpath, blk_size, blk_interval):
    df = read_csv(genom_fai, sep="\t", usecols=[0,1], names=["name","length"])
    filt_df = df[df['length'] > blk_size]
    block_name_list = []
    for n in filt_df['name']:
        k = filt_df[filt_df['name']==n]['length']
        nblock = floor(k/blk_interval)
        for i in range(nblock+1):
            if int(k)-(i*blk_interval) < blk_size:
                break
            block_name = f"{n}.{i*blk_interval}-{(i*blk_interval)+blk_size}"
            block_name_list.append(block_name)
    with open(os.path.join(txt_outpath,"block_names.txt"), "w") as txt:
        for r in block_name_list:
            txt.write(r+"\n")
    print(f"Done! A total of {len(block_name_list)} genome blocks from all scaffolds with size >{blk_size} bp were found.")


def merge_seq_per_samp(sample_dir, blk_name_txt, indir, outdir):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    outpath = os.path.join(outdir, "full_seq_aln.fasta")
    with open(blk_name_txt) as txt:
        blk_list = [line.strip() for line in txt.readlines()]
    sample_list = os.listdir(sample_dir)
    record_all_samp = []
    for samp in sample_list:
        full_samp_seq = ""
        for blk in blk_list:
            abspath = os.path.join(indir, blk+".aligned.fasta")
            for seq in SeqIO.parse(abspath, "fasta"):
                if samp in seq.id:
                    full_samp_seq += seq.seq
                else:
                    pass
        record = SeqRecord(Seq(str(full_samp_seq)), id=samp, description="")
        record_all_samp.append(record)
        print(f"\nSequences for sample {samp} succesfully merged and written to {outpath}! Length: {len(full_samp_seq)}.")
    SeqIO.write(record_all_samp, outpath, "fasta")
    print("\nDone!")



if __name__ == "__main__":
    main()