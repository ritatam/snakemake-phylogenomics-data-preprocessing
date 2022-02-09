import os
import glob
configfile: "genom_blks_config.yaml"

SAMPLES = os.listdir(config["INPUT_BAM_DIR"])


def grep_block_name_list(genome_block_txt):
    with open(genome_block_txt) as txt:
        return [line.strip() for line in txt.readlines()]


def print_samp_block_name(hap, sample_list, block_name_list, basepath):
    blkdict = {}
    for block in block_name_list:
        blkdict[block] = []
        for samp in sample_list:
            filename = os.path.join(basepath, samp, "fasta/seq_blocks_fasta", f"{hap}.{samp}.{block}.fasta")
            blkdict[block].append(filename)
    return blkdict


def get_unaligned_fa(wildcards):
    checkpoint_output = checkpoints.concat_sample_per_block.get(**wildcards).output[0]
    file_names = expand(os.path.join(config["ALIGNED_DIR"], "{block_name}.aligned.fasta"), \
                        block_name=glob_wildcards(os.path.join(checkpoint_output, "{block_name}.fasta")).block_name)
    return file_names
    


rule all:
    input:
        expand("{path}/{sample}/fasta/{hap}.{sample}.fa", path=config["INPUT_BAM_DIR"], sample=SAMPLES, hap=config["HAPLOTYPE"]),
        expand("{path}/{sample}/fasta/seq_blocks_fasta", path=config["INPUT_BAM_DIR"], sample=SAMPLES),
        config["BASE_DIR"] + "/block_names.txt",
        config["UNALIGNED_DIR"],
        get_unaligned_fa,
        config["FULL_ALN_DIR"] + "/full_seq_aln.fasta"
        
        

rule angsd_bam2fasta:
    params:
        env = "pycodes",
        prefix = config["INPUT_BAM_DIR"] + "/{sample}/fasta/" + config["HAPLOTYPE"] + ".{sample}",
        fagz = config["INPUT_BAM_DIR"] + "/{sample}/fasta/" + config["HAPLOTYPE"] + ".{sample}.fa.gz"
    input:
        bam = config["INPUT_BAM_DIR"] + "/{sample}/" + config["HAPLOTYPE"] + ".{sample}.bam"
    output:
        fa = config["INPUT_BAM_DIR"] + "/{sample}/fasta/" + config["HAPLOTYPE"] + ".{sample}.fa"
    shell:
        "set +eu "
        " && . $(conda info --base)/etc/profile.d/conda.sh "
        " && conda activate {params.env} "
        " && echo $CONDA_PREFIX; "
        
        "angsd -i {input.bam} -doFasta 2 -doCounts 1 -out {params.prefix}; "
        "gunzip {params.fagz}; "


rule extract_genome_blocks_fa:
    params:
        env = "pycodes",
        prefix = config["HAPLOTYPE"],
        sample = "{sample}"
    input:
        fa = config["INPUT_BAM_DIR"] + "/{sample}/fasta/" + config["HAPLOTYPE"] + ".{sample}.fa"
    output:
        outdir = directory(config["INPUT_BAM_DIR"] + "/{sample}/fasta/seq_blocks_fasta")
    shell:
        "set +eu "
        " && . $(conda info --base)/etc/profile.d/conda.sh "
        " && conda activate {params.env} "
        " && echo $CONDA_PREFIX; "
        
        "python {config[extract_genom_blks_path]} get-blocks -p {params.prefix} -n {params.sample} -i {config[INPUT_BAM_DIR]} \
        -S {config[genome_block_size]} -I {config[genome_block_interval]}"


rule genome_block_names_txt:
    params:
        env = "pycodes"
    output:
        config["BASE_DIR"] + "/block_names.txt"
    shell:
        "set +eu "
        " && . $(conda info --base)/etc/profile.d/conda.sh "
        " && conda activate {params.env} "
        " && echo $CONDA_PREFIX; "
        
        "python {config[extract_genom_blks_path]} block-txt -f {config[ref_fai_path]} -o {config[BASE_DIR]} -S {config[genome_block_size]} -I {config[genome_block_interval]}"


checkpoint concat_sample_per_block:
    input:
        blk_txt = config["BASE_DIR"] + "/block_names.txt",
        indir = expand(config["INPUT_BAM_DIR"] + "/{sample}/fasta/seq_blocks_fasta", sample=SAMPLES)
    output:
        outdir = directory(config["UNALIGNED_DIR"])
    run:
        if not os.path.exists(output.outdir):
            os.mkdir(output.outdir)
        block_name_list = grep_block_name_list(config["BASE_DIR"] + "/block_names.txt")
        block_dict = print_samp_block_name(config["HAPLOTYPE"], SAMPLES, block_name_list, config["INPUT_BAM_DIR"])
        for b in block_name_list:
            concat_block_path = os.path.join(output.outdir, b+".fasta")
            shell(f"cat {' '.join(block_dict[b])} >> {concat_block_path}")


def aggregate_block_name(wildcards):
    return glob.glob(config["UNALIGNED_DIR"] + "/{wildcards.block_name}.fasta")


rule mafft:
    params:
        env = "mafft"
    input:
        aggregate_block_name,
        unaligned_fa = config["UNALIGNED_DIR"] + "/{block_name}.fasta"
    output:
        aligned_fa = config["ALIGNED_DIR"] + "/{block_name}.aligned.fasta"
    shell:
        "set +eu "
        " && . $(conda info --base)/etc/profile.d/conda.sh "
        " && conda activate {params.env} "
        " && echo $CONDA_PREFIX; "

        "mafft --retree 1 --maxiterate 0 {input.unaligned_fa} > {output.aligned_fa}"


rule create_full_seq_aln:
    params:
        env = "pycodes"
    input:
        get_unaligned_fa,
        blk_txt = config["BASE_DIR"] + "/block_names.txt",
    output:
        full_aln = config["FULL_ALN_DIR"] + "/full_seq_aln.fasta"
    shell:
        "set +eu "
        " && . $(conda info --base)/etc/profile.d/conda.sh "
        " && conda activate {params.env} "
        " && echo $CONDA_PREFIX; "
        
        "python {config[extract_genom_blks_path]} full-seq-aln \
        -s {config[INPUT_BAM_DIR]} -b {input.blk_txt} -i {config[ALIGNED_DIR]} -o {config[FULL_ALN_DIR]}"
