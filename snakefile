import sys, glob, os

SAMPLE,=glob_wildcards("raw/{sample}_L001_R1_001.fastq.gz")
READ_ID=["1","2"]

## ---------------------------------------------------------------------------##
## Rule All
## ---------------------------------------------------------------------------##
rule all:
    input:
        # "results/mafft/rnaviralHMM.mafft.aln",
        # "results/mafft/rnaviral.mafft.aln",
        # "results/mafft/megahit.mafft.aln"
        expand("data/assembly/final/{sample}.rnaviral.fasta", sample=SAMPLE),
        expand("data/assembly/final/{sample}.rnaviralhmm.fasta",sample=SAMPLE),
        expand("data/assembly/final/{sample}.megahit.fasta",sample=SAMPLE)


## ---------------------------------------------------------------------------##
## Concatenating reads
## ---------------------------------------------------------------------------##
rule concatenation:
    input:
        L1="raw/{sample}_L001_R{read_id}_001.fastq.gz",
        L2="raw/{sample}_L002_R{read_id}_001.fastq.gz"
    output:
        temporary("raw/merged/{sample}_R{read_id}.cat.fastq.gz")
    shell:"""
        cat {input.L1} {input.L2} > {output}
    """

## ---------------------------------------------------------------------------##
## Read trimming
## ---------------------------------------------------------------------------##
rule bbduk:
    input:
        read1="raw/merged/{sample}_R1.cat.fastq.gz",
        read2="raw/merged/{sample}_R2.cat.fastq.gz",
        adapters="data/adapters.fa"
    output:
        out1="raw/trimmed/{sample}_R1.trimmed.cat.fastq.gz",
        out2="raw/trimmed/{sample}_R2.trimmed.cat.fastq.gz"
    params:
        stats="data/QC/{sample}/bbduk_{sample}_stats.txt"
    shell:"""
        bbduk.sh -Xmx2g in1={input.read1} \
        in2={input.read2} out1={output.out1} \
        out2={output.out2} stats={params.stats} \
        ref={input.adapters} \
        ktrim=r k=23 mink=11 hdist=1 qtrim=rl trimq=10 maq=10 minlen=50 \
        threads=48 tpe tbo
    """

## ---------------------------------------------------------------------------##
## rnavrialspades Assembly
## ---------------------------------------------------------------------------##
rule rna_spades:
    input:
        reads1="raw/trimmed/{sample}_R1.trimmed.cat.fastq.gz",
        reads2="raw/trimmed/{sample}_R2.trimmed.cat.fastq.gz"
    output:
        "data/assembly/rnaviral_spades/{sample}/contigs.fasta"
    params:
        outdir="data/assembly/rnaviral_spades/{sample}/"
    shell:"""
    rnaviralspades.py -1 {input.reads1} -2 {input.reads2} -o {params.outdir} -t 32
    """

## ---------------------------------------------------------------------------##
## rnavrialspades Assembly Custom HMMs
## ---------------------------------------------------------------------------##
rule rnaHMM_spades:
    input:
        reads1="raw/trimmed/{sample}_R1.trimmed.cat.fastq.gz",
        reads2="raw/trimmed/{sample}_R2.trimmed.cat.fastq.gz"
    output:
        "data/assembly/rnaviralHMM_spades/{sample}/scaffolds.fasta"
    params:
        outdir="data/assembly/rnaviralHMM_spades/{sample}/"
    shell:"""
    rnaviralspades.py -1 {input.reads1} -2 {input.reads2} -o {params.outdir} --custom-hmms HMM/hmms/ -t 32
    """
## ---------------------------------------------------------------------------##
## Megahit Assembly
## ---------------------------------------------------------------------------##
rule Megahit:
    input:
        reads1="raw/trimmed/{sample}_R1.trimmed.cat.fastq.gz",
        reads2="raw/trimmed/{sample}_R2.trimmed.cat.fastq.gz"
    output:
        "data/assembly/megahit/{sample}/final.contigs.fa"
    params:
        outdir="data/assembly/megahit/{sample}"
    shell:"""
    megahit -1 {input.reads1} -2 {input.reads2} -o {params.outdir} --force
    """

## ---------------------------------------------------------------------------##
## Assembly rename
## ---------------------------------------------------------------------------##
rule rename:
    input:
        rnaviralspades="data/assembly/rnaviral_spades/{sample}/contigs.fasta",
        rnaviralHMM_spades="data/assembly/rnaviralHMM_spades/{sample}/scaffolds.fasta",
        megahit="data/assembly/megahit/{sample}/final.contigs.fa"
    output:
        rnaviralspades="data/assembly/final/{sample}.rnaviral.fasta",
        rnaviralHMM_spades="data/assembly/final/{sample}.rnaviralhmm.fasta",
        megahit="data/assembly/final/{sample}.megahit.fasta"
    shell:"""
    cp {input.rnaviralspades} {output.rnaviralspades}
    cp {input.rnaviralHMM_spades} {output.rnaviralHMM_spades}
    cp {input.megahit} {output.megahit}
    """

# ## ---------------------------------------------------------------------------##
# ## concat refs and assemblies
# ## ---------------------------------------------------------------------------##
# rule assembly_cat:
#     input:
#         expand("data/assembly/final/{sample}.rnaviralhmm.fasta",sample=SAMPLE),
#         expand("data/assembly/final/{sample}.rnaviral.fasta",sample=SAMPLE),
#         expand("data/assembly/final/{sample}.megahit.fasta",sample=SAMPLE),
#         refs="data/WMV_refs.fasta"
#     output:
#         temporary("data/assembly/rnaviralHMM.cat.fasta"),
#         temporary("data/assembly/rnaviral.cat.fasta"),
#         temporary("data/assembly/megahit.cat.fasta")
#     params:
#         input_dir="data/assembly/final/"
#     shell:"""
#     cat {params.input_dir}*.rnaviralhmm.fasta {input.refs} > {output}
#     cat {params.input_dir}*.rnaviral.fasta {input.refs} > {output}
#     cat {params.input_dir}*.megahit.fasta {input.refs} > {output}
#     """

# ## ---------------------------------------------------------------------------##
# ## mafft spades HMM
# ## ---------------------------------------------------------------------------##
# rule mafft:
#     input:
#         hmm="data/assembly/rnaviralHMM.cat.fasta",
#         spades="data/assembly/rnaviral.cat.fasta",
#         megahit="data/assembly/megahit.cat.fasta"
#     output:
#         hmm="results/mafft/rnaviralHMM.mafft.aln",
#         spades="results/mafft/rnaviral.mafft.aln",
#         megahit="results/mafft/megahit.mafft.aln"
#     shell:"""
#     mafft {input.hmm} > {output.hmm}
#     mafft {input.spades} > {output.spades}
#     mafft {input.megahit} > {output.megahit}
#     """
