rule split_fasta:
    input:
        fasta=lambda wc: config["samples"][wc.sample],
    output:
        fasta=temp(
            scatter.fasta(
                "results/{{sample}}/RepeatMasker/{scatteritem}/{scatteritem}.fa"
            )
        ),
    resources:
        mem=config.get("mem", 8),
    threads: 1
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/fasta/split.log",
    script:
        "../scripts/split_fasta.py"


rule unzip_fasta:
    input:
        fasta=lambda wc: config["samples"][wc.sample],
    output:
        fasta=temp("results/unzipped/{sample}.fasta"),
        fai=temp("results/unzipped/{sample}.fasta.fai"),
    resources:
        mem=config.get("mem", 8),
    threads: 1
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/fasta/unzip.log",
    shell:
        """
        gunzip -c {input.fasta} | seqtk seq -l 60 > {output.fasta}
        samtools faidx {output.fasta}
        """
