rule split_fasta:
    input:
        fasta=lambda wc: config["samples"][wc.sample],
    output:
        fasta=temp(
            scatter.fasta(
                os.path.join(OUT_DIR, "{{sample}}/RepeatMasker/{scatteritem}/{scatteritem}.fa")
            )
        ),
    resources:
        mem=config.get("mem", 8),
    threads: 1
    conda:
        "../envs/env.yml"
    log:
        os.path.join(LOG_DIR, "{sample}/fasta/split.log"),
    script:
        "../scripts/split_fasta.py"


rule unzip_fasta:
    input:
        fasta=lambda wc: config["samples"][wc.sample],
    output:
        fasta=temp(os.path.join(OUT_DIR, "unzipped/{sample}.fasta")),
        fai=temp(os.path.join(OUT_DIR, "unzipped/{sample}.fasta.fai")),
    resources:
        mem=config.get("mem", 8),
    threads: 1
    conda:
        "../envs/env.yml"
    log:
        os.path.join(LOG_DIR, "{sample}/fasta/unzip.log"),
    shell:
        """
        seqtk seq -l 60 {input.fasta} > {output.fasta}
        samtools faidx {output.fasta}
        """
