
rule run_split_gaps:
    input:
        fasta=rules.run_split_RepeatMasker.input.fasta,
    output:
        bed=temp("results/{sample}/gaps/{scatteritem}/{scatteritem}.bed"),
    resources:
        mem=config.get("mem", 16),
    threads: 1
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/gaps/{scatteritem}.log",
    script:
        "../scripts/HardMaskToBed.py"


rule gaps:
    input:
        bed=gather.fasta(rules.run_split_gaps.output.bed, allow_missing=True),
        fai=lambda wc: f'{config["samples"][wc.sample]}.fai',
    output:
        bed="results/{sample}/gaps/gaps.bed.gz",
    resources:
        mem=config.get("mem", 8),
    threads: 1
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/gaps.log",
    shell:
        """
        cat {input.bed} \
            | bedtools sort \
                -g {input.fai} -i - \
            | gzip -c \
            > {output.bed}
        """
