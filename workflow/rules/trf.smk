"""
TRF pipeline
"""

"""
TRF options 
Where: (all weights, penalties, and scores are positive)
  File = sequences input file
  Match  = matching weight
  Mismatch  = mismatching penalty
  Delta = indel penalty
  PM = match probability (whole number)
  PI = indel probability (whole number)
  Minscore = minimum alignment score to report
  MaxPeriod = maximum period size to report

recomended from TRF:
trf yourfile.fa 2 5 7 80 10 50 2000

previous trf options:
trf_opts=config.get("trf_opts", "2 7 7 80 10 50 15 -l 25"),
"""


rule run_split_trf:
    input:
        fasta=rules.run_split_RepeatMasker.input.fasta,
    output:
        dat=temp("results/{sample}/trf/{scatteritem}/{scatteritem}.dat"),
    threads: 1
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/trf/{scatteritem}.log",
    resources:
        mem=config.get("mem", 8),
    params:
        trf_opts=config.get("trf_opts", "2 5 7 80 10 50 500 -l 50"),
    shell:
        """
        trf {input.fasta} {params.trf_opts} -h -ngs > {output.dat}
        """


rule trf_bed:
    input:
        dat=gather.fasta(rules.run_split_trf.output.dat, allow_missing=True),
    output:
        bed=temp("temp/{sample}/trf/trf.bed"),
    threads: 1
    resources:
        mem=config.get("mem", 8),
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/trf.log",
    script:
        "../scripts/trf_to_bed.py"


rule trf:
    input:
        bed=rules.trf_bed.output.bed,
        dat=gather.fasta(rules.run_split_trf.output.dat, allow_missing=True),
        fai=lambda wc: f'{config["samples"][wc.sample]}.fai',
    output:
        dat="results/{sample}/trf/trf.dat",
        bed="results/{sample}/trf/trf.bed.gz",
    threads: 1
    resources:
        mem=config.get("mem", 8),
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/trf.log",
    shell:
        """
        cat {input.dat} > {output.dat}
        cat {input.bed} \
            | bedtools sort -header -g {input.fai} -i - \
            | gzip -c \
            > {output.bed}
        """
