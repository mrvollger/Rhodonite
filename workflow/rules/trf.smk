rule run_split_trf:
    input:
        fasta=rules.run_split_RepeatMasker.input.fasta,
    output:
        dat=temp("results/{sample}/trf/{scatteritem}/{scatteritem}.dat"),
    resources:
        mem=config.get("mem", 16),
    threads: 1
    conda:
        "envs/env.yml"
    log:
        "logs/{sample}/trf/{scatteritem}.log",
    params:
        trf_opts=config.get("trf_opts", "2 7 7 80 10 50 15 -l 25"),
    shell:
        """
        trf {input.fasta} {params.trf_opts} -h -ngs > {output.dat}
        """


rule trf:
    input:
        dat=gather.fasta(rules.run_split_trf.output.dat, allow_missing=True),
        fai=lambda wc: f'{config[wc.sample]["ref"]}.fai',
    output:
        dat="results/{sample}/trf/trf.dat",
        bed="results/{sample}/trf/trf.bed.gz",
    threads: 1
    conda:
        "envs/env.yml"
    log:
        "logs/{sample}/trf.log",
    params:
        s_dir=S_DIR,
    shell:
        """
        cat {input.dat} > {output.dat}

        {params.s_dir}/scripts/trf_to_bed.py {input.dat} \
            | bedtools sort -header -g {input.fai} -i - \
            | gzip -c \
            > {output.bed}
        """
