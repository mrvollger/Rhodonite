include: "RepeatMasker.smk"


rule run_DupMasker_step_1:
    input:
        fasta=rules.run_split_RepeatMasker.input.fasta,
        out=rules.run_split_RepeatMasker.output.out,
    output:
        dup=temp("temp/{sample}/RepeatMasker/{scatteritem}/{scatteritem}.fa.dupout"),
    resources:
        mem=config.get("mem", 16),
    threads: config.get("threads", 16)
    conda:
        "envs/env.yml"
    log:
        "logs/{sample}/RepeatMasker/dup_masker_1.{scatteritem}.log",
    params:
        s_dir=S_DIR,
    shell:
        """
        {params.s_dir}/scripts/DupMaskerParallel \
            -pa {threads} \
            -dupout \
            -engine ncbi \
            {input.fasta}
        """


rule run_DupMasker_step_2:
    input:
        fasta=rules.run_split_RepeatMasker.input.fasta,
        out=rules.run_split_RepeatMasker.output.out,
    output:
        dup=temp(
            "results/{sample}/RepeatMasker/{scatteritem}/{scatteritem}.fa.duplicons"
        ),
        dup_out=temp(
            "results/{sample}/RepeatMasker/{scatteritem}/{scatteritem}.fa.dupout"
        ),
        dup_all=temp(
            "results/{sample}/RepeatMasker/{scatteritem}/{scatteritem}.fa.dup.tmpmask.cat.all"
        ),
    resources:
        mem=config.get("mem", 8),
    threads: config.get("threads", 4)
    conda:
        "envs/env.yml"
    log:
        "logs/{sample}/RepeatMasker/dup_masker_2.{scatteritem}.log",
    params:
        s_dir=S_DIR,
    shell:
        """
        DupMasker \
            -engine ncbi \
            {input.fasta}
        touch {output}
        """


rule run_DupMasker_step_3:
    input:
        dup=rules.run_DupMasker_step_2.output.dup,
    output:
        extra=temp(
            "results/{sample}/RepeatMasker/{scatteritem}/{scatteritem}.fa.duplicons.extra"
        ),
    resources:
        mem=4,
    threads: 1
    conda:
        "envs/env.yml"
    log:
        "logs/{sample}/RepeatMasker/dup_masker_3.{scatteritem}.log",
    params:
        s_dir=S_DIR,
    shell:
        """
        {params.s_dir}/scripts/DupMask_parserV6.pl -i {input.dup} -E -o {output.extra}
        """


rule DupMasker:
    input:
        extra=gather.fasta(
            rules.run_DupMasker_step_3.output.extra, allow_missing=True
        ),
        fai=lambda wc: f'{config[wc.sample]["ref"]}.fai',
    output:
        extra="results/{sample}/RepeatMasker/duplicons.extra",
        bed9="results/{sample}/RepeatMasker/duplicons.bed.gz",
    resources:
        mem=4,
    threads: 1
    conda:
        "envs/env.yml"
    log:
        "logs/{sample}/RepeatMasker/DupMasker.log",
    params:
        s_dir=S_DIR,
    shell:
        """
        head -n 1 {input.extra[0]} > {output.extra}
        tail -q -n+2 {input.extra} >> {output.extra}

        {params.s_dir}/scripts/DupMasker_bed9.py {output.extra} {output.bed9}
        """
