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
        "../envs/env.yml"
    log:
        "logs/{sample}/RepeatMasker/dup_masker_1.{scatteritem}.log",
    params:
        dupmasker=workflow.source_path("../scripts/DupMaskerParallel"),
        libs=workflow.source_path("../scripts/Libs/Libraries/dupliconlib.fa"),
    shell:
        """
        ln -s {params.libs} \
            $(dirname $(realpath $(which DupMasker)))/Libraries/. \
            || echo "duplicon lib already in place"

        {params.dupmasker} \
            -pa {threads} \
            -dupout \
            -engine ncbi \
            {input.fasta}
        """


rule setup_DupMasker:
    output:
        ready=temp("results/{sample}/RepeatMasker/DupMasker_ready.txt"),
    resources:
        mem=1,
    threads: 1
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/RepeatMasker/dup_masker_setup.log",
    params:
        libs=workflow.source_path("../scripts/Libs/Libraries/dupliconlib.fa"),
    shell:
        """
        cp -n {params.libs} \
            $(dirname $(realpath $(which DupMasker)))/Libraries/dupliconlib.fa \
            || echo "duplicon lib already in place"
        touch {output}
        """


rule run_DupMasker_step_2:
    input:
        ready=rules.setup_DupMasker.output.ready,
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
        "../envs/env.yml"
    log:
        "logs/{sample}/RepeatMasker/dup_masker_2.{scatteritem}.log",
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
        "../envs/env.yml"
    log:
        "logs/{sample}/RepeatMasker/dup_masker_3.{scatteritem}.log",
    params:
        DupMask_parserV6=workflow.source_path("../scripts/DupMask_parserV6.pl"),
        colors=workflow.source_path(
            "../scripts/DupMask_parserV6/r.all.repeat.10K.nr.Color"
        ),
        chrs=workflow.source_path("../scripts/DupMask_parserV6/bd35ChrSize"),
    shell:
        """
        perl {params.DupMask_parserV6} -c {params.colors} -L {params.chrs} -i {input.dup} -E -o {output.extra}
        """


rule DupMasker:
    input:
        extra=gather.fasta(rules.run_DupMasker_step_3.output.extra, allow_missing=True),
        fai=lambda wc: f'{config["samples"][wc.sample]}.fai',
    output:
        extra="results/{sample}/RepeatMasker/duplicons.extra",
        bed="results/{sample}/RepeatMasker/duplicons.bed.gz",
    resources:
        mem=4,
    threads: 1
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/RepeatMasker/DupMasker.log",
    params:
        DupMasker_bed9=workflow.source_path("../scripts/DupMasker_bed9.py"),
    shell:
        """
        head -n 1 {input.extra[0]} > {output.extra}
        tail -q -n+2 {input.extra} >> {output.extra}

        python {params.DupMasker_bed9} {output.extra} \
            | bedtools sort -header -g {input.fai} -i - \
            | gzip -c \
            > {output.bed}
        """
