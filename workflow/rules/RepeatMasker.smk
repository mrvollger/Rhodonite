rule run_split_RepeatMasker:
    input:
        fasta="results/{sample}/RepeatMasker/{scatteritem}/{scatteritem}.fa",
    output:
        msk=temp(
            "results/{sample}/RepeatMasker/{scatteritem}/{scatteritem}.fa.masked"
        ),
        out=temp("results/{sample}/RepeatMasker/{scatteritem}/{scatteritem}.fa.out"),
        cat=temp("results/{sample}/RepeatMasker/{scatteritem}/{scatteritem}.fa.cat"),
        cat_all=temp(
            "results/{sample}/RepeatMasker/{scatteritem}/{scatteritem}.fa.cat.all"
        ),
        tbl=temp("results/{sample}/RepeatMasker/{scatteritem}/{scatteritem}.fa.tbl"),
    resources:
        mem=config.get("mem", 8),
    threads: config.get("threads", 8)
    conda:
        "envs/env.yml"
    log:
        "logs/{sample}/RepeatMasker/{scatteritem}.log",
    params:
        opts=config.get("RepeatMaskerOptions", "-s -xsmall -e ncbi"),
        species=config.get("RepeatMaskerSpecies", "human"),
    shell:
        """
        RepeatMasker \
            {params.opts} \
            -species {params.species} \
            -pa {threads} \
            {input.fasta}  2> {log}

        if [ -f "{output.msk}" ]; then
            echo "masked fasta exists"
        else 
            echo "No repeats found, linking unmasked fasta to masked fasta"
            ln {input.fasta} {output.msk}
        fi

        # in case of no repeats touch the output
        touch {output}
        """


rule make_RepeatMasker_bed:
    input:
        out=rules.run_split_RepeatMasker.output.out,
    output:
        bed=temp(
            "results/{sample}/RepeatMasker/{scatteritem}/{scatteritem}.fa.rm.bed"
        ),
    conda:
        "envs/env.yml"
    log:
        "logs/{sample}/RepeatMasker/{scatteritem}.bed.log",
    threads: 1
    shell:
        """
        RM2Bed.py -d $(dirname {input.out}) {input.out} {output.bed}
        """


rule RepeatMasker:
    input:
        out=gather.fasta(rules.run_split_RepeatMasker.output.out, allow_missing=True),
        bed=gather.fasta(rules.make_RepeatMasker_bed.output.bed, allow_missing=True),
        fai=lambda wc: f'{config[wc.sample]["ref"]}.fai',
    output:
        out="results/{sample}/RepeatMasker/RM.out",
        bed="results/{sample}/RepeatMasker/RM.bed.gz",
    threads: 1
    conda:
        "envs/env.yml"
    log:
        "logs/{sample}/RepeatMasker.log",
    shell:
        """
        cat {input.out} \
            | grep -v "^There" \
            | head -n 3\
            > {output.out}
        tail -q -n +4 {input.out} \
            | grep -v "^There" \
            >> {output.out}


        cat {input.bed} \
            | bedtools sort \
                -g {input.fai} -i - \
            | gzip -c \
            > {output.bed}
        """
