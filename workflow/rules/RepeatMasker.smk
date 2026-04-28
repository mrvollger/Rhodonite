# this rule is required so that mutliple repeatmasker setup rules will not try and run at once.
rule setup_RepeatMasker:
    output:
        build=temp(os.path.join(OUT_DIR, "RepeatMasker.software.built.txt")),
    resources:
        mem=config.get("mem", 8),
    threads: config.get("threads", 8)
    conda:
        "../envs/env.yml"
    log:
        os.path.join(LOG_DIR, "RepeatMasker.build.log"),
    params:
        opts=config.get("RepeatMaskerOptions", "-s -xsmall -e ncbi"),
        species=config.get("RepeatMaskerSpecies", "human"),
    shell:
        """
        printf \
            ">n\naaaaaaAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n" \
            > {resources.tmpdir}/build.fa

        RepeatMasker \
            {params.opts} \
            -species {params.species} \
            -pa {threads} \
            {resources.tmpdir}/build.fa  2> {log}
        touch {output.build}
        """


rule run_split_RepeatMasker:
    input:
        build=rules.setup_RepeatMasker.output.build,
        fasta=os.path.join(OUT_DIR, "{sample}/RepeatMasker/{scatteritem}/{scatteritem}.fa"),
    output:
        msk=temp(os.path.join(OUT_DIR, "{sample}/RepeatMasker/{scatteritem}/{scatteritem}.fa.masked")),
        out=temp(os.path.join(OUT_DIR, "{sample}/RepeatMasker/{scatteritem}/{scatteritem}.fa.out")),
        cat=temp(os.path.join(OUT_DIR, "{sample}/RepeatMasker/{scatteritem}/{scatteritem}.fa.cat")),
        cat_all=temp(
            os.path.join(OUT_DIR, "{sample}/RepeatMasker/{scatteritem}/{scatteritem}.fa.cat.all")
        ),
        tbl=temp(os.path.join(OUT_DIR, "{sample}/RepeatMasker/{scatteritem}/{scatteritem}.fa.tbl")),
    resources:
        mem=config.get("mem", 8),
    threads: config.get("threads", 8)
    conda:
        "../envs/env.yml"
    log:
        os.path.join(LOG_DIR, "{sample}/RepeatMasker/{scatteritem}.log"),
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
        bed=temp(os.path.join(OUT_DIR, "{sample}/RepeatMasker/{scatteritem}/{scatteritem}.fa.rm.bed")),
    resources:
        mem=config.get("mem", 8),
    conda:
        "../envs/env.yml"
    log:
        os.path.join(LOG_DIR, "{sample}/RepeatMasker/{scatteritem}.bed.log"),
    threads: 1
    shell:
        """
        RM2Bed.py -d $(dirname {input.out}) {input.out} {output.bed}
        """


rule RepeatMasker:
    input:
        out=gather.fasta(rules.run_split_RepeatMasker.output.out, allow_missing=True),
        bed=gather.fasta(rules.make_RepeatMasker_bed.output.bed, allow_missing=True),
        fai=lambda wc: f'{config["samples"][wc.sample]}.fai',
    output:
        out=os.path.join(OUT_DIR, "{sample}/RepeatMasker/RM.out"),
        bed=os.path.join(OUT_DIR, "{sample}/RepeatMasker/RM.bed.gz"),
    resources:
        mem=config.get("mem", 8),
    threads: 1
    conda:
        "../envs/env.yml"
    log:
        os.path.join(LOG_DIR, "{sample}/RepeatMasker.log"),
    shell:
        """
        set +o pipefail
        grep -v "^There" {input.out} | head -n 3 > {output.out}
        tail -q -n +4 {input.out} | grep -v "^There" >> {output.out}

        echo "made .out for RepeatMasker"
        set -o pipefail
        cat {input.bed} | bedtools sort -g {input.fai} -i - | gzip > {output.bed}
        """
