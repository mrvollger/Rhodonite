rule make_dna_brnn:
    output:
        dna_brnn="temp/dna-nn/dna-brnn",
        model="temp/dna-nn/models/attcc-alpha.knm",
    resources:
        mem=1,
    threads: 1
    conda:
        "../envs/env.yml"
    log:
        "logs/make-dna-brnn.log",
    shell:
        """
        rm -rf temp/temp/dna-nn
        git clone https://github.com/lh3/dna-nn temp/temp/dna-nn
        pushd temp/temp/dna-nn
        make 
        popd
        cp temp/temp/dna-nn/dna-brnn {output.dna_brnn}
        cp temp/temp/dna-nn/models/attcc-alpha.knm {output.model}
        chmod +x {output.dna_brnn}
        rm -rf temp/temp/dna-nn
        """


rule run_split_dna_brnn:
    input:
        fasta=rules.run_split_RepeatMasker.input.fasta,
        dna_brnn=rules.make_dna_brnn.output.dna_brnn,
        model=rules.make_dna_brnn.output.model,
    output:
        bed=temp("results/{sample}/dna-brnn/{scatteritem}.bed"),
    resources:
        mem=config.get("mem", 8),
    threads: config.get("threads", 4)
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/dna-brnn/{scatteritem}.log",
    shell:
        """
        {input.dna_brnn} \
            -Ai {input.model} \
            -t {threads} \
            {input.fasta} \
        > {output.bed} 
        """


rule dna_brnn:
    input:
        bed=gather.fasta(rules.run_split_dna_brnn.output.bed, allow_missing=True),
        fai=lambda wc: f'{config["samples"][wc.sample]}.fai',
    output:
        bed="results/{sample}/dna-brnn/dna-brnn.bed.gz",
        bed9="results/{sample}/dna-brnn/dna-brnn.bed9.gz",
    resources:
        mem=config.get("mem", 8),
    threads: 1
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/dna-brnn.log",
    params:
        red="255,0,0",
        blue="0,0,255",
    shell:
        """
        cat {input.bed} | bedtools sort -g {input.fai} -i - | bgzip > {output.bed}

        ( 
            printf '#ct\\tst\\ten\\tname\\tscore\\tstrand\\ttst\\tten\\tcolor\\n'; \
            bgzip -dc {output.bed} \
                | sed 's/2$/{params.red}/g' \
                | sed 's/1$/{params.blue}/g' \
                | awk -v OFS=$'\\t' '{{print $1,$2,$3,"n","score",".",$2,$3,$4 }}' \
        ) \
            | gzip -c > {output.bed9}
        """
