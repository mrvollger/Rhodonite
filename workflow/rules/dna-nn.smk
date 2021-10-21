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
    threads: config.get("threads", 8)
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
    resources:
        mem=config.get("mem", 8),
    threads: 1
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/dna-brnn.log",
    shell:
        """
        cat {input.bed} | bedtools sort -g {input.fai} -i - | bgzip > {output.bed}
        """
