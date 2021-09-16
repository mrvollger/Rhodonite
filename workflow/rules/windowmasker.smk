rule dust_count:
    input:
        ref=rules.unzip_fasta.output.fasta,
    output:
        counts=temp("results/{sample}/windowmasker/dust.counts"),
    resources:
        mem=config.get("mem", 16),
    threads: 1
    conda:
        "../envs/blast.yml"
    log:
        "logs/{sample}/windowmasker/dust_count.log",
    params:
        mem=(config.get("mem", 16) - 2) * 1000,
    shell:
        """
        windowmasker -mem {params.mem} -mk_counts -in {input.ref} -out {output.counts}  
        """


rule run_windowmasker:
    input:
        counts=rules.dust_count.output.counts,
        ref=rules.unzip_fasta.output.fasta,
        fai=rules.unzip_fasta.output.fai,
    output:
        intervals=temp("temp/{sample}/windowmasker/dust.intervals"),
    resources:
        mem=config.get("mem", 16),
    threads: 1
    conda:
        "../envs/blast.yml"
    log:
        "logs/{sample}/windowmasker/intervals.log",
    shell:
        """
        windowmasker -ustat {input.counts} -dust true -in {input.ref} \
            -out {output.intervals} 
        """


rule run_windowmasker_bed:
    input:
        intervals=rules.run_windowmasker.output.intervals,
    output:
        bed=temp("temp/{sample}/windowmasker/dust.bed.gz"),
    resources:
        mem=config.get("mem", 16),
    threads: 1
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/windowmasker/bed.log",
    script:
        "../scripts/dust_to_bed.py"


rule windowmasker:
    input:
        intervals=rules.run_windowmasker.output.intervals,
        bed=rules.run_windowmasker_bed.output.bed,
        fai=lambda wc: f'{config["samples"][wc.sample]}.fai',
    output:
        bed="results/{sample}/windowmasker/dust.bed.gz",
        intervals="results/{sample}/windowmasker/dust.intervals",
    resources:
        mem=config.get("mem", 16),
    threads: 1
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/windowmasker/done.log",
    shell:
        """
        cp {input.intervals} {output.intervals}

        cat {input.bed} \
            | bedtools sort -header -g {input.fai} -i - \
            | gzip -c \
        > {output.bed} 
        """
