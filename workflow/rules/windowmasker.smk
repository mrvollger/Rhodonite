rule dust_count:
    input:
        ref=rules.unzip_fasta.output.fasta,
    output:
        counts=temp("results/{sample}/windowmasker/dust.counts"),
    resources:
        mem=config.get("mem", 16),
    threads: 1
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/windowmasker/dust_count.log",
    params:
        mem=(config.get("mem", 16) - 2) * 1000,
    shell:
        """
        windowmasker -mem {params.mem} -mk_counts -in {input.ref} -out {output.counts}  
        """


rule windowmasker:
    input:
        counts=rules.dust_count.output.counts,
        ref=rules.unzip_fasta.output.fasta,
        fai=rules.unzip_fasta.output.fai,
    output:
        intervals=temp("results/{sample}/windowmasker/dust.intervals"),
        bed=temp("results/{sample}/windowmasker/dust.bed.gz"),
    resources:
        mem=config.get("mem", 16),
    threads: 1
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/windowmasker/bed.log",
    params:
        s_dir=S_DIR,
    shell:
        """
        windowmasker -ustat {input.counts} -dust true -in {input.ref} \
            -out {output.intervals} 

        {params.s_dir}/scripts/dust_to_bed.py -i {output.intervals} \
            | bedtools sort -header -g {input.fai} -i - \
            | gzip -c \
            > {output.bed} 
        """
