from snakemake.utils import min_version

min_version("6.0")

# define your configuration
config = {
    "samples": {
        "test": "ref.fasta",
        "test2": "ref2.fasta",
        "test3": "ref3.fasta",
        "test4": "ref4.fasta",
    }
}

# or define a yaml e.g. see
# configfile: "config/config.yaml"


module Rhodonite:
    snakefile:
        "https://github.com/mrvollger/Rhodonite/raw/master/workflow/Snakefile"
    config:
        config


# use RepeatMasker rules from Rhodonite
use rule RepeatMasker from Rhodonite with:
    output:
        # You rename the output to anything you want
        # but maintain the order and keep "{sample}" in the name.
        # Every rule will have the standard output of the program
        # and a bed output that is sorted in the order of the ref.
        out="{sample}.rm.out",
        bed="{sample}.rm.bed.gz",


# use trf rules from Rhodonite
use rule trf from Rhodonite with:
    output:
        dat="{sample}.trf.dat",
        bed="{sample}.trf.bed.gz",


# use windowmasker rules from Rhodonite
use rule windowmasker from Rhodonite with:
    output:
        intervals="{sample}.windowmasker.dat",
        bed="{sample}.windowmasker.bed.gz",


# use DupMasker rules from Rhodonite
use rule DupMasker from Rhodonite with:
    output:
        extra="{sample}.duplicons.extra",
        bed="{sample}.duplicons.bed.gz",


rule all:
    input:
        expand(rules.RepeatMasker.output, sample=config["samples"].keys()),
        expand(rules.trf.output, sample=config["samples"].keys()),
        expand(rules.windowmasker.output, sample=config["samples"].keys()),
        expand(rules.DupMasker.output, sample=config["samples"].keys()),
