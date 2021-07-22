from snakemake.utils import min_version

min_version("6.0")

# define your configuration
config = {
    "samples": {
        "test": ".test/test.fa.gz",
        "test2": ".test/test.fa.gz",
    }
}
# or define a yaml e.g. see
# configfile: "config/config.yaml"


module Rhodonite:
    snakefile:
        "../workflow/Snakefile"
        #"https://github.com/mrvollger/Rhodonite/raw/master/workflow/Snakefile"
    config:
        config


# import the rules from Rhodonite
use rule * from Rhodonite as Rhodonite_*


# use RepeatMasker rules from Rhodonite
use rule RepeatMasker from Rhodonite as Rhodonite_RepeatMasker with:
    output:
        # You rename the output to anything you want
        # but maintain the order and keep "{sample}" in the name.
        # Every rule will have the standard output of the program
        # and a bed output that is sorted in the order of the ref.
        out="{sample}.rm.out",
        bed="{sample}.rm.bed.gz",


# use trf rules from Rhodonite
use rule trf from Rhodonite as Rhodonite_trf with:
    output:
        dat="{sample}.trf.dat",
        bed="{sample}.trf.bed.gz",


# use windowmasker rules from Rhodonite
use rule windowmasker from Rhodonite as Rhodonite_windowmasker with:
    output:
        intervals="{sample}.windowmasker.dat",
        bed="{sample}.windowmasker.bed.gz",


# use DupMasker rules from Rhodonite
use rule DupMasker from Rhodonite as Rhodonite_DupMasker with:
    output:
        extra="{sample}.duplicons.extra",
        bed="{sample}.duplicons.bed.gz",


rule all:
    input:
        expand(rules.Rhodonite_RepeatMasker.output, sample=config["samples"].keys()),
        expand(rules.Rhodonite_trf.output, sample=config["samples"].keys()),
        expand(rules.Rhodonite_windowmasker.output, sample=config["samples"].keys()),
        expand(rules.Rhodonite_DupMasker.output, sample=config["samples"].keys()),
