#
# calls all output of the pipeline.
#
include: "../Snakefile"


rule test:
    input:
        fasta=gather.fasta(
            rules.split_fasta.output.fasta, sample=config["samples"].keys()
        ),
        trf=expand(rules.trf.output, sample=config["samples"].keys()),
        RepeatMasker=expand(
            rules.RepeatMasker.output, sample=config["samples"].keys()
        ),
        DupMasker=expand(rules.DupMasker.output, sample=config["samples"].keys()),
        windowmasker=expand(
            rules.windowmasker.output, sample=config["samples"].keys()
        ),
