# Rhodonite: modular repeat masking

[![DOI](https://zenodo.org/badge/388233694.svg)](https://zenodo.org/badge/latestdoi/388233694)
[![Actions Status](https://github.com/mrvollger/Rhodonite/workflows/Linting/badge.svg)](https://github.com/mrvollger/Rhodonite/actions)
[![Actions Status](https://github.com/mrvollger/Rhodonite/workflows/black/badge.svg)](https://github.com/mrvollger/Rhodonite/actions)

## Modules

- [![Actions Status](https://github.com/mrvollger/Rhodonite/workflows/RepeatMasker/badge.svg)](https://github.com/mrvollger/Rhodonite/actions)
- [![Actions Status](https://github.com/mrvollger/Rhodonite/workflows/trf/badge.svg)](https://github.com/mrvollger/Rhodonite/actions)
- [![Actions Status](https://github.com/mrvollger/Rhodonite/workflows/windowmasker/badge.svg)](https://github.com/mrvollger/Rhodonite/actions)
- [![Actions Status](https://github.com/mrvollger/Rhodonite/workflows/DupMasker/badge.svg)](https://github.com/mrvollger/Rhodonite/actions)
- [![Actions Status](https://github.com/mrvollger/Rhodonite/workflows/dna-brnn/badge.svg)](https://github.com/mrvollger/dna-brnn/actions)

## Description

This is a modular Snakemake workflow for various repeat masking tools.

## Loading in **Snakemake**

```python
from snakemake.utils import min_version

min_version("6.0")

# define your configuration via python
# or define a yaml e.g.
# configfile: "config/config.yaml"
config = {
    "samples": {
        "test": ".test/test.fa.gz",
        "test2": ".test/test.fa.gz",
    }
}


module Rhodonite:
    snakefile:
        "https://github.com/mrvollger/Rhodonite/raw/v0.3-alpha/workflow/Snakefile"
    config:
        config


# import the rules from Rhodonite
use rule * from Rhodonite as Rhodonite_*


# use RepeatMasker rules from Rhodonite
use rule RepeatMasker from Rhodonite as Rhodonite_RepeatMasker with:
    output:
        # You rename the output to anything you want
        # but maintain the order and keep "{sample}" in the name and the (.gz).
        # Every rule will have the standard output of the program
        # and a bed output that is sorted in the order of the ref.
        out="{sample}.rm.out",
        bed="{sample}.rm.bed.gz",


# use trf rules from Rhodonite
use rule trf from Rhodonite as Rhodonite_trf with:
    output:
        dat="{sample}.trf.dat",
        bed="{sample}.trf.bed.gz",


# use DupMasker rules from Rhodonite
use rule DupMasker from Rhodonite as Rhodonite_DupMasker with:
    output:
        extra="{sample}.duplicons.extra",
        bed="{sample}.duplicons.bed.gz",


rule all:
    input:
        expand(rules.Rhodonite_RepeatMasker.output, sample=config["samples"].keys()),
        expand(rules.Rhodonite_trf.output, sample=config["samples"].keys()),
        expand(rules.Rhodonite_DupMasker.output, sample=config["samples"].keys()),

```

For a **complete** documented example see `examples/example.smk`.
