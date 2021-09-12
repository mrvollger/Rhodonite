# Rhodonite: modular repeat masking
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5501705.svg)](https://doi.org/10.5281/zenodo.5501705)
[![Actions Status](https://github.com/mrvollger/Rhodonite/workflows/Linting/badge.svg)](https://github.com/mrvollger/Rhodonite/actions)
[![Actions Status](https://github.com/mrvollger/Rhodonite/workflows/black/badge.svg)](https://github.com/mrvollger/Rhodonite/actions)

## Modules

- [![Actions Status](https://github.com/mrvollger/Rhodonite/workflows/RepeatMasker/badge.svg)](https://github.com/mrvollger/Rhodonite/actions)
- [![Actions Status](https://github.com/mrvollger/Rhodonite/workflows/trf/badge.svg)](https://github.com/mrvollger/Rhodonite/actions)
- [![Actions Status](https://github.com/mrvollger/Rhodonite/workflows/windowmasker/badge.svg)](https://github.com/mrvollger/Rhodonite/actions)
- [![Actions Status](https://github.com/mrvollger/Rhodonite/workflows/DupMasker/badge.svg)](https://github.com/mrvollger/Rhodonite/actions)

## Description

This is a modular Snakemake workflow for various repeat masking tools.

## Loading in **Snakemake**

```python
from snakemake.utils import min_version

min_version("6.0")

config = {
    "samples": {
        "test": ".test/test.fa.gz",
    }
}


module Rhodonite:
    snakefile:
        "https://github.com/mrvollger/Rhodonite/raw/master/workflow/Snakefile"
    config:
        config


# import the rules from Rhodonite
use rule * from Rhodonite as Rhodonite_*


# use RepeatMasker rules from Rhodonite
use rule RepeatMasker from Rhodonite as Rhodonite_RepeatMasker with:
    output:
        out="{sample}.rm.out",
        bed="{sample}.rm.bed.gz",
```

For a **complete** documented example see `examples/example.smk`.

## Loading in **Python**

Requires: [Snakedeploy](https://snakedeploy.readthedocs.io/en/latest/index.html)

```python
from snakedeploy.deploy import deploy

deploy("https://github.com/mrvollger/Rhodonite",
		dest_path=".",
		name="Rhodonite",
		tag="master",
		force=True
	  )
```
