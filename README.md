# Rhodonite: modular repeat masking

[![Actions Status](https://github.com/mrvollger/Rhodonite/workflows/Linting/badge.svg)](https://github.com/mrvollger/Rhodonite/actions)
[![Actions Status](https://github.com/mrvollger/Rhodonite/workflows/black/badge.svg)](https://github.com/mrvollger/Rhodonite/actions)

## Modules

- [![Actions Status](https://github.com/mrvollger/Rhodonite/workflows/RepeatMasker/badge.svg)](https://github.com/mrvollger/Rhodonite/actions)
- [![Actions Status](https://github.com/mrvollger/Rhodonite/workflows/trf/badge.svg)](https://github.com/mrvollger/Rhodonite/actions)
- [![Actions Status](https://github.com/mrvollger/Rhodonite/workflows/windowmasker/badge.svg)](https://github.com/mrvollger/Rhodonite/actions)
- [![Actions Status](https://github.com/mrvollger/Rhodonite/workflows/DupMasker/badge.svg)](https://github.com/mrvollger/Rhodonite/actions)

## Description

This is a modular Snakemake for various repeat masking steps. The `Snakefile` is under `workflow`.

## Loading in **Snakemake**

```python
configfile: "config/config.yaml"


module Rhodonite:
    snakefile:
		"https://github.com/mrvollger/Rhodonite/raw/master/workflow/Snakefile"
    config:
        config


# use all rules from Rhodonite
use rule * from Rhodonite
```

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
