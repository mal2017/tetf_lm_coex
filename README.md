# TETF_LM_COEX: processing of DGRP expression data into TF/gene coexpression relationships for the Ellison lab TETF project

## dependencies

- `snakemake` (tested with v7.21.0 on Ubuntu 22 LTS and 7.18.2 on CentOS 7)
- `mamba` (tested with version 1.0.0 on Ubuntu 22 LTS and 0.9.1 on CentOS 7). `conda` alone should work just fine.
- `singularity` (tested with v3.1.0-1 on CentOS 7) or `apptainer` (tested with v1.1.3 on Ubuntu 22 LTS)

## usage

```
snakemake --use-conda --use-singularity --cores <threads>
```

A test run when upstream pipelines have also been run in `test` mode, a smaller test run can be performed by adding `--config SALMON_REPLICATES=1 RUN_TYPE=test`

## Description

This workflow fits linear models relating expression of each TE to each protein codig gene. It depends on several other `tetf_*` repositories developed for the project and is intended to be run as a submodule/subworkflow of the parent repo: `tetf`.

Specifically, it assumes that directory `tetf_expression_quant/` exists with the same parent dir as this working dir.
