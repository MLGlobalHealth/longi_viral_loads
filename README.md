## Introduction
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![medRxiv link](https://img.shields.io/badge/medRxiv-link%20to%20paper-blue)](TODO)
[![Zeondo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10955672.svg)](https://doi.org/10.5281/zenodo.10955672)


**Hello!** This repository contains the code for the analyses presented in the paper *Age and gender profiles of HIV infection burden and viraemia: novel metrics for HIV epidemic control in African populations with high antiretroviral therapy coverage* by A Brizzi et al.

## Table of contents
 
- [Quick Start](#rocket-quick-start)
    - [System Requirements](#system-requirements)
    - [Installation](#installation)
    - [Reproducing our Analyses](#reproducing-our-analyses)
- [License](#pagefacingup-license)
- [Warranty](#shield-warranty)
- [Citation](#books-citation)
- [Acknowledgments](#acknowledgments)
- [Funding](#funding)

## :rocket: Quick Start

### System Requirements

- UNIX or MacOs, the code was developed on Ubuntu 22.04 LTS.
- [R](https://www.r-project.org/) version `≥4.1.2`.
- [cmdstan](https://mc-stan.org/users/interfaces/cmdstan) version `2.27.0`.

We additionally suggests [conda or miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) to reproduce our package environment.

### Installation

You can follow the steps below to clone our codebase locally and rerun the analyses. 

1. Clone the github repository in a directory of your choice.
2. Install the dependencies by running `sh dependencies_cmdstan.sh`.
3. Download the data necessary to run the analyses from [Zenodo](https://doi.org/10.5281/zenodo.10955672).
4. Specify the path to downloaded data directory in `R/paths.R` at line 22: \
    `indir.zenodo <- "PATH-TO-DOWNLOADED-ZENODO-DATA"`
 
>  **Note** A common error with new installations of `cmdstanr` is the inability to compile stan files due to a `need to set TBB_CXX_TYPE` (see [here](https://bytemeta.vip/repo/stan-dev/cmdstanpy/issues/374?page=1) for a solution).

### Reproducing our Analyses

You are now in the position to run our code, reproduce our analyses, or explore new possibilities.
We suggest to save the results of your runs in the `results` subdirectory of our zenodo repository.
As such, we suggest to define the following environment variable on the command line:

```{sh}
ZENODO_DIR="PATH-TO-DOWNLOADED-ZENODO-DATA"
```

To run the code with the same settings as in our paper, we then run from the root directory of the repository the [`run_stan`](run_stan) executable. 

```{sh}
./run_stan OUTDIR="$ZENODO_DIR/results" ROUND=19 MODELS="run-gp-prevl" LOCAL=TRUE REFIT=TRUE
```

Note that: 
1. The above can also be run as a bash script through `bash run_stan ...`.
2. alternative output directories can be used, as in `OUTDIR="path/to/other/directory"`.
3. stan argumebts can be updated in the [`stan/binomial_gp_config.yml`](stan/binomial_gp_config.yml) file.
 
Once all jobs are run, it is possible to run our postprocessing code, involving reproducibility of figures, tables and an `html` report summarising the runs.

```{sh}
./run_postprocessing OUTDIR="$ZENODO_DIR/results" LOCAL=TRUE
```

All the `R` scripts called by the two executable can be found in the [`scripts`](scripts).

> **Note on flags.** The two executables [`run_stan`](run_stan) and [`run_postprocessing`](run_postprocessing) take a number of flags whose documentation can be read through the `-h` flag (eg: `./run_stan -h`). The most important are:
> * `OUTDIR`: indicates the directory where we want to save the results of our runs.
> * `LOCAL`: indicates we want to run jobs locally, instead of submitting them to a job scheduler with `qsub`.
> * `MODELS` indicates which models we want to run, and can be one, or a combination of, the following: `run-gp-prevl`, `run-gp-supp-pop`, `run-gp-supp-hiv`.
> * `ROUND` flag indicates the survey round for which we want to run our analyses 
> * `REFIT` flag indicates whether we want to refit the model and overwrite them in case a previous version already exists.
> 
> Every flag should be specified as `FLAG=value` or `FLAG="value with spaces"`, leaving no white spaces before and after the `=` sign.

## :page_facing_up: License 

The code in this repository is licensed under [CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg).

## :shield: Warranty 

Imperial makes no representation or warranty about the accuracy or completeness of the data nor that the results will not constitute in infringement of third-party rights. Imperial accepts no liability or responsibility for any use which may be made of any results, for the results, nor for any reliance which may be placed on any such work or results.

## :books: Citation

Please cite this work as *TODO* (leave a citation.bib in the repo?)

## Acknowledgments

We thank all contributors, program staff and participants to the Rakai Community Cohort Study; the [Rakai Health Sciences Program](https://www.rhsp.org/index.php).

We also extend our gratitude to the [Imperial College Research Computing Service](https://doi.org/10.14469/hpc/2232) for providing the computational resources to perform this study. 

## Funding

> This study was supported by the National Institute of Allergy and Infectious Diseases [U01AI075115 to RHG, R01AI087409 to RHG, U01AI100031 to RHG, ZIAAI001040 to TCQ]; the National Institute of Mental Health [F31MH095649 to Dr Jennifer Wagman, R01MH099733 to Ned Sacktor and MJW, R01MH107275 to LWC]; the Division of Intramural Research of the National Institute for Allergy and Infectious Diseases [TCQ, OL, SJR], NIAID [K01AA024068 to Dr Jennifer Wagman]; the Johns Hopkins University Center for AIDS Research [2P30AI094189 to Dr Richard Chaisson]; the U.S. President’s Emergency Plan for AIDS Relief (PEPFAR) through the Centers for Disease Control and Prevention [NU2GGH000817 to RHSP]; the Engineering and Physical Sciences Research Council through the EPSRC Centre for Doctoral Training in Modern Statistics and Statistical Machine Learning at Imperial and Oxford [EP/S023151/1 to Prof Axel Gandy]. The funders had no role in study design, data collection and analysis, decision to publish or preparation of the manuscript. The findings and conclusions in this report are those of the author(s) and do not necessarily represent the official position of the Centers for Disease Control and Prevention. 
