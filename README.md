## Introduction
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![medRxiv link](https://img.shields.io/badge/medRxiv-link%20to%20paper-blue)](TODO)


**Hello!** This repository contains the code for the analyses presented in the paper *TITLE* by A Brizzi et al.
(I would put the citation here)

## Table of contents
 
- [License](#license)
- [Quick Start](#quick-start)
    - [System Requirements](#system-requirements)
    - [Installation](#installation)
    - [Reproducing our Analyses](#reproducing-our-analyses)
- [Warranty](#warranty)
- [Citation](#citation)
- [Acknowledgments](#acknowledgments)
- [Funding](#funding)


## License 

The code in this repository is licensed under [CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg).


## Quick Start

### System Requirements

- UNIX or MacOs, the code was developed on Ubuntu 22.04 LTS.
- [R](https://www.r-project.org/) version `≥4.1.2`.
- [cmdstan](https://mc-stan.org/users/interfaces/cmdstan) version `2.27.0`.

We additionally [conda]() or [miniconda]() to reproduce our package environment.

### Installation

You can follow the steps below to clone our codebase locally and rerun the analyses. 

1. Clone the github repository in a directory of your choice.
2. Install the dependencies by running `sh dependencies_cmdstan.sh`.
3. Download the data necessary to run the analyses from TODO and specify the path to that directory in TODO . 
 
>  **Note** A common error with new installations of `cmdstanr` is the inability to compile stan files due to a `need to set TBB_CXX_TYPE` (see [here](https://bytemeta.vip/repo/stan-dev/cmdstanpy/issues/374?page=1) for a solution).

### Reproducing our Analyses

You are now in the position to run our code, reproduce our analyses, or explore new possibilities.
To run the code with the same settings as in our paper, run from the root directory of the repository:

```{sh}
# Run ./run_stan -h for more help on the available flags
./run_stan OUTDIR="path/to/desired/output/directory" ROUND=19 MODELS="run-gp-prevl" LOCAL=TRUE
```

> **Note on flags**
> * `LOCAL`: indicates we want to run jobs locally, instead of submitting them to a job scheduler with `qsub`.
> * `MODELS` indicates which models we want to run, and can be one, or a combination of, the following: `run-gp-prevl`, `run-gp-supp-pop`, `run-gp-supp-hiv`.
> * `ROUND` flag indicates the survey round for which we want to run our analyses

Once all jobs are run, it is possible to run our postprocessing code, involving reproducibility of figures, tables and an `html` report summarising the runs.

```{sh}
# Run ./run_postprocessing -h for more help on the available flags
./run_postprocessing OUTDIR="path/to/desired/output/directory" LOCAL=TRUE
```

## Warranty 

Imperial makes no representation or warranty about the accuracy or completeness of the data nor that the results will not constitute in infringement of third-party rights. Imperial accepts no liability or responsibility for any use which may be made of any results, for the results, nor for any reliance which may be placed on any such work or results.

## Citation:

Please cite this work as *TODO* (leave a citation.bib in the repo?)

## Acknowledgments

We thank all contributors, program staff and participants to the Rakai Community Cohort Study; all members of the PANGEA-HIV consortium, the [Rakai Health Sciences Program](https://www.rhsp.org/index.php).

We also extend our gratitude to the [Imperial College Research Computing Service](https://doi.org/10.14469/hpc/2232) for providing the computational resources to perform this study. 

## Funding

TOCOPY FROM PAPER.
