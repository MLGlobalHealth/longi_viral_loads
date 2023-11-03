[#](#) Introduction
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![medRxiv link](https://img.shields.io/badge/medRxiv-link%20to%20paper-blue)](TODO)


**Hello!** This repository contains the code for the analyses presented in the paper *TITLE* by A Brizzi et al.

## Table of contents


## License 

The code in this repository is licensed under [CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg).

## Warranty 

Imperial makes no representation or warranty about the accuracy or completeness of the data nor that the results will not constitute in infringement of third-party rights. Imperial accepts no liability or responsibility for any use which may be made of any results, for the results, nor for any reliance which may be placed on any such work or results.

## Citation:

Please cite this work as *TODO* (leave a citation.bib in the repo?)

## Acknowledgments

We thank all contributors, program staff and participants to the Rakai Community Cohort Study; all members of the PANGEA-HIV consortium, the [Rakai Health Sciences Program](https://www.rhsp.org/index.php).

We also extend our gratitude to the [Imperial College Research Computing Service](https://doi.org/10.14469/hpc/2232) for providing the computational resources to perform this study. 

## Funding

TOCOPY FROM PAPER.


## Quick Start

### System Requirements

- UNIX or MacOs, the code was developed on Ubuntu 22.04 LTS.
- [R](https://www.r-project.org/) version `≥4.1.2`.
- [cmdstan](https://mc-stan.org/users/interfaces/cmdstan) version TODO....

We additionally [conda]() or [miniconda]() to reproduce our package environment.

### Installation

You can follow the steps below to clone our codebase locally and rerun the analyses. 

1. Clone the github repository in a directory of your choice.
2. Install the dependencies by running `sh dependencies_cmdstan.sh`.
 
> ** Note ** A common error with new installations of `cmdstanr` is the inability to compile stan files due to a `need to set TBB_CXX_TYPE` (see [here](https://bytemeta.vip/repo/stan-dev/cmdstanpy/issues/374?page=1) for a solution).

### Reproducing our Analyses

You are now in the position to run our code, reproduce our analyses, or explore new possibilities.
To run the code with the same settings as in our paper, run from the root directory of the repository:

```{sh}
# Run ./run_stan -h for more help on the available flags
./run_stan OUTDIR="path/to/desired/output/directory"
```

[ MAYBE ADD FLAGS `ROUND` `MODEL` AND `LOCAL` TO NOT RUN TOO MANY THINGS AT ONCE??? ]
The above prepares a different bash job for every combination of the three prevalence measures of interest and each of the four survey rounds. The resulting `.pbs` files are saved in the `OUTDIR` and can be either submitted through a `PBS` job scheduler, or they can be run as bash scripts individually.
Note that running the script redirects you to the specified `OUTDIR`.

Once all jobs are run, it is possible to run our postprocessing code, involving reproducibility of figures, tables and an `html` report summarising the runs.

```{sh}
# Run ./run_postprocessing -h for more help on the available flags
./run_postprocessing OUTDIR="path/to/desired/output/directory"
```

### Data.

Data used exclusively for this analysis is stored on the HPC, at `$DEEPDATA/RCCS_R15_R20`. 
The data consists of: 

- `Quest_R015_R020_VOIs_May062022.csv`, containing Viral Loads measurements for HIV positive participants in rounds 15 to 20. As of May 2022, round 20 still hasn't completed, so the data are not definitive. Also, some data are missing for round 15.
Concerning ARV: JS: "the main interest was to capture the round-specific interview dates where arvmed/cuarvmed==1...Ideally if you are to come up with a binary variable still those 2s and 8s will be 0s": the problem is that I cannot distinguish those participants that reporterted not-adering to treatment, and those participants who didn't disclose that information.

- DESCRIBE: `new questionnaire data`

- The `RCCS_census_eligible_individuals_221209.csv` I have in gitdir.data was produced with Melodies `misc/get_census_eligible_count.R` but additionaly including round 19. NOTE: should find where this is

## Code 

The code to reproduce the analyses is split into multiple parts:
1. Data cleaning and preparation
2. Running statistical models
3. Analysis of the results and postprocessing.


### Data cleaning and preparation

Data preprocessing is done in scripts within the `src/preprocess` directory. 

### Running statistical models

Statistical models are run through the main script `src/run_stan.R`.
Reproducing the model analyses consists in multiple steps.
First, we need to run the models on both the set of all participants, or the set of first-time participants. 
These jobs can be run in parallel, and can be simultaneously be submitted to the HPC:
```{bash}
./run_stan.sh 1000
./run_stan.sh 1000 firstpart 
```


### Postprocessing

To estimate prevalences at the census-eligible population level, we assume that prevalences are homogeneous out-of-study and within first-time participants.
As such, we need both sets of results to obtain the results descibed in the paper.


