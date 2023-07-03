# Introduction

This repository contains the code used to analyse Viral Load data obtained from Rakai, with the objective of answering some questions posed in the [LONGVIEW](TODO:LINK) grant.

# Installation:

There are 2 installation scripts ...

Concerning cmdstanr, you may get an error such as: "Need to set TBB_CXX_TYPE" (see [here](https://bytemeta.vip/repo/stan-dev/cmdstanpy/issues/374?page=1) )
In this case, it is necessary to locat the path to the cmdstanr installed repository and add the following lines to the `make/local` file.
```{bash make/local}
CXX=g++
TBB_CXX_TYPE=gcc
```

# Data.

Data used exclusively for this analysis is stored on the HPC, at `$DEEPDATA/RCCS_R15_R20`. 
The data consists of: 

- `Quest_R015_R020_VOIs_May062022.csv`, containing Viral Loads measurements for HIV positive participants in rounds 15 to 20. As of May 2022, round 20 still hasn't completed, so the data are not definitive. Also, some data are missing for round 15.
Concerning ARV: JS: "the main interest was to capture the round-specific interview dates where arvmed/cuarvmed==1...Ideally if you are to come up with a binary variable still those 2s and 8s will be 0s": the problem is that I cannot distinguish those participants that reporterted not-adering to treatment, and those participants who didn't disclose that information.

- DESCRIBE: `new questionnaire data`

- The `RCCS_census_eligible_individuals_221209.csv` I have in gitdir.data was produced with Melodies `misc/get_census_eligible_count.R` but additionaly including round 19. NOTE: should find where this is

# Code 

The code to reproduce the analyses is split into multiple parts:
1. Data cleaning and preparation
2. Running statistical models
3. Analysis of the results and postprocessing.


## Data cleaning and preparation

Data preprocessing is done in scripts within the `src/preprocess` directory. 

## Running statistical models

Statistical models are run through the main script `src/run_stan.R`.
Reproducing the model analyses consists in multiple steps.
First, we need to run the models on both the set of all participants, or the set of first-time participants. 
These jobs can be run in parallel, and can be simultaneously be submitted to the HPC:
```{bash}
./run_stan.sh 1000
./run_stan.sh 1000 firstpart 
```


## Postprocessing

To estimate prevalences at the census-eligible population level, we assume that prevalences are homogeneous out-of-study and within first-time participants.
As such, we need both sets of results to obtain the results descibed in the paper.


