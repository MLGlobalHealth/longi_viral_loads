# Introduction

This repository contains the code used to analyse Viral Load data obtained from Rakai, with the objective of answering some questions posed in the [LONGVIEW](TODO:LINK) grant.


# Data.

Data used exclusively for this analysis is stored on the HPC, at `$DEEPDATA/RCCS_R15_R20`. 
The data consists of: 

- `Quest_R015_R020_VOIs_May062022.csv`, containing Viral Loads measurements for HIV positive participants in rounds 15 to 20. As of May 2022, round 20 still hasn't completed, so the data are not definitive. Also, some data are missing for round 15.

# Code

# Objective
Look at `concept_analysis.md` in `Documents/HIV_rccs`

- Complex trajectories f(treatment, prevention) cannot be represented as averages.
- $\gamma\_i$ used for VL trajectory of $i$th participant

- Estimate Group trajectories that are exhaustively representative of individual traj.
- Probabilities $\pi\_{ij}$ that ind i falls in group j.

- "Our prior work indicates individual trajectories will be categorized as
durably suppressed, newly HIV suppressed, intermittently viremic and
persistently viremic, though the Bayesian latent class analysis may
reveal further group trajectories.""


# Whatever else

