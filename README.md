# Introduction

This repository contains the code used to analyse Viral Load data obtained from Rakai, with the objective of answering some questions posed in the [LONGVIEW](TODO:LINK) grant.

# Data.

Data used exclusively for this analysis is stored on the HPC, at `$DEEPDATA/RCCS_R15_R20`. 
The data consists of: 

- `Quest_R015_R020_VOIs_May062022.csv`, containing Viral Loads measurements for HIV positive participants in rounds 15 to 20. As of May 2022, round 20 still hasn't completed, so the data are not definitive. Also, some data are missing for round 15.
Concerning ARV: JS: "the main interest was to capture the round-specific interview dates where arvmed/cuarvmed==1...Ideally if you are to come up with a binary variable still those 2s and 8s will be 0s": the problem is that I cannot distinguish those participants that reporterted not-adering to treatment, and those participants who didn't disclose that information.

- DESCRIBE: `new questionnaire data`

- The `RCCS_census_eligible_individuals_221209.csv` I have in gitdir.data was produced with Melodies `misc/get_census_eligible_count.R` but additionaly including round 19. NOTE: should find where this is

# Ideas:

Community level viral suppression is, for the most part, determined by HIV positive individuals aware of their status.

These are exactly the people who do not benefit from extra testing: they benefit from accessing treatment.

Thus, at the participant level, I would not be shocked to see that suppressed individuals are as likely not to test as suppressed ones.

Based on this, I would maybe look at the relationship between community level suppression vs « testing rates » among negatives and new positives.
