#!/bin/Rscript

# non-interactive debugging
if( ! interactive() )
    options(error=dump.frames)

# AIMS:
# - adapt older code from Oli to multi-round, longitudinal settings
# TODO: discuss: we are removing individuals with missing VLs: they are very little

{
    library(data.table) |> suppressPackageStartupMessages()
    library(ggplot2)    |> suppressPackageStartupMessages()
    library(Hmisc)      |> suppressPackageStartupMessages()
    library(cmdstanr)   |> suppressPackageStartupMessages()
    library(optparse)   |> suppressPackageStartupMessages()
    library(here)       |> suppressPackageStartupMessages()
}

################
#    PATHS     #
################


# automatically finding the github directory may be complicated 
# if script is called outside from it.
self_relative_path <- 'scripts/VL_run_cmdstan.R'
if( interactive() )
{
    here::i_am(self_relative_path)
    gitdir <- here::here()
} else {
    cmd <- commandArgs()
    cmd <- cmd[cmd %like% 'file']
    gitdir <- gsub(paste0("--file=(.*)/", self_relative_path, '$'), "\\1", cmd)
}

# helpers
source(file.path(gitdir, "R/paths.R"))
source(file.path(gitdir.functions, "phsc_vl_helpers.R"))
source(file.path(gitdir.functions, "phsc_vl_cmdstan_helpers.R"))

# options (automatically sourced in R/options.R)
args_stan <- args[names(args) %like% "^iter.|chains"]
args <- args[names(args) %like% '^run|viral.load|jobname|indir|out.dir|refit|round|^only.firstparticipants$']
if(interactive()){ # testing
    args$only.firstparticipants <- TRUE 
    args$run.gp.prevl <- TRUE
    args$out.dir.exact <- "/home/andrea/HPC/ab1820/home/projects/2022/longvl//cmdstan_vl_1000/run-gp-prevl"
} 
print(args); print(args_stan)

# parallel backend (use multiple cores if local)
parallelise <- FALSE
if (parallelise) 
    source(file.path(gitdir.R, "local_cores_parallelisation.R"))

file.exists(path.hivstatusvl.r1520) |>
    all() |>
    stopifnot()


################
#     MAIN     #
################

# check study round exists
stopifnot(all(args$round %in% 16:19))

# set viral load thresholds
VL_DETECTABLE <- args$vl.detectable
VIREMIC_VIRAL_LOAD <- args$viremic.viral.load

# specify and create output directories as func(vl, jobname)
if(! is.na(args$out.dir.exact) ){
    vl.out.dir <- args$out.dir.exact
    dir.create(vl.out.dir, recursive = TRUE)
}else{
    stopifnot(dir.exists(args$out.dir.prefix))
    out.dir <- file.path(args$out.dir.prefix)
    suffix <- make.suffix(args, cmdstan=TRUE)
    vl.out.dir <- file.path( out.dir, suffix)
}
cat('vl.out.dir specified as:\n ', vl.out.dir, '\n')
dir.create(vl.out.dir, showWarnings = FALSE)
stopifnot(dir.exists(vl.out.dir))

# get data
dall <- get.dall(path = path.hivstatusvl.r1520, only_firstpart = args$only.firstparticipants)
dall <- subset(dall, ROUND %in% args$round)


# Estimate HIV prevalence
# ________________________

if (args$run.gp.prevl) {
    vl.prevalence.by.gender.loc.age.gp.cmdstan(dall, refit = args$refit)
}

# Estimate mean viral load
# ________________________

if (args$run.icar.mean.vl) {
    vl.meanviralload.by.gender.loc.age.icar.cmdstan(dall, refit = args$refit)
}

# Estimate suppressed pop
# _______________________

if (args$run.gp.supp.hiv) { # Among HIV positive
    vl.suppofinfected.by.gender.loc.age.gp.cmdstan(dall, refit = args$refit)
}


if (args$run.gp.supp.pop) { # Among Entire population
    vl.suppofpop.by.gender.loc.age.gp.cmdstan(dall, refit = args$refit)
}
