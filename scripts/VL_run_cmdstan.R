#!/bin/Rscript

self_relative_path <- "scripts/VL_run_cmdstan.R"

########################
cat("\nStart of:", self_relative_path, "\n")
########################

# non-interactive debugging
if( ! interactive() )
    options(error=dump.frames)

# AIMS:
# - adapt older code from Oli to multi-round, longitudinal settings
# TODO: discuss: we are removing individuals with missing VLs: they are very little

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(Hmisc)
    library(cmdstanr)
    library(optparse)
    library(here)
})

################
#    PATHS     #
################

# automatically finding the github directory may be complicated 
# if script is called outside from it.
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
source(file.path(gitdir.functions, "preprocessing_tables.R"))
source(file.path(gitdir.functions, "phsc_vl_cmdstan_helpers.R"))

# options (automatically sourced in R/options.R)
args_stan <- args[names(args) %like% "^iter.|chains"]
args <- args[names(args) %like% '^run|viral.load|jobname|indir|out.dir|refit|round|^only.firstparticipants$|^stan.alpha$|^shared.hyper$|confidential']
if(interactive()){ # testing
    args$only.firstparticipants <- FALSE 
    args$run.gp.supp.pop <- TRUE
    args$round <- 19
    args$stan.alpha <- 1.00
    args$shared.hyper <- TRUE
} 
stopifnot(!(args$shared.hyper & args$only.firstparticipants))
print(args); print(args_stan)

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
VL_DETECTABLE <- args$vl.detectable.viral.load
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
if( args$confidential ){

    dall <- get.dall(path = path.hivstatusvl.r1520, only_firstpart = args$only.firstparticipants)
    dall <- subset(dall, ROUND %in% args$round)
    # if(local()){
    #     tmp <- make.table.firstparticipant.NPhiv.NPunsupp(dall)
    #     write.to.googlesheets(tmp, sheet='SuppTable1')
    # }

    tmp <- .preprocess.ds.oli(dall)
    cols <- c("N", "HIV_N", "VLNS_N", "ARV_N")
    vla_all <- .preprocess.make.vla(tmp, select = cols )
    vla_ftp <- .preprocess.make.vla(tmp[FIRST_PARTICIPATION == 1], select = cols )
    vla <- rbind(vla_all[, PTYPE := "all"], vla_ftp[, PTYPE := "ftp"])
    if( ! file.exists( path.aggregated.nums.denoms.r1619) ){
        sprintf( "Saving file: %s\n", path.aggregated.nums.denoms.r1619) |> cat()
        fwrite(x=vla, file=path.aggregated.nums.denoms.r1619 )
    }

}else{
    # load aggregated data...
    vla <- fread( path.aggregated.nums.denoms.r1619 )
}

# Estimate HIV prevalence
# ________________________

if(args$shared.hyper){
    dpartrates <- readRDS(path.participation.rates) |>
        subset(select = c("ROUND", "FC", "SEX", "AGEYRS", "PARTRATE_SMOOTH.25")) |>
        setnames(c("FC", "PARTRATE_SMOOTH.25"), c("LOC", "PARTRATE"))
    source(file.path(gitdir.functions, "phsc_vl_cmdstan_helpers_sharedhyper.R"))
}

if (args$run.gp.prevl) {
    if(args$shared.hyper){
        vl.prevalence.by.gender.loc.age.gp.cmdstan.hyper(vla, refit = args$refit, alpha_hyper = args$stan.alpha)
    }else{
        vl.prevalence.by.gender.loc.age.gp.cmdstan(vla, refit = args$refit, alpha_hyper = args$stan.alpha)
    }
}

# Estimate suppressed pop
# _______________________

if (args$run.gp.supp.hiv) { # Among HIV positive
    if(args$shared.hyper){
        vl.suppofinfected.by.gender.loc.age.gp.cmdstan.hyper(vla, refit = args$refit, alpha_hyper = args$stan.alpha)
    }else{
        vl.suppofinfected.by.gender.loc.age.gp.cmdstan(vla, refit = args$refit, alpha_hyper = args$stan.alpha)
    }
}


if (args$run.gp.supp.pop) { # Among Entire population
    if(args$shared.hyper){
        vl.suppofpop.by.gender.loc.age.gp.cmdstan.hyper(vla, refit = args$refit, alpha_hyper = args$stan.alpha)
    }else{
        vl.suppofpop.by.gender.loc.age.gp.cmdstan(vla, refit = args$refit, alpha_hyper = args$stan.alpha)
    }
}
