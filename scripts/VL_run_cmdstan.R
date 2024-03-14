#!/bin/Rscript

self_relative_path <- "scripts/VL_run_cmdstan.R"

########################
cat("\nStart of:", self_relative_path, "\n")
########################

# non-interactive debugging
if( ! interactive() )
    options(error=dump.frames)

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
    dall[, mean(!is.na(HIV_STATUS))]
    dall[HIV_STATUS == TRUE, mean(!is.na(VL_COPIES))]
    dall <- subset(dall, 
        ROUND %in% args$round & ! is.na(HIV_STATUS)
    )
    dall[, ARVMED := NA]

    tmp <- .preprocess.ds.oli(dall, rm.na.vl=FALSE)
    # tmp[ HIV_STATUS == 1 & is.na(VL_COPIES)]

    cols <- c("N", "HIV_N", "VLNS_N", "VLNA_N")
    vla_all <- .preprocess.make.vla(tmp, select = cols )
    # vla_all[ VLNA_N != 0]

    vla_ftp <- .preprocess.make.vla(tmp[FIRST_PARTICIPATION == 1], select = cols )
    vla <- rbind(vla_all[, PTYPE := "all"], vla_ftp[, PTYPE := "ftp"])
    setkeyv(vla, c("ROUND", "LOC_LABEL", "SEX_LABEL", "PTYPE", "AGE_LABEL"))

    if( ! file.exists( path.aggregated.nums.denoms.r1619) ){
        sprintf( "Saving file: %s\n", path.aggregated.nums.denoms.r1619) |> cat()
        fwrite(x=vla, file=path.aggregated.nums.denoms.r1619 )
    }

}else{
    # load aggregated data...
    vla <- fread( path.aggregated.nums.denoms.r1619 )
}

# Estimate HIV seroprevalence
# ________________________

if(args$shared.hyper){
    dpartrates <- readRDS(path.participation.rates) |>
        subset(select = c("ROUND", "FC", "SEX", "AGEYRS", "PARTRATE_SMOOTH.25")) |>
        setnames(c("FC", "PARTRATE_SMOOTH.25"), c("LOC", "PARTRATE"))
    source(file.path(gitdir.functions, "phsc_vl_cmdstan_helpers_sharedhyper.R"))
}

if (args$run.gp.prevl) {
    # Number of HIV positive out of individuals with known HIV status
    if(args$shared.hyper){
        stopifnot("args$only.firstparticipants must not be TRUE"= ! args$only.firstparticipants)
        vl.prevalence.by.gender.loc.age.gp.cmdstan.hyper(vla, refit = args$refit, alpha_hyper = args$stan.alpha)
    }else{
        vla <- if(args$only.firstparticipants){vla[PTYPE == "ftp"]}else{vla[PTYPE=="all"]}
        vl.prevalence.by.gender.loc.age.gp.cmdstan(vla, refit = args$refit, alpha_hyper = args$stan.alpha)
    }
}

# Estimate suppressed pop
# _______________________

if (args$run.gp.supp.hiv) { # Among HIV positive
    # Number of Virally supppressed out of PLHIV with known VL
    if(args$shared.hyper){
        stopifnot("args$only.firstparticipants must not be TRUE"= ! args$only.firstparticipants)
        vl.suppofinfected.by.gender.loc.age.gp.cmdstan.hyper(vla, refit = args$refit, alpha_hyper = args$stan.alpha)
    }else{
        vla <- if(args$only.firstparticipants){vla[PTYPE == "ftp"]}else{vla[PTYPE=="all"]}
        vl.suppofinfected.by.gender.loc.age.gp.cmdstan(vla, refit = args$refit, alpha_hyper = args$stan.alpha)
    }
}


if (args$run.gp.supp.pop) { # Among Entire population
    # Number of virally unsuppressed out of entire pop
    if(args$shared.hyper){
        stopifnot("args$only.firstparticipants must not be TRUE"= ! args$only.firstparticipants)
        vl.suppofpop.by.gender.loc.age.gp.cmdstan.hyper(vla, refit = args$refit, alpha_hyper = args$stan.alpha)
    }else{
        vla <- if(args$only.firstparticipants){vla[PTYPE == "ftp"]}else{vla[PTYPE=="all"]}
        vl.suppofpop.by.gender.loc.age.gp.cmdstan(vla, refit = args$refit, alpha_hyper = args$stan.alpha)
    }
}
