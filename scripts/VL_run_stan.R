#!/bin/Rscript

# AIMS:
# - adapt older code from Oli to multi-round, longitudinal settings
# TODO: discuss: we are removing individuals with missing VLs: they are very little

{
    library(data.table) |> suppressPackageStartupMessages()
    library(ggplot2)    |> suppressPackageStartupMessages()
    library(Hmisc)      |> suppressPackageStartupMessages()
    library(rstan)      |> suppressPackageStartupMessages()
    library(optparse)   |> suppressPackageStartupMessages()
    library(here)       |> suppressPackageStartupMessages()
}

################
#    PATHS     #
################


# automatically finding the github directory may be complicated 
# if script is called outside from it.
self_relative_path <- 'scripts/VL_run_stan.R'
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

# options (automatically sourced in R/options.R)
args <- args[names(args) %like% '^run|viral.load|jobname|indir|out.dir|refit|round|^only.firstparticipants']
if(interactive()){ # testing
    args$only.firstparticipants <- TRUE 
    args$run.gp.prevl <- TRUE
} 
print(args)

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
    suffix <- make.suffix(args)
    vl.out.dir <- file.path( out.dir, suffix)
}
cat('vl.out.dir specified as:\n ', vl.out.dir, '\n')
dir.create(vl.out.dir, showWarnings = FALSE)
stopifnot(dir.exists(vl.out.dir))

# get data
dall <- get.dall(path = path.hivstatusvl.r1520, only_firstpart = args$only.firstparticipants)

if (0) { # Info for introduction to results

    .mean2 <- function(x) paste0(round(100 * mean(x), 2), "%")

    # proportions of viraemic measurements
    dall[VL_COPIES > VIREMIC_VIRAL_LOAD, .mean2(SEX == "F"), ]
    dall[VL_COPIES > VIREMIC_VIRAL_LOAD, .mean2(SEX == "F"), by = "ROUND"]
    dall[VL_COPIES > VIREMIC_VIRAL_LOAD, .mean2(FC == "inland"), by = "ROUND"]

    # ARVs
    dall[HIV_AND_VL == 1, .mean2(is.na(ARVMED)), by = "ROUND"]
    cat("mean log10 VL across people not reporting ARVMED\n")
    dall[HIV_AND_VL == 1 & is.na(ARVMED), .(VLmean = mean(log(VL_COPIES + 1, 10))), by = "ROUND"]
    cat("mean log10 VL across people reporting ARVMED\n")
    dall[HIV_AND_VL == 1 & !is.na(ARVMED), .(VLmean = mean(log(VL_COPIES + 1, 10))), by = "ROUND"]
    #
}

# Make some plots? 
# ________________



# community analysis
# __________________

if (args$run.comm.analysis) {
    source(file.path(gitdir.scripts, 'community_analysis.R'))
}

# Estimate HIV prevalence
# ________________________

if (args$run.gp.prevl) {
    vl.prevalence.by.gender.loc.age.gp(dall, refit = args$refit)
}

# Estimate mean viral load
# ________________________

if (args$run.icar.mean.vl) {
    vl.meanviralload.by.gender.loc.age.icar(dall, refit = args$refit)
}

# Estimate suppressed pop
# _______________________

if (args$run.gp.supp.hiv) { # Among HIV positive
    vl.suppofinfected.by.gender.loc.age.gp(dall, refit = args$refit)
}


if (args$run.gp.supp.pop) { # Among Entire population
    vl.suppofpop.by.gender.loc.age.gp(dall, refit = args$refit)
}


if (0) {
    # GET POSTERIORS ON SUPP AMONG POP
    # ________________________________

    .f <- function(file) {
        tmp <- new.env()
        round <- as.integer(gsub("^.*?round([0-9][0-9]).*?$", "\\1", file))
        load(file, envir = tmp)
        ls(tmp)
        tmp <- tmp$nspop.by.age
        tmp[, ROUND := round]
        tmp
    }
    dsupp <- list.files(vl.out.dir, "220729f_suppAmongPop.*rda", full.names = T)
    dsupp <- lapply(dsupp, .f)
    dsupp <- rbindlist(dsupp)



    # CAN WE GET INCIDENCE NOW? (why?)
    # ________________________

    dinc <- file.path(git.repository, "data", "RCCS_1518_incidence.csv")
    dinc <- fread(dinc)

    dinc[, .(INCIDENCE * PY, NEWINF)]

    tmp <- grep("ROUND|COMM|AGEYRS|SEX|INCIDEN|PREVALEN", names(dinc), value = TRUE)
    dinc <- dinc[, ..tmp]

    p <- ggplot(dinc, aes(x = AGEYRS, colour = ROUND)) +
        geom_line(aes(y = INCIDENCE, group = ROUND)) +
        facet_grid(COMM ~ SEX) +
        viridis::scale_color_viridis() +
        theme_bw()


    # Want to
    # 1. get Mean Viral Load by location
    # 2. run an analysis GLM like
    # 3. Using what as a predictor?
}
