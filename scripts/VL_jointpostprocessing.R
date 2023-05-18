################
# DEPENDENCIES #
################

library(data.table)
library(ggplot2)
library(ggtext)
library(ggpubr)
library(knitr)
library(Hmisc)
library(xtable)
library(here)
library(optparse)

################
#    PATHS     #
################

gitdir <- here()
source(file.path(gitdir, "R/paths.R"))

file.exists(
        path.hivstatusvl.r1520,
        path.census.eligible
) |> all() |> stopifnot()

# command line options, stored in args. Then subset
opts_vec <- c('viremic.viral.load', 'detectable.viral.load', 'out.dir.prefix', 'indir','round', 'jobname', 'only.firstparticipants')
args <- args[ names(args) %in% opts_vec]

source( file.path(gitdir.functions,'plotting_main_figures.R') )
source( file.path(gitdir.functions,'postprocessing_helpers.R') )
source( file.path(gitdir.functions,'phsc_vl_helpers.R') )
naturemed_reqs()


# output directories with vl files

make_paper_numbers <- TRUE
if(make_paper_numbers)
    ppr_numbers <- list()

naturemed_reqs()
VL_DETECTABLE = args$vl.detectable
VIREMIC_VIRAL_LOAD = args$viremic.viral.load

############
#   MAIN   #
############

# get model fits for both scenarios: all participants and first participants
args2 <- copy(args)
args2$only.firstparticipants <- ! args$only.firstparticipants
files1 <- list.files.from.output.directory('.rda', args=args, rounds=16:19)
files2 <- list.files.from.output.directory('.rda', args=args2, rounds=16:19)

# get paths of models
dfiles <- data.table(F = c(files1, files2))
dfiles[, `:=` (
    D = dirname(F),
    F = basename(F),
    MODEL = basename(dirname(F)),
    ROUND = gsub('^.*round([0-9]+).*$','\\1',F) |> as.integer(),
    IDX = basename(dirname(dirname(F)))
)]
dfiles[,  c('VL', 'FTP', 'JOB') := fetch.args.from.suffix(.BY), by=IDX ]
dfiles[, IDX := NULL ]
stopifnot(dfiles[, .N, by='F'][, all(N == 2)])

# load proportions of census eligibles
dpartrates <- readRDS(path.participation.rates)

# by round and participant type
tmp <- new.env()
dfiles[, load(file.path(D[1], F[1]), env=tmp)]

# TODO: probably store this somewhere so easy to read? 
# Create a list of environments containing objects required for plots 

env_list <- store.rda.environments.in.list.by.round.ftpstatus(dfiles[ROUND == 16])

rbind.ftpstatus.datatables.by.round <- function(DTname, round, envir_list=env_list){

    if(is.numeric(round)) 
        round <- as.character(round)

    DT.allp <- get(DTname,envir=envir_list[[round]][['allp']])
    DT.ftp <- get(DTname,envir=envir_list[[round]][['ftp']])

    DT.allp[, FTP_LAB := "All participants"]
    DT.ftp[, FTP_LAB := "First-time participants"]
    
    rbind(DT.allp, DT.ftp)
}

# plot HIV prevalence; estimates are similar but ofc more wiggly in ftp
prev.hiv.by.age <- rbind.ftpstatus.datatables.by.round('prev.hiv.by.age', 16, envir_list=env_list)
prev.hiv.by.age |> plot.comparison.ftptype.colftp(ylab="HIV prevalence")
prev.hiv.by.age |> plot.comparison.ftptype.colsex(ylab="HIV prevalence")


# suppression among infected: again similar, no significant differences except for 'olde'
nsinf.by.age <- rbind.ftpstatus.datatables.by.round('nsinf.by.age', 16, envir_list=env_list)
nsinf.by.age |> plot.comparison.ftptype.colftp(ylab="Viral suppression among HIV positives")
nsinf.by.age |> plot.comparison.ftptype.colsex(ylab="Viral suppression among HIV positives")

# viraemia among all
nspop.by.age <- rbind.ftpstatus.datatables.by.round('nspop.by.age', 16, envir_list=env_list)
nspop.by.age |> plot.comparison.ftptype.colftp(ylab="Prevalence of viraemia")
nspop.by.age |> plot.comparison.ftptype.colsex(ylab="Prevalence of viraemia")

