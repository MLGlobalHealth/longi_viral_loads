#!/bin/Rscript

########################
cat("\nStart of: VL_jointpostprocessing.R\n")
########################


################
# DEPENDENCIES #
################
{
    library(data.table)
    library(ggplot2)
    library(ggtext)
    library(ggpubr)
    library(knitr)
    library(Hmisc)
    library(xtable)
    library(here)
    library(optparse)
    library(posterior)
}

################
#    PATHS     #
################

# NOTE: contributions by age/group do not make sense in the "among hiv", if we merge by dcens.
# instead, we should be merging by the N of hiv+

# TODO: 
# - 2. [X] make table with improvements over time in tot # usuppressed by loc and gender
# - 3. [] I am not saving the sex comparison in the  plot.comparison.ftptype.colsex section: do it

gitdir <- here::here()
source(file.path(gitdir, "R/paths.R"))

file.exists(
    path.hivstatusvl.r1520,
    path.census.eligible
) |>
    all() |>
    stopifnot()

# command line options, stored in args. Then subset
opts_vec <- c("viremic.viral.load", "detectable.viral.load", "out.dir.prefix", "indir", "round", "jobname", "only.firstparticipants")
args <- args[names(args) %in% opts_vec]

source(file.path(gitdir.functions, "plotting_main_figures.R"))
source(file.path(gitdir.functions, "postprocessing_helpers.R"))
source(file.path(gitdir.functions, "postprocessing_tables.R"))
source(file.path(gitdir.functions, "phsc_vl_helpers.R"))
naturemed_reqs()

overwrite <- !interactive()
make_plots <- TRUE
make_tables <- TRUE

make_paper_numbers <- TRUE
if (make_paper_numbers) {
    ppr_numbers <- list()
}

VL_DETECTABLE <- args$vl.detectable
VIREMIC_VIRAL_LOAD <- args$viremic.viral.load

# output directories
out.dir <- args$out.dir.prefix
out.dir <- file.path(out.dir, paste0("vl_", VIREMIC_VIRAL_LOAD, "_joint"))
out.dir.figures <- file.path(out.dir, "figures")
out.dir.tables <- file.path(out.dir, "tables")
dir.create(out.dir.tables) |> suppressWarnings()
dir.create(out.dir.figures) |> suppressWarnings()

############
#   MAIN   #
############

# load first participation rates
dfirst_prop <- get.first.participant.rates()

# get census eligible
dcens <- get.census.eligible() |>
    setnames(c("AGE_LABEL", "SEX_LABEL", "LOC_LABEL"), c("AGEYRS", "SEX", "LOC"))
dcens[, AGEGROUP :=  split.agegroup(AGEYRS)]

# load number of census eligible individuals (.50 too rough)
dpartrates <- readRDS(path.participation.rates) |>
    subset(select = c("ROUND", "FC", "SEX", "AGEYRS", "PARTRATE_SMOOTH.25")) |>
    setnames(c("FC", "PARTRATE_SMOOTH.25"), c("LOC", "PARTRATE"))

# get model fits for both scenarios: all participants and first participants

####################################################
catn("=== Compare model fits among FTP and ALL ===")
####################################################

args2 <- copy(args)
args2$only.firstparticipants <- !args$only.firstparticipants
files1 <- list.files.from.output.directory(".rda", args = args, rounds = 16:19)
files2 <- list.files.from.output.directory(".rda", args = args2, rounds = 16:19)

# get paths of models
dfiles_rda <- data.table(F = c(files1, files2))
dfiles_rda[, `:=`(
    D = dirname(F),
    F = basename(F),
    MODEL = basename(dirname(F)),
    ROUND = gsub("^.*round([0-9]+).*$", "\\1", F) |> as.integer(),
    IDX = basename(dirname(dirname(F)))
)]
dfiles_rda[, c("VL", "FTP", "JOB") := fetch.args.from.suffix(.BY), by = IDX]
dfiles_rda[, IDX := NULL]
stopifnot(dfiles_rda[, .N, by = "F"][, all(N == 2)])


# by round and participant type
env_list <- store.rda.environments.in.list.by.round.ftpstatus(dfiles_rda)

if (make_plots) {
    for (round in 16:19) {
        catn("Round:", round)
        # plot HIV prevalence; estimates are similar but ofc more wiggly in ftp
        prev.hiv.by.age <- rbind.ftpstatus.datatables.by.round("prev.hiv.by.age", round, envir_list = env_list)
        prev.hiv.by.age |> plot.comparison.ftptype.colsex(ylab = "HIV prevalence")
        p1 <- prev.hiv.by.age |> plot.comparison.ftptype.colftp(ylab = "HIV prevalence")
        filename <- paste0("fit_hivprev_byftpstatus_round", round, ".pdf")
        ggsave2(p = p1, file = filename, LALA = out.dir.figures, w = 9, h = 8)


        # suppression among infected: again similar, no significant differences except for 'olde'
        nsinf.by.age <- rbind.ftpstatus.datatables.by.round("nsinf.by.age", round, envir_list = env_list)
        nsinf.by.age |> plot.comparison.ftptype.colsex(ylab = "Viral suppression among HIV positives")
        p2 <- nsinf.by.age |> plot.comparison.ftptype.colftp(ylab = "Viral suppression among HIV positives")
        filename <- paste0("fit_suppofhiv_byftpstatus_round", round, ".pdf")
        ggsave2(p = p2, file = filename, LALA = out.dir.figures, w = 9, h = 8)

        # viraemia among all
        nspop.by.age <- rbind.ftpstatus.datatables.by.round("nspop.by.age", round, envir_list = env_list)
        nspop.by.age |> plot.comparison.ftptype.colsex(ylab = "Prevalence of viraemia")
        p3 <- nspop.by.age |> plot.comparison.ftptype.colftp(ylab = "Prevalence of viraemia")
        filename <- paste0("fit_suppofpop_byftpstatus_round", round, ".pdf")
        ggsave2(p = p3, file = filename, LALA = out.dir.figures, w = 9, h = 8)
    }
    rm(round)
}

##################################################
catn("=== Get posterior draws from rds files ===")
##################################################

args2 <- copy(args)
args2$only.firstparticipants <- !args$only.firstparticipants
files1 <- list.files.from.output.directory(".rds", args = args, rounds = 16:19)
files2 <- list.files.from.output.directory(".rds", args = args2, rounds = 16:19)

dfiles_rds <- data.table(F = c(files1, files2))
dfiles_rds[, `:=`(
    D = dirname(F),
    F = basename(F),
    MODEL = basename(dirname(F)),
    ROUND = gsub("^.*round([0-9]+).*$", "\\1", F) |> as.integer(),
    IDX = basename(dirname(dirname(F)))
)]
dfiles_rds[, c("VL", "FTP", "JOB") := fetch.args.from.suffix(.BY), by = IDX]
dfiles_rds[, IDX := NULL]
stopifnot(dfiles_rds[, .N, by = "F"][, all(N == 2)])

catn("Load all previous results and plot") 
#_________________________________________ 

# get posteriors for proportions among entire pop, as weighted averages of FTP and non.
filename_rds <- file.path(out.dir.tables, "fullpop_posteriorquantiles_by_agesexround.rds")

if (file.exists(filename_rds) & !overwrite) {
    djoint <- readRDS(filename_rds)
} else {
    djoint <- dfiles_rds[,
        {
            stopifnot(.N == 2)
            paths <- file.path(D, F)
            idx.all <- which(FTP == FALSE)
            idx.ftp <- which(FTP == TRUE)
            cat(paths[idx.all], "\n")
            get.weighted.average.p_predict(
                readRDS(paths[idx.all]),
                readRDS(paths[idx.ftp]),
                round = unique(ROUND)
            )
        },
        by = c("MODEL", "ROUND")
    ]
    cat("\nSaving ",filename_rds,"\n")
    saveRDS(object = djoint, file = filename_rds)
}

if (make_plots) {
    # set dimensions for all plots below
    .w <- 10; .h <- 12

    p_hiv <- plot.fit.weighted.by.ftpstatus(djoint, "run-gp-prevl")
    p_supp <- plot.fit.weighted.by.ftpstatus(djoint, "run-gp-supp-hiv")
    p_vir <- plot.fit.weighted.by.ftpstatus(djoint, "run-gp-supp-pop")

    .fnm <- function(lab) {
        paste("fit", lab, "byroundcommgender.pdf", sep = "_")
    }

    ggsave2(p = p_hiv, file = .fnm("hivprev"), LALA = out.dir.figures, .w, .h)
    ggsave2(p = p_supp, file = .fnm("suppofhiv"), LALA = out.dir.figures, .w, .h)
    ggsave2(p = p_vir, file = .fnm("suppofpop"), LALA = out.dir.figures, .w, .h)
    
    rm(.w, .h)
}


catn("Get quantiles for population prevalences by agegroup") 
#___________________________________________________________

filename_rds <- file.path(out.dir.figures, "posterior_quantiles_agegroups.rds")

if( file.exists(filename_rds) & !overwrite ){
    djoint_agegroup <- readRDS(filename_rds)
} else {
    djoint_agegroup <- dfiles_rds[ MODEL != 'run-gp-supp-hiv' ,
        {
            stopifnot(.N == 2)
            paths <- file.path(D, F)
            idx.all <- which(FTP == FALSE)
            idx.ftp <- which(FTP == TRUE)
            cat(paths[idx.all], "\n")
            get.weighted.average.p_predict(
                readRDS(paths[idx.all]),
                readRDS(paths[idx.ftp]),
                round = unique(ROUND),
                expression_prereturn = {
                    tmp <- merge(draws_all, dcens, by = c("LOC", "SEX", "AGEYRS", "ROUND"))
                    tmp <- tmp[, 
                            list( joint =  sum(joint * ELIGIBLE_SMOOTH) ), 
                        by = c(dot.cols, "LOC", "SEX", "AGEGROUP")
                    ] 
                    tmp[, quantile2(joint), by = c("SEX", "LOC", "AGEGROUP")]
                }
            )
        },
        by = c("MODEL", "ROUND")
    ]
    saveRDS(object= djoint_agegroup, filename_rds)
}

catn("Get quantiles for population prevalences aggregated over age") 
#___________________________________________________________

filename_rds <- file.path(out.dir.figures, "posterior_quantiles_ageaggregate.rds")

if( file.exists(filename_rds) & !overwrite ){
    joint_ageagrr_list <- readRDS(filename_rds)
} else {
    djoint_ageaggr <- dfiles_rds[ MODEL != 'run-gp-supp-hiv' ,
        {
            stopifnot(.N == 2)
            paths <- file.path(D, F)
            idx.all <- which(FTP == FALSE)
            idx.ftp <- which(FTP == TRUE)
            cat(paths[idx.all], "\n")
            get.weighted.average.p_predict(
                readRDS(paths[idx.all]),
                readRDS(paths[idx.ftp]),
                round = unique(ROUND),
                expression_prereturn = {
                    tmp <- merge(draws_all, dcens, by = c("LOC", "SEX", "AGEYRS", "ROUND"))
                    tmp <- tmp[, 
                            list( joint =  sum(joint * ELIGIBLE_SMOOTH) ), 
                        by = c(dot.cols, "LOC", "SEX")
                    ] 
                }
            )
        },
        by = c("MODEL", "ROUND")
    ]
    .f <- function(old, new)
        (old - new)/old

    tmp <- dcast(djoint_ageaggr, 
        MODEL + .chain + .iteration + .draw + LOC + SEX ~ paste0("R",ROUND), 
        value.var='joint' )[, 
        list(
        MODEL=MODEL, LOC=LOC, SEX=SEX,
        redR17=.f(R16,R17),
        redR18=.f(R16,R18),
        redR19=.f(R16,R19),
        redR1718=.f(R17,R18),
        redR1819=.f(R18,R19)
    )] |> 
        melt( 
            id.vars = c('MODEL','LOC', 'SEX'),
            measure.vars = c("redR17" , "redR18" ,"redR19", "redR1718", "redR1819")
        )

    joint_ageagrr_list <- list(
        percent_reduction = tmp[, quantile2(value), by=c("MODEL","SEX", 'LOC', 'variable')],
        round_totals = djoint_ageaggr[, quantile2(joint), by=c("MODEL","SEX", "LOC", "ROUND") ]
    )
    saveRDS(object= joint_ageagrr_list, filename_rds)
    rm(djoint_ageaggr)
}

if(make_tables){

    tab2 <- tablify.posterior.Nunsuppressed(joint_ageagrr_list)

    filename_tex <- file.path(out.dir.tables, 'table_aggregatedNunsuppressed.tex')
    write.to.tex(tab2, file=filename_tex)
    filename <- 'table_aggregatedNunsuppressed.pdf'
    p <- gridExtra::tableGrob(tab2) |> gridExtra::grid.arrange()
    ggsave2(p=p, file=filename, LALA=out.dir.tables, w=8, h=5.5)
    
}


catn("Get quantiles for contributions by age") 
#_____________________________________________

filename_rds <- file.path(out.dir.tables, "fullpop_allcontributions.rds")

if (file.exists(filename_rds) & !overwrite) {
    dcontrib <- readRDS(filename_rds)
} else {
    dcontrib <- dfiles_rds[,
        {
            stopifnot(.N == 2)
            paths <- file.path(D, F)
            idx.all <- which(FTP == FALSE)
            idx.ftp <- which(FTP == TRUE)

            cat(paths[idx.all], "\n")
            get.weighted.average.p_predict(
                readRDS(paths[idx.all]),
                readRDS(paths[idx.ftp]),
                round = unique(ROUND),
                expression_prereturn = {
                    tmp <- merge(draws_all, dcens, by = c("LOC", "SEX", "AGEYRS", "ROUND"))
                    tmp <- tmp[,
                        {
                            z <- joint * ELIGIBLE_SMOOTH 
                            list(
                                AGEYRS = AGEYRS,
                                SEX = SEX,
                                CONTRIBUTION = z / sum(z)
                            )
                        },
                        by = c(dot.cols, "LOC")
                    ]
                    tmp[, quantile2(CONTRIBUTION), by = c("SEX", "LOC", "AGEYRS")]
                }
            )
        },
        by = c("MODEL", "ROUND")
    ]
    saveRDS(object = dcontrib, file = filename_rds)
}

check_median_contr_approx1(dcontrib)

if (make_plots) {

    .w <- 10; .h <- 12

    p_contrib_prevl <- plot.agesex.contributions.by.roundcomm(dcontrib, label = "run-gp-prevl", include_baseline =TRUE)
    p_contrib_supph <- plot.agesex.contributions.by.roundcomm(dcontrib, label = "run-gp-supp-hiv", include_baseline = TRUE)
    p_contrib_suppp <- plot.agesex.contributions.by.roundcomm(dcontrib, label = "run-gp-supp-pop", include_baseline = TRUE)

    .fnm <- function(lab)
        paste("contrib_agegender", lab, "byroundcomm.pdf", sep = "_")

    ggsave2(p = p_contrib_prevl, file = .fnm("prevl"), LALA = out.dir.figures, .w, .h)
    ggsave2(p = p_contrib_supph, file = .fnm("suppofhiv"), LALA = out.dir.figures, .w, .h)
    ggsave2(p = p_contrib_suppp, file = .fnm("suppofpop"), LALA = out.dir.figures, .w, .h)

}


catn("Get quantiles for contributions by agegroup") 
#__________________________________________________

filename_rds <- file.path(out.dir.tables, "fullpop_allcontributions_byagegroup.rds")
overwrite <- FALSE

if (file.exists(filename_rds) & !overwrite) {
    dcontrib_agegroup <- readRDS(filename_rds)
} else {
    dcontrib_agegroup <- dfiles_rds[,
        {
            stopifnot(.N == 2)
            paths <- file.path(D, F)
            idx.all <- which(FTP == FALSE)
            idx.ftp <- which(FTP == TRUE)

            cat(paths[idx.all], "\n")
            get.weighted.average.p_predict(
                readRDS(paths[idx.all]),
                readRDS(paths[idx.ftp]),
                round = unique(ROUND),
                expression_prereturn = {
                    tmp <- merge(draws_all, dcens, by = c("LOC", "SEX", "AGEYRS", "ROUND"))
                    tmp <- tmp[,
                        {
                            z <- joint * ELIGIBLE_SMOOTH
                            list( z = sum(z))
                        },
                        by = c(dot.cols, "LOC", "SEX", "AGEGROUP")
                    ][, list(
                        AGEGROUP = AGEGROUP,
                        SEX = SEX,
                        CONTRIBUTION = z/sum(z)
                    ),
                    by=c(dot.cols, "LOC")]
                    tmp[, quantile2(CONTRIBUTION), by = c("SEX", "LOC", "AGEGROUP")]
                }
            )
        },
        by = c("MODEL", "ROUND")
    ]
    saveRDS(object = dcontrib_agegroup, file = filename_rds)
}

check_median_contr_approx1(dcontrib_agegroup)

if (make_plots) {

    .w <- 10; .h <- 12

    p_contrib_prevl <- dcontrib_agegroup |> plot.agesex.contributions.by.roundcomm2(label = "run-gp-prevl")
    # p_contrib_supph <- dcontrib_agegroup |> plot.agesex.contributions.by.roundcomm2(label = "run-gp-supp-hiv")
    p_contrib_suppp <- dcontrib_agegroup |> plot.agesex.contributions.by.roundcomm2(label = "run-gp-supp-pop")

    .fnm <- function(lab)
        paste("contrib_agegroupgender", lab, "byroundcomm.pdf", sep = "_")

    ggsave2(p = p_contrib_prevl, file = .fnm("prevl"), LALA = out.dir.figures, .w, .h)
    # ggsave2(p = p_contrib_supph, file = .fnm("suppofhiv"), LALA = out.dir.figures, .w, .h)
    ggsave2(p = p_contrib_suppp, file = .fnm("suppofpop"), LALA = out.dir.figures, .w, .h)
}

if(make_tables){

    t_contrib_supp <- tablify.agecontributions(dcontrib_agegroup, model= 'run-gp-supp-pop')
    t_contrib_supp$ROUND_LAB <- labeller(ROUND_LAB = round_labs)(t_contrib_supp[, .(ROUND_LAB)])

    filename_tex <- file.path(out.dir.tables, 'table_contrib_supppop.tex')
    write.to.tex(t_contrib_supp, file=filename_tex)
    filename <- 'table_contrib_supppop.pdf'
    p <- gridExtra::tableGrob(t_contrib_supp) |> gridExtra::grid.arrange()
    ggsave2(p=p, file=filename, LALA=out.dir.tables, w=22.5, h=5.5)
}


#####################
catn("End of script")
#####################
