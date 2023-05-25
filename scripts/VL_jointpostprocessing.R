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
library(posterior)

################
#    PATHS     #
################

gitdir <- here()
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
source(file.path(gitdir.functions, "phsc_vl_helpers.R"))
naturemed_reqs()

make_paper_numbers <- TRUE
if (make_paper_numbers) {
    ppr_numbers <- list()
}

VL_DETECTABLE <- args$vl.detectable
VIREMIC_VIRAL_LOAD <- args$viremic.viral.load

debug <- 0

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

# load number of census eligible individuals (.25 too rough)
dpartrates <- readRDS(path.participation.rates) |>
    subset(select = c("ROUND", "FC", "SEX", "AGEYRS", "PARTRATE_SMOOTH.25")) |>
    setnames(c("FC", "PARTRATE_SMOOTH.25"), c("LOC", "PARTRATE"))

# get model fits for both scenarios: all participants and first participants

# let us first work with rda files
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

rbind.ftpstatus.datatables.by.round <- function(DTname, round, envir_list = env_list) {
    if (is.numeric(round)) {
        round <- as.character(round)
    }

    DT.allp <- get(DTname, envir = envir_list[[round]][["allp"]])
    DT.ftp <- get(DTname, envir = envir_list[[round]][["ftp"]])

    DT.allp[, FTP_LAB := "All participants"]
    DT.ftp[, FTP_LAB := "First-time participants"]

    rbind(DT.allp, DT.ftp)
}

round <- 19

# plot HIV prevalence; estimates are similar but ofc more wiggly in ftp
prev.hiv.by.age <- rbind.ftpstatus.datatables.by.round("prev.hiv.by.age", round, envir_list = env_list)
prev.hiv.by.age |> plot.comparison.ftptype.colsex(ylab = "HIV prevalence")
p1 <- prev.hiv.by.age |> plot.comparison.ftptype.colftp(ylab = "HIV prevalence")
filename <- paste0("fit_hivprev_byftpstatus_round", round, ".pdf")
ggsave2(p = p1, file = filename, LALA = out.dir.figures, w = 9, h = 8)


# suppression among infected: again similar, no significant differences except for 'olde'
nsinf.by.age <- rbind.ftpstatus.datatables.by.round("nsinf.by.age", round, envir_list = env_list)
nsinf.by.age |> plot.comparison.ftptype.colsex(ylab = "Viral suppression among HIV positives")
p1 <- nsinf.by.age |> plot.comparison.ftptype.colftp(ylab = "Viral suppression among HIV positives")
filename <- paste0("fit_suppofhiv_byftpstatus_round", round, ".pdf")
ggsave2(p = p1, file = filename, LALA = out.dir.figures, w = 9, h = 8)

# viraemia among all
nspop.by.age <- rbind.ftpstatus.datatables.by.round("nspop.by.age", round, envir_list = env_list)
nspop.by.age |> plot.comparison.ftptype.colsex(ylab = "Prevalence of viraemia")
p1 <- nspop.by.age |> plot.comparison.ftptype.colftp(ylab = "Prevalence of viraemia")
filename <- paste0("fit_suppofpop_byftpstatus_round", round, ".pdf")
ggsave2(p = p1, file = filename, LALA = out.dir.figures, w = 9, h = 8)

rm(round)

# RDS files
# _________

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

# get posteriors for proportions among entire pop, as weighted averages of FTP and non.
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

# test
if(0)
{
    dfiles_rds[ROUND==19 & MODEL == 'run-gp-prevl', {
        stopifnot(.N==2)
        paths <- file.path(D, F)
        idx.all <- which(FTP==FALSE)
        idx.ftp <- which(FTP==TRUE)
        fit1 <<- readRDS(paths[idx.all])
        fit2 <<- readRDS(paths[idx.ftp])
    }]
}


filename <- "fit_hivprev_byroundcommgender.pdf"
p_hiv <- plot.fit.weighted.by.ftpstatus(djoint, "run-gp-prevl")
ggsave2(p = p_hiv, file = filename, LALA = out.dir.figures, w = 10, h = 12)

filename <- "fit_suppofhiv_byroundcommgender.pdf"
p_supp <- plot.fit.weighted.by.ftpstatus(djoint, "run-gp-supp-hiv")
ggsave2(p = p_supp, file = filename, LALA = out.dir.figures, w = 10, h = 12)

filename <- "fit_suppofpop_byroundcommgender.pdf"
p_vir <- plot.fit.weighted.by.ftpstatus(djoint, "run-gp-supp-pop")
ggsave2(p = p_vir, file = filename, LALA = out.dir.figures, w = 10, h = 12)

# Now I want to get the raw numbers: what is the contribution to LALALA

cols <- c("ROUND", "LOC_LABEL", "SEX_LABEL", "AGE_LABEL")
dviraemic <- djoint[MODEL == "run-gp-supp-pop", list(
    ROUND = ROUND,
    SEX_LABEL = stan_dicts$INTtoSEX[as.character(SEX)],
    LOC_LABEL = stan_dicts$INTtoLOC[as.character(LOC)],
    AGE_LABEL = AGE_LABEL,
    CL, IL, M, UL, CU
)] |>
    unique() |>
    merge(dcens, by = cols)

dplot <- dviraemic[,
    lapply(.SD, \(x) x * ELIGIBLE),
    .SDcols = c("CL", "M", "CU"),
    by = cols
]
dplot |> prettify_labels()

filename <- "bar_enum_viraemic_censel.pdf"
p1 <- plot.estimated.number.viraemic.among.census.eligible(dplot)
ggsave2(p = p1, file = filename, LALA = out.dir.figures, w = 10, h = 12)

# NOTE: contribution should be calculated starting from the draws...
# this is slightly annoying but it's doable nonetheless if we really want it
# also, results here are counterintuitive. Am I plotting the right thing???
dplot2 <- dplot[, .(AGE_LABEL, SEX_LAB, P = M / sum(M)), by = c("ROUND_LAB", "LOC_LAB")]
filename <- "line_econtr_viraemic_censel.pdf"
p2 <- plot.estimated.contribution.viraemic.among.census.eligible(dplot2)
ggsave2(p = p2, file = filename, LALA = out.dir.figures, w = 10, h = 12)

dplot2[, sum(P), by = c("ROUND_LAB", "LOC_LAB", "SEX_LAB")] |>
    ggplot(aes(x = ROUND_LAB, y = V1, fill = SEX_LAB)) +
    geom_col() +
    facet_grid(~LOC_LAB) +
    scale_fill_manual(values = palettes$sex) +
    scale_y_percentage +
    theme_default() +
    my_labs(y = "Sex contribution to population viraemia")
