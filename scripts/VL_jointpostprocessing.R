#!/bin/Rscript

if (!interactive()) {
    options(error = dump.frames)
}

self_relative_path <- "scripts/VL_jointpostprocessing.R"
########################
cat("\nStart of:", self_relative_path, "\n")
########################


################
# DEPENDENCIES #
################

suppressPackageStartupMessages({
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
})

################
#    PATHS     #
################

# NOTE: contributions by age/group do not make sense in the "among hiv", if we merge by dcens.
# instead, we should be merging by the N of hiv+

if (interactive()) {
    gitdir <- here::here()
} else {
    cmd <- commandArgs()
    cmd <- cmd[cmd %like% "file"]
    gitdir <- gsub(paste0(".*--file=(.*)/", self_relative_path, "$"), "\\1", cmd)
}

{
    source(file.path(gitdir, "R/paths.R"))
    source(file.path(gitdir.functions, "plotting_main_figures.R"))
    source(file.path(gitdir.functions, "postprocessing_helpers.R"))
    source(file.path(gitdir.functions, "postprocessing_tables.R"))
    source(file.path(gitdir.functions, "phsc_vl_helpers.R"))
    source(file.path(gitdir.functions, "paper_statements.R"))
    naturemed_reqs()
}

# command line options, stored in args. Then subset
opts_vec <- c(
    "viremic.viral.load",
    "detectable.viral.load",
    "out.dir.prefix",
    "out.dir.exact",
    "round",
    "jobname",
    "shared.hyper"
)
args <- args[names(args) %in% opts_vec]
# testing
if (interactive() & usr=="andrea" ) {
    # args$jobname <- "vl_1000_firstpart"
    args$jobname <- "cmdstan_alpha100sharedhyper_vl_1000"
    args$out.dir.exact <- "/home/andrea/HPC/ab1820/home/projects/2022/longvl/cmdstan_alpha100sharedhyper_vl_1000"
    args$shared.hyper <- TRUE
}
print(args)

overwrite <- !interactive()
make_plots <- make_tables <- TRUE
# Get out.dir etc...
fetch.postprocessing.settings.from.args(args)

####################
catn("=== MAIN ===")
####################

# dfirst_prop <- get.first.participant.rates()
# plot.first.participant.rates(add_loess=FALSE)
# plot.first.participant.rates(add_loess=TRUE)

# get census-eligible
if (file.exists(path.census.eligible.aggregated)){
    dcens <- fread(path.census.eligible.aggregated)
}else{
    dcens <- get.census.eligible(path=path.census.eligible) |>
        setnames(c("AGE_LABEL", "SEX_LABEL", "LOC_LABEL"), c("AGEYRS", "SEX", "LOC"))
}
dcens[, AGEGROUP := split.agegroup(AGEYRS)]

# load number of census-eligible individuals (.50 too rough)
# not sure if needed ...
dpartrates <- readRDS(path.participation.rates) |>
    subset(select = c("ROUND", "FC", "SEX", "AGEYRS", "PARTRATE_SMOOTH.25")) |>
    setnames(c("FC", "PARTRATE_SMOOTH.25"), c("LOC", "PARTRATE"))

####################################################
catn("=== Compare model fits among FTP and ALL ===")
####################################################

dfiles_rda <- get.output.paths.ftp.and.all(dir.shared = out.dir, regex = ".rda$")

if (!args$shared.hyper) {
    stopifnot(dfiles_rda[, .N, by = "F"][, all(N == 2)])
    env_list <- store.rda.environments.in.list.by.round.ftpstatus(dfiles_rda)
}

if (make_plots & !args$shared.hyper) {
    # for supplementary figure
    require(patchwork)
    p16 <- rbind.ftpstatus.datatables.by.round("nsinf.by.age", 16, envir_list = env_list) |>
        plot.comparison.ftptype.colftp(ylab = "Viral suppression among PLHIV") +
        ggtitle("Round 16")
    p19 <- rbind.ftpstatus.datatables.by.round("nsinf.by.age", 19, envir_list = env_list) |>
        plot.comparison.ftptype.colftp(ylab = "Viral suppression among PLHIV") +
        ggtitle("Round 19")
    p1619 <- p16 + p19
    cmd <- ggsave2(p = p1619, file = "fit_nsinf_byftpstatus_round1619.pdf", LALA = out.dir.figures, w = 14, h = 10)
    # system(zathura2gthumb(cmd))

    for (round in 16:19) {
        catn("Round:", round)
        # plot HIV prevalence; estimates are similar but ofc more wiggly in ftp
        prev.hiv.by.age <- rbind.ftpstatus.datatables.by.round("prev.hiv.by.age", round, envir_list = env_list)
        # g1 <- prev.hiv.by.age |> plot.comparison.ftptype.colsex(ylab = "HIV prevalence")
        p1 <- prev.hiv.by.age |> plot.comparison.ftptype.colftp(ylab = "HIV prevalence")
        filename <- paste0("fit_hivprev_byftpstatus_round", round, ".pdf")
        ggsave2(p = p1, file = filename, LALA = out.dir.figures, w = 9, h = 8)
        # suppression among infected: again similar, no significant differences except for 'olde'
        nsinf.by.age <- rbind.ftpstatus.datatables.by.round("nsinf.by.age", round, envir_list = env_list)
        # g2 <- nsinf.by.age |> plot.comparison.ftptype.colsex(ylab = "Viral suppression among PLHIV")
        p2 <- nsinf.by.age |> plot.comparison.ftptype.colftp(ylab = "Viral suppression among PLHIV")
        filename <- paste0("fit_suppofhiv_byftpstatus_round", round, ".pdf")
        ggsave2(p = p2, file = filename, LALA = out.dir.figures, w = 9, h = 8)
        # viraemia among all
        nspop.by.age <- rbind.ftpstatus.datatables.by.round("nspop.by.age", round, envir_list = env_list)
        # g3 <- nspop.by.age |> plot.comparison.ftptype.colsex(ylab = "Prevalence of viraemia")
        p3 <- nspop.by.age |> plot.comparison.ftptype.colftp(ylab = "Prevalence of viraemia")
        filename <- paste0("fit_suppofpop_byftpstatus_round", round, ".pdf")
        ggsave2(p = p3, file = filename, LALA = out.dir.figures, w = 9, h = 8)
    }
    rm(round)
}

# particular focus on suppression among HIV positivss
if (make_plots & !args$shared.hyper) {
    x <- as.character(c(16, 19))
    tmp <- lapply(x, rbind.ftpstatus.datatables.by.round, DTname = "nsinf.by.age")
    names(tmp) <- x
    dplot <- rbindlist(tmp, idcol = "ROUND", use.names = TRUE) |>
        set(j = c("SEX", "LOC", "Var2", "IL", "IU", "STAT"), value = NULL)
    prettify_labels(dplot)

    raw <- lapply(x, rbind.ftpstatus.datatables.by.round, DTname = "DT") |>
        rbindlist() |>
        prettify_labels()
    p1 <- ggplot(dplot, aes(x = AGE_LABEL, y = M, ymin = CL, ymax = CU, fill = FTP_LAB)) +
        geom_line() +
        geom_ribbon(alpha = .3) +
        theme_default() +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_percentage +
        scale_color_manual(values = palettes$ftp) +
        scale_fill_manual(values = palettes$ftp) +
        facet_grid(LOC_LAB + SEX_LAB ~ ROUND_LAB) +
        my_labs(y = "Prevalence of suppression among HIV positives") +
        NULL
    p2 <- raw |>
        subset(ROUND == 19 & LOC_LABEL == "fishing" & SEX_LABEL == "M") |>
        ggplot(aes(x = AGE_LABEL, y = 1 - VLNS_N / HIV_N, color = FTP_LAB, size = HIV_N)) +
        # geom_point() +
        geom_label(aes(label = HIV_N)) +
        scale_color_manual(values = palettes$ftp) +
        my_labs(size = "Denominator") +
        theme_default() +
        NULL
    require(patchwork)
    p1 + p2

    # if(MODEL == 'run-gp-supp-hiv'){
    #     geom_texthline(
    #         yintercept=.95^3, color='red', linetype='dashed',
    #         label="UNAIDS 95-95-95", vjust=0.5, hjust=1)
    # }

    p <- ggplot(dplot, aes(x = AGE_LABEL, y = M, linetype = FTP_LAB, color = SEX_LAB)) +
        geom_line() +
        theme_default() +
        facet_grid(ROUND_LAB ~ LOC_LAB,
            labeller = labeller(
                LOC_LAB = community_dictionary$longest2,
                ROUND_LAB = round_labs
            )
        ) +
        scale_x_continuous(expand = c(0, 0), limits = c(15, 49), breaks = seq(15, 50, by = 5)) +
        scale_y_percentage +
        scale_color_manual(values = palettes$sex, labels = sex_dictionary2) +
        scale_fill_manual(values = palettes$sex, labels = sex_dictionary2) +
        my_labs(y = "Prevalence of suppression among HIV positives", fill = "Gender", linetype = "") +
        NULL
    p2 <- p + geom_ribbon(alpha = .1, aes(ymin = CL, ymax = CU, fill = SEX_LAB), color = NA)
    filenames <- paste0("fit_suppofhiv_compare_ftpvsall", c("", "CrIs"), "_r1619.pdf")
    cmd <- ggsave2(p = p, file = filenames[1], LALA = out.dir.figures, w = 16, h = 14)
    cmd <- ggsave2(p = p2, file = filenames[2], LALA = out.dir.figures, w = 16, h = 14)
}


##################################################
catn("=== Get posterior draws from rds files ===")
##################################################

dfiles_rds <- get.output.paths.ftp.and.all("round1[0-9].rds$|220729.rds$") |> 
    unique() |>
    subset(MODEL != "tables")
if (!args$shared.hyper) {
    stopifnot(dfiles_rds[, .N, by = "F"][, all(N == 2)])
    stopifnot(dfiles_rds[, .N, by = "MODEL"][, all(N == 8)])
} else {
    stopifnot(dfiles_rds[, .N, by = "F"][, all(N == 1)])
    stopifnot(dfiles_rds[, .N, by = "MODEL"][, all(N == 4)])
}

catn("Load all previous results and plot")
# _________________________________________

# get posteriors for proportions among entire pop, as weighted averages of FTP and non.
filename_rds  <- file.path(out.dir.tables, "fullpop_posteriorquantiles_by_agesexround.rds")

if (file.exists(filename_rds) & !overwrite) {
    djoint <- readRDS(filename_rds)
} else {
    djoint <- dfiles_rds[,
        {
            eval(expr_setup_ftp_all)
            get.weighted.average.p_predict(
                fit1 = readRDS(paths[idx.all]),
                fit2 = readRDS(paths[idx.ftp]),
                fit_all = readRDS(paths),
                round = unique(ROUND)
            )
        },
        by = c("MODEL", "ROUND")
    ]
    cat("\nSaving ", filename_rds, "\n")
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

    # need to compare suppression in Round 19 in Inland and fishing communities...
    .w <- 10; .h <- 7
    p1_vir <- plot.comparison.prevalence.fishinginland.oneround(DT = djoint, model = "run-gp-supp-pop", round = 19)
    p1_supp <- plot.comparison.prevalence.fishinginland.oneround(DT = djoint, model = "run-gp-supp-hiv", round = 19, ylim = 1)
    p1_hiv <- plot.comparison.prevalence.fishinginland.oneround(DT = djoint, model = "run-gp-prevl", round = 19)
    .fnm <- function(lab) {
        paste("fit", lab, "comparefishinland_bygender_round19.pdf", sep = "_")
    }
    ggsave2(p = p1_hiv, file = .fnm("hivprev"), LALA = out.dir.figures, .w, .h)
    ggsave2(p = p1_supp, file = .fnm("suppofhiv"), LALA = out.dir.figures, .w, .h)
    ggsave2(p = p1_vir, file = .fnm("suppofpop"), LALA = out.dir.figures, .w, .h)

    if (exists("env_list")) {
        p <- plot.comparison.ftp.nonftp.and.all(env = env_list, DT = djoint, model = "run-gp-supp-hiv")
        filename <- "fit_suppofhiv_compare_ftpallandjoint.pdf"
        ggsave2(p = p, file = filename, LALA = out.dir.figures, w = 9, h = 8)
    }

    .w <- 18; .h <- 16
    p1 <- plot_propofpop_of_viraemic_byagesex_stratbycommround(djoint, colorby = "ROUND", cri = TRUE)
    p2 <- plot_propofpop_of_viraemic_byagesex_stratbycommround(djoint, colorby = "ROUND", cri = FALSE)
    filename <- paste0("propofpop_of_viraemic_byagesex_stratbycommround", c("_cri", ""), ".pdf")
    cmd <- ggsave2(p = p1, file = filename[1], LALA = out.dir.figures, w = .w, h = .w, u = "cm")
    cmd <- ggsave2(p = p2, file = filename[2], LALA = out.dir.figures, w = .w, h = .w, u = "cm")

    rm(.w, .h)
}

if (make_tables) {
    # .null <- paper_statements_viraemic_among_hiv(negate=FALSE)
    .null <- paper_statements_viraemic_among_hiv(negate = TRUE, age = 30, rounds = c(16, 19))
    .null <- paper_statements_viraemic_among_hiv(negate = TRUE, age = 25, rounds = 19)
}

catn("Get quantiles for population prevalences by agegroup")
# ___________________________________________________________

filename_rds <- file.path(out.dir.tables, "posterior_quantiles_agegroups.rds")

if (file.exists(filename_rds) & !overwrite) {
    djoint_agegroup <- readRDS(filename_rds)
} else {
    djoint_agegroup <- dfiles_rds[MODEL != "run-gp-supp-hiv",
        {
            eval(expr_setup_ftp_all)
            get.weighted.average.p_predict(
                fit1 = readRDS(paths[idx.all]),
                fit2 = readRDS(paths[idx.ftp]),
                fit_all = readRDS(paths),
                round = unique(ROUND),
                expression_prereturn = {
                    # use grouping sets instead of standard 'by' to allow for sex-totals
                    by_cols <- c(dot.cols, "LOC", "SEX", "AGEGROUP")
                    tmp <- merge(draws_all, dcens, by = c("LOC", "SEX", "AGEYRS", "ROUND"))
                    tmp <- groupingsets(tmp,
                        by = by_cols,
                        j = list(N = sum(joint * ELIGIBLE_SMOOTH), NE = sum(ELIGIBLE_SMOOTH)),
                        sets = list(by_cols, setdiff(by_cols, "AGEGROUP"), setdiff(by_cols, c("AGEGROUP", "SEX")))
                    )
                    tmp[, quantile2(N / NE), by = c("LOC", "SEX", "AGEGROUP")]
                }
            )
        }, by = c("MODEL", "ROUND")
    ]
    djoint_agegroup[is.na(AGEGROUP), AGEGROUP := "Total"]
    djoint_agegroup[is.na(SEX), SEX := "Total"]
    saveRDS(object = djoint_agegroup, filename_rds)
}

if (make_tables) {
    .null <- paper_statements_overall_prevalence(round = 19)
    .null <- paper_statements_female_prevalence(djoint_agegroup)
    .null <- paper_statements_prevalence_viraemia(djoint_agegroup, model="run-gp-supp-pop")
    # .null <- paper_statements_prevalence_viraemia(djoint_agegroup, model="run-gp-supp-hiv")
    .null <- paper_statements_prevalence_viraemia2(model = "run-gp-supp-pop")
    # NOT WORK .null <- paper_statements_prevalence_viraemia2(model = 'run-gp-supp-hiv')
    rm(.null)

    # For CROI:
    tmp <- djoint_agegroup[AGEGROUP == "Total" & SEX != "Total" & ROUND %in% c(16, 19)]
    tmp[,  CELL := prettify_cell(M*100, CL*100, CU*100, percent=TRUE)]
    tmp[, dcast(.SD, ROUND ~ SEX, value.var="CELL") , by=c( "MODEL", "LOC")]
}

catn("Get quantiles for population prevalences aggregated over age")
# ___________________________________________________________

filename_rds  <- file.path(out.dir.tables, "posterior_quantiles_ageaggregate.rds")
filename_rds2 <- file.path(out.dir.tables, "posterior_pvalues_virred_comparison.rds")

if (file.exists(filename_rds) & !overwrite) {
    joint_ageagrr_list <- readRDS(filename_rds)
    pval_red <- readRDS(filename_rds2)
} else {
    djoint_ageaggr <- dfiles_rds[MODEL != "run-gp-supp-hiv",
        {
            eval(expr_setup_ftp_all)
            get.weighted.average.p_predict(
                fit1 = readRDS(paths[idx.all]),
                fit2 = readRDS(paths[idx.ftp]),
                fit_all = readRDS(paths),
                round = unique(ROUND),
                expression_prereturn = {
                    tmp <- merge(draws_all, dcens, by = c("LOC", "SEX", "AGEYRS", "ROUND"))
                    tmp <- tmp[,
                        list(joint = sum(joint * ELIGIBLE_SMOOTH)),
                        by = c(dot.cols, "LOC", "SEX")
                    ]
                }
            )
        },
        by = c("MODEL", "ROUND")
    ]
    .percent.reduction <- function(old, new) {
        (old - new) / old
    }
    red <- dcast(djoint_ageaggr,
        MODEL + .chain + .iteration + .draw + LOC + SEX ~ paste0("R", ROUND),
        value.var = "joint"
    )[j=list(
            MODEL = MODEL, LOC = LOC, SEX = SEX,
            redR17 = .percent.reduction(R16, R17),
            redR18 = .percent.reduction(R16, R18),
            redR19 = .percent.reduction(R16, R19),
            redR1718 = .percent.reduction(R17, R18),
            redR1819 = .percent.reduction(R18, R19)
        )
    ] |> melt(
        id.vars = c("MODEL", "LOC", "SEX"),
        measure.vars = c("redR17", "redR18", "redR19", "redR1718", "redR1819")
    )
    pval_red <- red[MODEL %like% "run-gp-supp-pop" & variable == "redR19", {
        iter_F <- as.logical(seq_len(.N) %% 2)
        list(
            type = "Viraemia reduction larger in women rather than men",
            p_val = mean(value[iter_F] > value[!iter_F])
        ) }, by="LOC"]
    joint_ageagrr_list <- list(
        percent_reduction = red[, quantile2(value), by = c("MODEL", "SEX", "LOC", "variable")],
        round_totals = djoint_ageaggr[, quantile2(joint), by = c("MODEL", "SEX", "LOC", "ROUND")]
    )
    saveRDS(object = joint_ageagrr_list, filename_rds)
    saveRDS(object = pval_red, filename_rds2)
    rm(djoint_ageaggr)
}

pval_red

catn("Get quantiles for suppression levels in agegroups")
# _______________________________________________________

filename_rds <- file.path(out.dir.tables, "posterior_quantiles_suppression_agegroup.rds")

if (file.exists(filename_rds) & !overwrite) {
    dsupp_agegroup <- readRDS(filename_rds)
} else {
    dsupp_agegroup <- dfiles_rds[,
        {
            paths <- file.path(D, F)
            .check <- function(x) {
                stopifnot(length(x) == 1 | args$shared.hyper)
                return(x)
            }
            eval(expr_setup_hivprev_supphiv_ftpall)
            cat("Round ", unique(ROUND), "\n")
            draws_prev <- get.weighted.average.p_predict(
                fit1 = readRDS(paths[idx.hivprev.all]),
                fit2 = readRDS(paths[idx.hivprev.ftp]),
                fit_all = readRDS(paths[idx.hivprev]),
                round = unique(ROUND),
                expression_prereturn = {
                    # find composition of prevalnce by age group
                    by_cols <- c(dot.cols, "LOC", "SEX", "AGEGROUP")
                    tmp <- merge(draws_all, dcens, by = c("LOC", "SEX", "AGEYRS", "ROUND"))
                    tmp <- tmp[, N_HIV := joint * ELIGIBLE_SMOOTH] |>
                        subset(select = c(dot.cols, "LOC", "SEX", "AGEGROUP", "AGEYRS", "N_HIV"))
                    return(tmp)
                }
            )
            cat("prevalence done\n")
            draws_supp <- get.weighted.average.p_predict(
                fit1 = readRDS(paths[idx.supphiv.all]),
                fit2 = readRDS(paths[idx.supphiv.ftp]),
                fit_all = readRDS(paths[idx.supphiv]),
                round = unique(ROUND),
                expression_prereturn = {
                    draws_all
                }
            )
            cat("suppression done\n")
            dot.cols <- c(".chain", ".iteration", ".draw")
            tmp <- merge(draws_supp, draws_prev, by = c(dot.cols, "LOC", "SEX", "AGEYRS"))
            .aggr.hiv.prev <- function(DT, by_cols) {
                DT[, .(S = sum(joint * proportions(N_HIV))), by = c(dot.cols, by_cols)][, quantile2(S), by = by_cols]
            }
            list(
                .aggr.hiv.prev(tmp, by_cols = c("LOC", "SEX", "AGEGROUP")),
                .aggr.hiv.prev(tmp, by_cols = c("LOC", "SEX"))[, AGEGROUP := "Total"],
                .aggr.hiv.prev(tmp, by_cols = c("LOC"))[, `:=`(SEX = "Total", AGEGROUP = "Total")]
            ) |> rbindlist(use.names = TRUE)
        },
        by = c("ROUND")
    ]
    saveRDS(object = dsupp_agegroup, filename_rds)
}

filename_rds2 <- file.path(out.dir.tables, "posterior_pvalues_virofplhiv_comparison.rds")
if (file.exists(filename_rds2) & !overwrite) {
    NULL
} else {
    pvalues_draws <- dfiles_rds[ ROUND %in% c(19),
        {
            paths <- file.path(D, F)
            .check <- function(x) {
                stopifnot(length(x) == 1 | args$shared.hyper)
                return(x)
            }
            eval(expr_setup_hivprev_supphiv_ftpall)
            cat("Round ", unique(ROUND), "\n")
            draws_prev <- get.weighted.average.p_predict(
                fit1 = readRDS(paths[idx.hivprev.all]),
                fit2 = readRDS(paths[idx.hivprev.ftp]),
                fit_all = readRDS(paths[idx.hivprev]),
                round = unique(ROUND),
                expression_prereturn = {
                    # find composition of prevalnce by age group
                    by_cols <- c(dot.cols, "LOC", "SEX", "AGEGROUP")
                    tmp <- merge(draws_all, dcens, by = c("LOC", "SEX", "AGEYRS", "ROUND"))
                    tmp <- tmp[, N_HIV := joint * ELIGIBLE_SMOOTH] |>
                        subset(select = c(dot.cols, "LOC", "SEX", "AGEGROUP", "AGEYRS", "N_HIV"))
                    return(tmp)
                }
            )
            cat("prevalence done\n")
            draws_supp <- get.weighted.average.p_predict(
                fit1 = readRDS(paths[idx.supphiv.all]),
                fit2 = readRDS(paths[idx.supphiv.ftp]),
                fit_all = readRDS(paths[idx.supphiv]),
                round = unique(ROUND),
                expression_prereturn = {
                    draws_all
                }
            )
            cat("suppression done\n")
            dot.cols <- c(".chain", ".iteration", ".draw")
            tmp <- merge(draws_supp, draws_prev, by = c(dot.cols, "LOC", "SEX", "AGEYRS"))
            tmp[, .(
                S = sum(joint * proportions(N_HIV))
            ), by = c(dot.cols, c("LOC", "SEX"))]
        },
        by = c("ROUND")
    ]

    pvalues <- pvalues_draws[ROUND == 19, {
        # stopifnot(.N == 2)
        length(S[SEX == "F"])
        z <- mean(S[SEX == "F"] > S[SEX == "M"])
        if (z == 1){
            sprintf(">= 1 - 1/16,000")
        }else{ as.character(z) }
    }, by = c("LOC", "ROUND")]
    pvalues[, TYPE := "Difference in suppression among PLHIV (pd women - men)"]

    saveRDS(object = pvalues, filename_rds2)
    rm(pvalues_draws)
}


if (make_plots & FALSE) {
    .w <- 14; .h <- 14
    p_list <- lapply(c(16:19), plot.prevalence.by.age.group, DT = dsupp_agegroup)
    filenames <- paste0("plot_suppofhiv_by_agegroup_round", c(16:19), ".pdf")
    # for (i in 1:length(p_list)) {
    #     # ggsave2(p=p_list[[i]], file=filenames[i], LALA=out.dir.figures, .w, .h)
    # }
}

if (make_tables) {
    # paper_statements_suppression_PLHIV_aggregated()
    .null <- paper_statements_suppression_PLHIV_aggregated(reverse = TRUE, round = 19)
    .null <- paper_statements_suppression_PLHIV_aggregated(reverse = TRUE, round = 16)
}

catn("Get quantiles for population prevalences by custom age group")
# __________________________________________________________________

dcens_custom <- copy(dcens)
dcens_custom[, AGEGROUP := split.agegroup(AGEYRS,
    breaks = c(15, 18, 23, 28, 33, 38, 43, 48, 50))]

filename_rds <- file.path(out.dir.tables, "posterior_quantiles_suppression_agegroup_custom.rds")

if (file.exists(filename_rds) & !overwrite) {
    dsupp_agegroup_custom <- readRDS(filename_rds)
} else {
    dsupp_agegroup_custom <- dfiles_rds[,
        {
            paths <- file.path(D, F)
            .check <- function(x) {
                stopifnot(length(x) == 1 | args$shared.hyper)
                return(x)
            }
            eval(expr_setup_hivprev_supphiv_ftpall)
            cat("Round ", unique(ROUND), "\n")
            draws_prev <- get.weighted.average.p_predict(
                fit1 = readRDS(paths[idx.hivprev.all]),
                fit2 = readRDS(paths[idx.hivprev.ftp]),
                fit_all = readRDS(paths[idx.hivprev]),
                round = unique(ROUND),
                expression_prereturn = {
                    # find composition of prevalnce by age group
                    by_cols <- c(dot.cols, "LOC", "SEX", "AGEGROUP")
                    tmp <- merge(draws_all, dcens_custom, by = c("LOC", "SEX", "AGEYRS", "ROUND"))
                    tmp <- tmp[, N_HIV := joint * ELIGIBLE_SMOOTH] |>
                        subset(select = c(dot.cols, "LOC", "SEX", "AGEGROUP", "AGEYRS", "N_HIV"))
                    return(tmp)
                }
            )
            cat("prevalence done\n")
            draws_supp <- get.weighted.average.p_predict(
                fit1 = readRDS(paths[idx.supphiv.all]),
                fit2 = readRDS(paths[idx.supphiv.ftp]),
                fit_all = readRDS(paths[idx.supphiv]),
                round = unique(ROUND),
                expression_prereturn = {
                    draws_all
                }
            )
            cat("suppression done\n")
            dot.cols <- c(".chain", ".iteration", ".draw")
            tmp <- merge(draws_supp, draws_prev, by = c(dot.cols, "LOC", "SEX", "AGEYRS"))
            .aggr.hiv.prev <- function(DT, by_cols) {
                d <- DT[, .(S = sum(joint * proportions(N_HIV))), by = c(dot.cols, by_cols)]
                d[, quantile2(S), by = by_cols]
            }
            list(
                .aggr.hiv.prev(tmp, by_cols = c("LOC", "SEX", "AGEGROUP"))
            ) |> rbindlist(use.names = TRUE)
        },
        by = c("ROUND")]
    saveRDS(object = dsupp_agegroup_custom, filename_rds)
}

# I would then need to get a round 19 histogram...
if (make_plots) {
    .w <- 12; .h <- 18
    # p_list <- lapply(c(16:19), plot.prevalence.by.age.group, DT = dsupp_agegroup)
    p <- hist_prevalence_by_age_group_custom(dsupp_agegroup_custom)
    filename <- "hist_suppofhiv_by_agegroup_custom_round19.pdf"
    cmd <- ggsave2(p=p, file=filename, LALA=out.dir.figures, .w, .h, u="cm")
    system(cmd)
}

if ( make_tables ){
    .null <- paper_statements_prevalence_viraemia_maximum(dsupp_agegroup_custom[SEX%like% "F"])
}


catn("Get quantiles for contributions by age")
    
# _____________________________________________

filename_rds <- file.path(out.dir.tables, "fullpop_allcontributions.rds")
if (file.exists(filename_rds) & !overwrite) {
    dcontrib <- readRDS(filename_rds)
} else {
    dcontrib <- dfiles_rds[,
        {
            eval(expr_setup_ftp_all)
            get.weighted.average.p_predict(
                fit1 = readRDS(paths[idx.all]),
                fit2 = readRDS(paths[idx.ftp]),
                fit_all = readRDS(paths),
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

# check_median_contr_approx1(dcontrib)

if (make_plots) {
    .w <- 10
    .h <- 12
    p_contrib_prevl <- plot.agesex.contributions.by.roundcomm(dcontrib, label = "run-gp-prevl", include_baseline = TRUE)
    p_contrib_supph <- plot.agesex.contributions.by.roundcomm(dcontrib, label = "run-gp-supp-hiv", include_baseline = TRUE)
    p_contrib_suppp <- plot.agesex.contributions.by.roundcomm(dcontrib, label = "run-gp-supp-pop", include_baseline = TRUE)
    .fnm <- function(lab) {
        paste("contrib_agegender", lab, "byroundcomm.pdf", sep = "_")
    }
    ggsave2(p = p_contrib_prevl, file = .fnm("prevl"), LALA = out.dir.figures, .w, .h)
    ggsave2(p = p_contrib_supph, file = .fnm("suppofhiv"), LALA = out.dir.figures, .w, .h)
    ggsave2(p = p_contrib_suppp, file = .fnm("suppofpop"), LALA = out.dir.figures, .w, .h)
}


catn("Get quantiles for sex-contributions by age")
# _____________________________________________

# This first paragraph is not strictly necessary, but it was used to determine the range of ages 
# which contribute to ~ 50 % of the PLHIV population
filename_rds <- file.path(out.dir.tables, "fullpop_allcontributions_sex.rds")
if (file.exists(filename_rds) & !overwrite) {
    dcontrib_sex <- readRDS(filename_rds)
} else {
    dcontrib_sex <- dfiles_rds[ MODEL != "run-gp-supp-hiv",
        {
            eval(expr_setup_ftp_all)
            get.weighted.average.p_predict(
                fit1 = readRDS(paths[idx.all]),
                fit2 = readRDS(paths[idx.ftp]),
                fit_all = readRDS(paths),
                round = unique(ROUND),
                expression_prereturn = {
                    tmp <- merge(draws_all, dcens, by = c("LOC", "SEX", "AGEYRS", "ROUND"))
                    tmp <- tmp[,
                        {
                            z <- joint * ELIGIBLE_SMOOTH
                            list(
                                AGEYRS = AGEYRS,
                                CONTRIBUTION = z / sum(z)
                            )
                        },
                        by = c(dot.cols, "LOC", "SEX")
                    ]
                    tmp[, quantile2(CONTRIBUTION), by = c("SEX", "LOC", "AGEYRS")]
                }
            )
        },
        by = c("MODEL", "ROUND")
    ]
    saveRDS(object = dcontrib_sex, file = filename_rds)
}

if(make_tables){
    # these were the results 
    plhiv_contributors <- data.table(
        SEX = c("F", "F", "M", "M"),
        LOC = c("fishing", "inland", "fishing", "inland"),
        MIN = c(27, 27, 30, 34),
        MAX = c(37, 39, 39, 45)
    )

    dcontrib_50p_PLHIV <- dfiles_rds[ ROUND == 19 & MODEL == "run-gp-prevl",
        {
            eval(expr_setup_ftp_all)
            get.weighted.average.p_predict(
                fit1 = readRDS(paths[idx.all]),
                fit2 = readRDS(paths[idx.ftp]),
                fit_all = readRDS(paths),
                round = unique(ROUND),
                expression_prereturn = {
                    tmp <- merge(draws_all, dcens, by = c("LOC", "SEX", "AGEYRS", "ROUND"))
                    tmp <- merge(tmp, plhiv_contributors, by=c("SEX", "LOC"))
                    tmp[, INAGEGROUP := (AGEYRS >= MIN & AGEYRS <= MAX)]
                    tmp <- tmp[, list(z = sum(joint * ELIGIBLE_SMOOTH)), by = c(dot.cols, "LOC", "SEX", "INAGEGROUP")]
                    tmp <- tmp[,
                        {
                            list(
                                INAGEGROUP = INAGEGROUP,
                                CONTRIBUTION = z / sum(z)
                            )
                        },
                        by = c(dot.cols, "LOC", "SEX")
                    ]
                    tmp[, quantile2(CONTRIBUTION), by = c("SEX", "LOC", "INAGEGROUP")] |>
                        subset(INAGEGROUP == TRUE)
                }
            )
        },
        by = c("MODEL", "ROUND")
    ]
    print_statements_half_plhiv(dcontrib_50p_PLHIV, plhiv_contributors)
}

catn("Get average age for PLHIV & viraemic population")
# _____________________________________________________

filename_rds <- file.path(out.dir.tables, "mean_ages_plhiv_viraemic.rds")

if (file.exists(filename_rds) & !overwrite) {
    dmeanage <- readRDS(filename_rds)
} else {

    .linearly_interpolate <- function(value, y, x){
        idx <- which(y >= value)[1]
        gradient = y[idx] - y[idx - 1]
        b = ( y[idx] - value ) / gradient
        return(x[idx] - b )
    }

    dmeanage <- dfiles_rds[MODEL %in% c("run-gp-prevl", "run-gp-supp-pop"),
        {
            eval(expr_setup_ftp_all)
            get.weighted.average.p_predict(
                fit1 = readRDS(paths[idx.all]),
                fit2 = readRDS(paths[idx.ftp]),
                fit_all = readRDS(paths),
                round = unique(ROUND),
                expression_prereturn = {
                    tmp <- merge(draws_all, dcens, by = c("LOC", "SEX", "AGEYRS", "ROUND"))
                    # Not sure if valid way to compute the posterior "mode"...
                    tmp[, cum_cdf := cumsum(joint * ELIGIBLE_SMOOTH) / sum(joint * ELIGIBLE_SMOOTH),
                        by = c(dot.cols, "LOC", "ROUND", "SEX")]
                    tmp1 <- tmp[
                        j = .(
                            AGEMODE = AGEYRS[which.max(joint * ELIGIBLE_SMOOTH )],
                            AGEMEAN = Hmisc::wtd.mean(AGEYRS, joint * ELIGIBLE_SMOOTH),
                            AGESTD = Hmisc::wtd.var(AGEYRS, joint * ELIGIBLE_SMOOTH) |> sqrt(),
                            AGE25 =  .linearly_interpolate(.25, y=cum_cdf, x=AGEYRS),
                            AGE50 =  .linearly_interpolate(.50, y=cum_cdf, x=AGEYRS),
                            AGE75 =  .linearly_interpolate(.75, y=cum_cdf, x=AGEYRS)
                        ),
                        by = c(dot.cols, "LOC", "ROUND", "SEX")
                    ] |> melt.data.table(
                            measure.vars = c("AGEMEAN", "AGESTD", "AGEMODE", "AGE25", "AGE50", "AGE75"),
                            variable.name = "TYPE",
                            value.name = "value"
                        )
                    tmp1[,j = quantile2(value), by = c("SEX", "LOC", "TYPE")]
                }
            )
        },
        by = c("MODEL", "ROUND")
    ]
    saveRDS(object = dmeanage, file = filename_rds)

    rm(.linearly_interpolate)
}

# TODO! (check that this worked.)
if (make_tables) {
    # dt <- subset(dmeanage, select=-AGEMODE)
    tmp1 <- rbind(
        paper_statements_meanage_population(DT = dmeanage, label = "run-gp-supp-pop", type = "AGEMEAN"),
        paper_statements_meanage_population(DT = dmeanage, label = "run-gp-supp-pop", type = "AGESTD")
    )
    rm(tmp1)

    tab <- make.supp.table.meanage(newline=TRUE, include_date=TRUE)
}

catn("Other contribution to viraemia quantiles for text")
# ________________________________________________________

filename_overleaf <- file.path(out.dir.tables, "overleaf_viraemiacontribution_custom.rds")

dcens_custom <- copy(dcens)
dcens_custom[, AGEGROUP := split.agegroup(AGEYRS, breaks = c(15, 25, 35 ,40, 50))]

if (file.exists(filename_overleaf) & !overwrite) {
    contrib_viraemia_custom <- readRDS(filename_overleaf)
} else {
    contrib_viraemia_custom <- dfiles_rds[MODEL == "run-gp-supp-pop",
        {
            eval(expr_setup_ftp_all)
            get.weighted.average.p_predict(
                fit1 = readRDS(paths[idx.all]),
                fit2 = readRDS(paths[idx.ftp]),
                fit_all = readRDS(paths),
                round = unique(ROUND),
                expression_prereturn = {
                    tmp <- merge(draws_all, dcens_custom, by = c("LOC", "SEX", "AGEYRS", "ROUND"))
                    tmp <- tmp[,
                        {
                            z <- joint * ELIGIBLE_SMOOTH
                            list(z = sum(z))
                        },
                        by = c(dot.cols, "LOC", "SEX", "AGEGROUP")
                    ][, list(
                        AGEGROUP = AGEGROUP,
                        SEX = SEX,
                        CONTRIBUTION = z / sum(z)
                    ),
                    by = c(dot.cols, "LOC")
                    ]
                    tmp[, quantile2(CONTRIBUTION), by = c("SEX", "LOC", "AGEGROUP")]
                }
            )
        },
        by = c("MODEL", "ROUND")
    ]


    # save
    saveRDS(object = contrib_viraemia_custom, file = filename_overleaf)

    # paper_statements_contributions_viraemia_round(round=16, agegroup="15-24")
    # paper_statements_contributions_viraemia_round(round=19, agegroup="15-24")
    paper_statements_contributions_viraemia_round(round=19, agegroup="25-34")
}

if (make_tables) {

    paper_statements_contributions_census_eligible(agegroup = '15-24', comm="inland", sex='M', smooth=FALSE)
    paper_statements_contributions_census_eligible(agegroup = '25-39', comm="inland", sex='M', smooth=FALSE)
    # paper_statements_contributions_census_eligible(agegroup = '25-34', comm="inland", sex='M', smooth=FALSE)

    # paper_statements_contributions_viraemia_round(round = 16)
    paper_statements_contributions_viraemia_round(round = 19, agegroup="25-39")
    paper_statements_contributions_viraemia_round(round = 19, agegroup="15-24")
    contrib_viraemia_custom |> plot_quantiles(x = AGEGROUP, color = SEX_LAB, facet = LOC_LAB ~ ROUND_LAB)
}

catn("Other contribution to PLHIV quantiles for text")
# ________________________________________________________

filename_overleaf <- file.path(out.dir.tables, "overleaf_PLHIVcontribution_custom.rds")

dcens_custom <- copy(dcens)
dcens_custom[, AGEGROUP := split.agegroup(AGEYRS, breaks = c(15, 30, 50))]

if (file.exists(filename_overleaf) & !overwrite) {
    contrib_plhiv_custom <- readRDS(filename_overleaf)
} else {
    contrib_plhiv_custom <- dfiles_rds[MODEL == "run-gp-prevl" & ROUND %in% c(16, 19),
        {
            eval(expr_setup_ftp_all)
            get.weighted.average.p_predict(
                fit1 = readRDS(paths[idx.all]),
                fit2 = readRDS(paths[idx.ftp]),
                fit_all = readRDS(paths),
                round = unique(ROUND),
                expression_prereturn = {
                    # removed "SEX" here
                    by_cols1 <- c(dot.cols, "LOC", "AGEGROUP")
                    by_cols2 <- c(dot.cols, "LOC", "SEX", "AGEGROUP")

                    tmp <- merge(draws_all, dcens_custom, by = c("LOC", "SEX", "AGEYRS", "ROUND"))

                    tmp1 <- tmp[,
                        {
                            z <- joint * ELIGIBLE_SMOOTH
                            list(z = sum(z))
                        },
                        by = by_cols1
                    ][,
                        .(AGEGROUP = AGEGROUP, CONTRIBUTION = z / sum(z)),
                        by = c(dot.cols, "LOC")
                    ]
                    tmp1[, SEX := "Total"]

                    tmp2 <- tmp[,
                        {
                            z <- joint * ELIGIBLE_SMOOTH
                            list(z = sum(z))
                        },
                        by = by_cols2
                    ][,
                        .(AGEGROUP = AGEGROUP, SEX = SEX, CONTRIBUTION = z / sum(z)),
                        by = c(dot.cols, "LOC")
                    ]

                    rbind(tmp1, tmp2)[, quantile2(CONTRIBUTION), by = c("LOC", "SEX", "AGEGROUP")]
                }
            )
        },
        by = c("MODEL", "ROUND")
    ]
    saveRDS(object = contrib_plhiv_custom, file = filename_overleaf)
    # print statements for paper
}

if (make_tables) {
    tmp <- paper_statements_contributions_PLHIV_custom()
}

catn("Within gender contributions to PLHIV")
# __________________________________________

filename_rds <- file.path(out.dir.tables, "contribhiv_sex_byagegroup.rds")

dcens_custom <- copy(dcens)
dcens_custom[, AGEGROUP := split.agegroup(AGEYRS, breaks = c(15, 25, 40, 50))]

# dcens_very_custom <- copy(dcens)
# dcens_very_custom[LOC == "fishing", AGEGROUP := split.agegroup(AGEYRS, breaks=c(15, 25, 35, 50)),]
# dcens_very_custom[LOC == "inland",  AGEGROUP := split.agegroup(AGEYRS, breaks=c(15, 30, 40, 50)),]

if (file.exists(filename_rds) & !overwrite) {
    dsexcontrib_agegroup <- readRDS(filename_rds)
} else {
    dsexcontrib_agegroup <- dfiles_rds[MODEL == "run-gp-prevl",
        {
            eval(expr_setup_ftp_all)
            get.weighted.average.p_predict(
                fit1 = readRDS(paths[idx.all]),
                fit2 = readRDS(paths[idx.ftp]),
                fit_all = readRDS(paths),
                round = unique(ROUND),
                expression_prereturn = {
                    tmp <- merge(draws_all, dcens_custom, by = c("LOC", "SEX", "AGEYRS", "ROUND"))
                    tmp <- tmp[, .(z = sum(joint * ELIGIBLE_SMOOTH)), by = c(dot.cols, "LOC", "SEX", "AGEGROUP")]
                    tmp <- tmp[, .(
                        AGEGROUP = AGEGROUP,
                        CONTRIBUTION = z / sum(z)
                    ), by = c(dot.cols, "LOC", "SEX")]

                    tmp[, quantile2(CONTRIBUTION), by = c("SEX", "LOC", "AGEGROUP")]
                }
            )
        },
        by = c("MODEL", "ROUND")
    ]
    saveRDS(object = dsexcontrib_agegroup, file = filename_rds)
}

if (make_tables) {
    dsexcontrib_agegroup |> plot_quantiles(facet = LOC ~ ROUND, color = SEX_LAB, x = AGEGROUP)
}

catn("Get quantiles for contributions by agegroup")
# __________________________________________________

filename_rds <- file.path(out.dir.tables, "fullpop_allcontributions_byagegroup.rds")

if (file.exists(filename_rds) & !overwrite) {
    dcontrib_agegroup <- readRDS(filename_rds)
} else {
    dcontrib_agegroup <- dfiles_rds[,
        {
            eval(expr_setup_ftp_all)
            get.weighted.average.p_predict(
                fit1 = readRDS(paths[idx.all]),
                fit2 = readRDS(paths[idx.ftp]),
                fit_all = readRDS(paths),
                round = unique(ROUND),
                expression_prereturn = {
                    tmp <- merge(draws_all, dcens, by = c("LOC", "SEX", "AGEYRS", "ROUND"))
                    tmp <- tmp[, list(z = sum(joint * ELIGIBLE_SMOOTH)), by = c(dot.cols, "LOC", "SEX", "AGEGROUP")]
                    tmp <- tmp[, list(
                        AGEGROUP = AGEGROUP,
                        SEX = SEX,
                        CONTRIBUTION = z / sum(z)
                    ), by = c(dot.cols, "LOC")]
                    tmp_totals <- rbind(
                        tmp[, .(AGEGROUP = "Total", CONTRIBUTION = sum(CONTRIBUTION)), by = c(dot.cols, "LOC", "SEX")],
                        tmp[, .(AGEGROUP = "Total", SEX = "Total", CONTRIBUTION = sum(CONTRIBUTION)), by = c(dot.cols, "LOC")]
                    )
                    tmp <- rbind(tmp, tmp_totals)
                    tmp[, quantile2(CONTRIBUTION), by = c("SEX", "LOC", "AGEGROUP")]
                }
            )
        },
        by = c("MODEL", "ROUND")
    ]
    saveRDS(object = dcontrib_agegroup, file = filename_rds)

    dcontrib_agegroup[
        AGEGROUP == "Total" & ROUND %in% c(16, 19) & LOC == "inland" & MODEL %like% "supp-pop", prettify_cell(100 * M, 100 * CL, 100 * CU),
        by = c("SEX", "ROUND")
    ]
}

check_median_contr_approx1(dcontrib_agegroup)

if (make_plots) {
    .w <- 10
    .h <- 12

    p_contrib_prevl <- dcontrib_agegroup |> plot.agesex.contributions.by.roundcomm2(label = "run-gp-prevl")
    # p_contrib_supph <- dcontrib_agegroup |> plot.agesex.contributions.by.roundcomm2(label = "run-gp-supp-hiv")
    p_contrib_suppp <- dcontrib_agegroup |> plot.agesex.contributions.by.roundcomm2(label = "run-gp-supp-pop")

    .fnm <- function(lab) {
        paste("contrib_agegroupgender", lab, "byroundcomm.pdf", sep = "_")
    }

    ggsave2(p = p_contrib_prevl, file = .fnm("prevl"), LALA = out.dir.figures, .w, .h)
    # ggsave2(p = p_contrib_supph, file = .fnm("suppofhiv"), LALA = out.dir.figures, .w, .h)
    ggsave2(p = p_contrib_suppp, file = .fnm("suppofpop"), LALA = out.dir.figures, .w, .h)
}

if (make_tables) {
    t_contrib_supp <- tablify.agecontributions(dcontrib_agegroup, model = "run-gp-supp-pop")
    t_contrib_supp$ROUND_LAB <- labeller(ROUND_LAB = round_labs)(t_contrib_supp[, .(ROUND_LAB)])

    filename_tex <- file.path(out.dir.tables, "table_contrib_supppop.tex")
    write.to.tex(t_contrib_supp, file = filename_tex)
    filename <- "table_contrib_supppop.pdf"
    p <- table.to.plot(t_contrib_supp)
    ggsave2(p = p, file = filename, LALA = out.dir.tables, w = 22.5, h = 5.5)

    t_contrib_hiv <- tablify.agecontributions(dcontrib_agegroup, model = "run-gp-prevl")
    t_contrib_hiv$ROUND_LAB <- labeller(ROUND_LAB = round_labs)(t_contrib_hiv[, .(ROUND_LAB)])

    filename_tex <- file.path(out.dir.tables, "table_contrib_hivpop.tex")
    write.to.tex(t_contrib_hiv, file = filename_tex)
    filename <- "table_contrib_hivpop.pdf"
    p <- table.to.plot(t_contrib_hiv)
    ggsave2(p = p, file = filename, LALA = out.dir.tables, w = 22.5, h = 5.5)

    # dcontrib_agegroup |>
    #     subset(MODEL %like% 'supp-pop' & ROUND == 19) |>
    #     plot_quantiles(color=SEX_LAB, facet=ROUND ~ LOC_LAB , x=AGEGROUP)
}

if (make_tables) {
    tmp <- paper_statements_female_contributions_prevalence(dcontrib_agegroup)
}

catn("Get log-ratio for suppression among FTP and non-FTP")
# __________________________________________________________

# TODO: take blue-black plots, and take posterior sample ratios.
# Report the ratio (or log-ratio) by 5 years age group, and whether any significantly > 1 (or > 0).

filename_rds <- file.path(out.dir.tables, "posterior_ftp_logratio_quantiles.rds")
filename_rds2 <- file.path(out.dir.tables, "posterior_ftp_difference_quantiles.rds")

# lo(parts) - log(ftp)
if (file.exists(filename_rds) & !overwrite) {
    dlogratio <- readRDS(filename_rds)
} else {
    dlogratio <- dfiles_rds[,
        {
            eval(expr_setup_ftp_all)
            get.posterior.logratios.ftp(
                fit1=readRDS(paths[idx.all]),
                fit2=readRDS(paths[idx.ftp]),
                fit_all = readRDS(paths),
                round = unique(ROUND)
            )
        },
        by = c("MODEL", "ROUND")
    ]
    saveRDS(object = dlogratio, filename_rds)
}

if (file.exists(filename_rds2) & !overwrite) {
    d_diff <- readRDS(filename_rds2)
} else {
    d_diff <- dfiles_rds[,
        {
            eval(expr_setup_ftp_all)
            get.posterior.diff.ftp(
                fit1=readRDS(paths[idx.all]),
                fit2=readRDS(paths[idx.ftp]),
                fit_all = readRDS(paths)
            )
        },
        by = c("MODEL", "ROUND")
    ]
    saveRDS(object = d_diff, filename_rds2)
}


if (make_plots) {
    # only need to change sligtly, nice!
    .w <- 10
    .h <- 12
    .out <- out.dir.figures

    MODELS <- c("run-gp-prevl", "run-gp-supp-hiv", "run-gp-supp-pop")
    cols <- c("M", "CL", "CU", "IL", "IU")
    dratio <- copy(dlogratio)
    dratio[, (cols) := lapply(.SD, exp), .SDcols = cols]

    p_diff <- lapply(MODELS, plot.diff.ftpvsnon, DT = d_diff)
    p_logp <- lapply(MODELS, plot.logratio.ftpvsnon, DT = dlogratio, log = TRUE)
    names(p_logp) <- names(p_diff) <- MODELS

    .fnm <- function(lab) {
        paste("logratio", lab, "ftpvsnnon_byroundcomm.pdf", sep = "_")
    }
    with(p_logp, {
        ggsave2(p = `run-gp-prevl`, file = .fnm("prevl"), LALA = .out, .w, .h)
        ggsave2(p = `run-gp-supp-hiv`, file = .fnm("suppofhiv"), LALA = .out, .w, .h)
        ggsave2(p = `run-gp-supp-hiv`, file = .fnm("suppofpop"), LALA = .out, .w, .h)
    })

    .fnm2 <- function(lab) {
        paste("diff", lab, "ftpvsnnon_byroundcomm.pdf", sep = "_")
    }
    with(p_diff, {
        ggsave2(p = `run-gp-prevl`, file = .fnm2("prevl"), LALA = .out, .w, .h)
        ggsave2(p = `run-gp-supp-hiv`, file = .fnm2("suppofhiv"), LALA = .out, .w, .h)
        ggsave2(p = `run-gp-supp-hiv`, file = .fnm2("suppofpop"), LALA = .out, .w, .h)
    })

    rm(.out, .w, .h)
}

if (make_tables & 0) { # age groups for which CrI does not include 0
    dlogratio[(CL > 0 | CU < 0) & MODEL == "run-gp-supp-pop",
        {
            if (.N > 0) {
                list(min = min(AGEYRS), max = max(AGEYRS))
            }
        },
        by = c("MODEL", "ROUND", "SEX", "LOC")
    ]
}



catn("Increases in suppression relative to round 16")
# ____________________________________________________

filename_rds <- file.path(out.dir.tables, "posterior_suppressionincrease_vsround16.rds")
filename_rds2 <- file.path(out.dir.tables, "posterior_suppressionincrease_diff_vsround16.rds")

if (file.exists(filename_rds2) & !overwrite) {
    dincreasessupp <- readRDS(filename_rds)
    dincreasessupp_diff <- readRDS(filename_rds2)
} else {
    # doesn't work cause there are 2
    draws16 <- dfiles_rds[MODEL == "run-gp-supp-hiv" & ROUND == 16, {
        paths <- file.path(D, F)
        idx.all <- which(FTP == FALSE)
        idx.ftp <- which(FTP == TRUE)
        cat(paths[idx.all], "\n")
        get.weighted.average.p_predict(
            fit1 = readRDS(paths[idx.all]),
            fit2 = readRDS(paths[idx.ftp]),
            fit_all = readRDS(paths),
            round = unique(ROUND),
            expression_prereturn = draws_all
        )
    }]
    draws16[, `:=`(joint16 = joint, joint = NULL, parts = NULL, ftp = NULL, ROUND = NULL, PARTRATE = NULL)]

    dincreasessupp <- dfiles_rds[MODEL == "run-gp-supp-hiv" & ROUND >= 17,
        {
            paths <- file.path(D, F)
            idx.all <- which(FTP == FALSE)
            idx.ftp <- which(FTP == TRUE)
            cat(paths[idx.all], "\n")
            get.weighted.average.p_predict(
                fit1 = readRDS(paths[idx.all]),
                fit2 = readRDS(paths[idx.ftp]),
                fit_all = readRDS(paths),
                round = unique(ROUND),
                expression_prereturn = {
                    draws_all <- merge(draws_all, draws16, by = c(demo.cols, dot.cols))
                    if ("joint16.x" %in% names(draws_all)) {
                        draws_all[, `:=`(joint16 = joint16.x, joint16.x = NULL, joint16.y = NULL)]
                    }
                    draws_all[, `:=`(joint = joint / joint16)]
                    return(draws_all[, quantile2(joint), by = c("SEX", "LOC", "AGEYRS")])
                }
            )
        },
        by = c("MODEL", "ROUND")
    ]

    saveRDS(object = dincreasessupp, filename_rds)

    dincreasessupp_diff <- dfiles_rds[MODEL == "run-gp-supp-hiv" & ROUND >= 17,
        {
            paths <- file.path(D, F)
            idx.all <- which(FTP == FALSE)
            idx.ftp <- which(FTP == TRUE)
            cat(paths[idx.all], "\n")
            get.weighted.average.p_predict(
                fit1 = readRDS(paths[idx.all]),
                fit2 = readRDS(paths[idx.ftp]),
                fit_all = readRDS(paths),
                round = unique(ROUND),
                expression_prereturn = {
                    draws_all <- merge(draws_all, draws16, by = c(demo.cols, dot.cols))
                    if ("joint16.x" %in% names(draws_all)) {
                        draws_all[, `:=`(joint16 = joint16.x, joint16.x = NULL, joint16.y = NULL)]
                    }
                    draws_all[, `:=`(joint = joint - joint16)]
                    return(draws_all[, quantile2(joint), by = c("SEX", "LOC", "AGEYRS")])
                }
            )
        },
        by = c("MODEL", "ROUND")
    ]

    saveRDS(object = dincreasessupp_diff, filename_rds2)
    rm(draws16)
}


if (make_plots) {
    p_supphiv_logprob16_ratio <- plot.relative.suppression.vs.round16.ratio(dincreasessupp)
    p_supphiv_logprob16_diff <- plot.relative.suppression.vs.round16.diff(dincreasessupp_diff)

    .w <- 9
    .h <- 8

    filename <- paste0("fit_supphiv_logprob_vs_baseline_by_locgenderage.pdf")
    ggsave2(p = p_supphiv_logprob16_ratio, file = filename, LALA = out.dir.figures, w = .w, h = .h)
    filename <- paste0("fit_supphiv_logprob_vs_baseline_by_locgenderage_diff.pdf")
    ggsave2(p = p_supphiv_logprob16_diff, file = filename, LALA = out.dir.figures, w = .w, h = .h)

    rm(.w, .h)
}

if (make_tables) {
    tmp <- dincreasessupp[ROUND == 19 & M > 5]
    prettify_labels(tmp)
    sprintf("Posterior medians for ratios were larger than 5 for:\n")
    tmp[,
        {
            sprintf("%s %s aged %s to %s\n", unique(LOC_LAB), unique(SEX_LAB), min(AGEYRS), max(AGEYRS)) |> cat()
            NULL
        },
        by = c("LOC", "SEX")
    ]

    tmp <- dincreasessupp[ROUND == 19 & CL > 3]
    prettify_labels(tmp)
    sprintf("Posterior medians for ratios were larger than 5 for:\n")
    tmp[,
        {
            sprintf("%s %s aged %s to %s\n", unique(LOC_LAB), unique(SEX_LAB), min(AGEYRS), max(AGEYRS)) |> cat()
            NULL
        },
        by = c("LOC", "SEX")
    ]

    # filename_overleaf <- file.path(out.dir.tables, "overleaf_")
    # saveRDS(object = tab, file = filename_overleaf)
}


############################################################
catn("=== M/F PLHIV suppression ratio custom agegroups ===")
############################################################

# dcens_custom <- copy(dcens)
# dcens_custom[, AGEGROUP := split.agegroup(AGEYRS, breaks = c(15, 25, 40, 50))]

filename_rds <- file.path(out.dir.tables, "posterior_mf_ratios_custom.rds")
filename_rds2 <- file.path(out.dir.tables, "posterior_mf_ratios_custom_pvals.rds")

if (file.exists(filename_rds) & !overwrite) {
    dmf_ratios <- readRDS(filename_rds)
} else {
    p_values <- data.table(
        LOC = as.character(),
        ROUND = as.character(),
        P = as.numeric(),
        AGEGROUP = as.character()
    )
    dmf_ratios <- dfiles_rds[,
        {
            paths <- file.path(D, F)
            .check <- function(x) {
                stopifnot(length(x) == 1 | args$shared.hyper)
                return(x)
            }
            # eval(expr_setup_ftp_all)
            eval(expr_setup_hivprev_supphiv_ftpall )
            cat("Round ", unique(ROUND), "\n")
            # load hiv prevalences to compute N HIV positive in age groups
            draws_prev <- get.weighted.average.p_predict(
                fit1 = readRDS(paths[idx.hivprev.all]),
                fit2 = readRDS(paths[idx.hivprev.ftp]),
                fit_all = readRDS(paths[idx.hivprev]),
                round = unique(ROUND),
                expression_prereturn = {
                    # find composition of prevalnce by age group
                    by_cols <- c(dot.cols, "LOC", "SEX", "AGEGROUP")
                    tmp <- merge(draws_all, dcens, by = c("LOC", "SEX", "AGEYRS", "ROUND"))
                    tmp <- tmp[, N_HIV := joint * ELIGIBLE_SMOOTH] |>
                        subset(select = c(dot.cols, "LOC", "SEX", "AGEGROUP", "AGEYRS", "N_HIV"))
                    return(tmp)
                }
            )
            cat("prevalence done - ")
            draws_supp <- get.weighted.average.p_predict(
                fit1 = readRDS(paths[idx.supphiv.all]),
                fit2 = readRDS(paths[idx.supphiv.ftp]),
                fit_all = readRDS(paths[idx.supphiv]),
                round = unique(ROUND),
                expression_prereturn = {
                    draws_all
                }
            )
            cat("suppression done\n")
            dot.cols <- c(".chain", ".iteration", ".draw")
            tmp <- merge(draws_supp, draws_prev, by = c(dot.cols, "LOC", "SEX", "AGEYRS"))
            .take.supp.ratio <- function(DT, by_cols) {
                by_cols_noage <- setdiff(by_cols, "AGEGROUP")
                by_cols_nosex <- setdiff(by_cols, "SEX")
                .formula <- as.formula(paste0(".draw +", paste(by_cols_nosex, collapse = " + "), " ~ SEX"))
                tmp1 <- tmp[, .(S = sum(joint * proportions(N_HIV))), by = c(dot.cols, by_cols)] |>
                    dcast.data.table(.formula, value.var = "S", drop = TRUE)
                tmp1[, `:=`(
                    RATIO_MF_SUPP = M / F,
                    RATIO_MF_VIR = (1 - M) / (1 - F),
                    DIFF_MF_SUPP = F - M
                )]
                out <- rbind(
                    cbind(tmp1[, quantile2(RATIO_MF_SUPP), by = by_cols_nosex], TYPE = "SUP"),
                    cbind(tmp1[, quantile2(RATIO_MF_VIR), by = by_cols_nosex], TYPE = "VIR"),
                    cbind(tmp1[, quantile2(DIFF_MF_SUPP), by = by_cols_nosex], TYPE = "SUP-DIFF")
                )
                return(list(Q=out, pv=tmp1[, .(P=mean(DIFF_MF_SUPP > 0)), by=by_cols_nosex]))
            }
            lala <- .take.supp.ratio(tmp, by_cols = c("LOC", "SEX", "AGEGROUP"))
            lala2 <- .take.supp.ratio(tmp, by_cols = c("LOC", "SEX"))

            p_values  <<- rbind(p_values, lala2$pv[, `:=` (ROUND=ROUND, AGEGROUP="Total")])
            list( lala$Q, lala2$Q[, AGEGROUP := "Total"]) |> 
                rbindlist(use.names = TRUE)
        },
        by = c("ROUND")
    ]
    saveRDS(object = dmf_ratios, filename_rds)
    p_values$type <- "Suppression difference is larger than 0"
    saveRDS(object = p_values, filename_rds2)
}

if (make_tables) {
    # paper_statements_malefemaleratio_suppression(DT = dmf_ratios, reverse = FALSE)
    tmp <- paper_statements_malefemaleratio_suppression2()

    # paper_statements_malefemaleratio_suppression()

    # dmf_ratios[TYPE == "VIR" & AGEGROUP != "Total" & ROUND == 19, ] |>
    #     plot_quantiles(x=AGEGROUP, facet=.~LOC_LAB, color=LOC_LAB)
    # dmf_ratios[TYPE == "VIR" & AGEGROUP != "Total" & ROUND == 19, AGEGROUP[which.max(M)], by=LOC]

    # djoint |>
    #     subset(MODEL %like% 'prevl') |>
    #     plot_quantiles(x=AGEYRS, facet=ROUND~LOC_LAB, color=SEX_LAB)
}

if (make_tables) {
    ## age-aggregated HIV prevalence by sex, round, loc
    # tab <- make.table.Nhivpositive (DT=joint_ageagrr_list$round_totals, DC=dcens)
    # filename_overleaf <- file.path(out.dir.tables, "overleaf_ageaggr_hivprev.rds")

    dsupp_amonghiv <- subset(dsupp_agegroup, AGEGROUP == "Total" & SEX != "Total") |>
        prettify_labels()

    tab_mf_diffs <- make.table.malefemale.diff.viraemia(DT = dmf_ratios)

    .dict <- dict_table_names$percent_reduction
    tab_eligible <- dcens[, 
        .(N_ELIGIBLE = comma(sum(ELIGIBLE))),
        by = c("ROUND", "LOC", "SEX")] |>
        prettify_labels() |>
        remove.nonpretty()

    tab_pvir <- djoint_agegroup[ MODEL %like% 'supp-pop' & AGEGROUP == "Total" & SEX != "Total",  .(
        LOC, SEX, ROUND,
        VIR_P = prettify_cell(M*100, CL*100, CU*100, percent = TRUE)
    )] |> prettify_labels()
    set(tab_pvir, j=c("LOC", "SEX", "ROUND"), value=NULL)

    tab_hiv <- tablify.posterior.Nunsuppressed(joint_ageagrr_list, 
        CELLname = "HIV",
        model = "run-gp-prevl")
    tab_unsupp <- tablify.posterior.Nunsuppressed(joint_ageagrr_list,
        CELLname = "UNSUPP",
        model = "run-gp-supp-pop")
    dsupp_amonghiv <- dsupp_amonghiv[, .(
        ROUND_LAB, LOC_LAB, SEX_LAB, 
        CELL_SUPPHIV = prettify_cell(M * 100, CL * 100, CU * 100, percent = TRUE)
    )] 

    by_cols <- c("LOC_LAB", "SEX_LAB", "ROUND_LAB")
    tab_merge <- merge(tab_eligible, tab_hiv) |>
        merge(tab_unsupp, by=by_cols) |>
        merge(dsupp_amonghiv, by=by_cols) |>
        merge(tab_pvir, by=by_cols) |>
        merge(tab_mf_diffs, by=by_cols)
    tab_merge[, SEX_LAB := unname(sex_dictionary2[SEX_LAB])]
    tab_merge <- delete.repeated.table.values(tab_merge) |>
        setnames(names(.dict), unname(.dict), skip_absent = TRUE)
    setcolorder(tab_merge, unname(.dict))

    tab_merge <- subset(tab_merge, select = - CELL_HIV_P)
    .sd <- names(tab_merge)[-(1:2)]
    tab_merge[, (.sd) := lapply(.SD, function(x) { x[x == "" | is.na(x)] <- "--"; x}), .SDcols = .sd]
    colnames(tab_merge)

    if (interactive()) {
        write.to.googlesheets(tab_merge, sheet = "Table2")
    }
    filename_table <- file.path(out.dir.tables, "table_reductionHIVandUNSUPP.rds")
    saveRDS(tab_merge, file = filename_table)
    # view_xl(tab_merge)
}



###################
catn("Assess fits")
###################

# necessary for later...
source(file.path(gitdir.functions,"phsc_vl_cmdstan_helpers.R"))
group_codes <- expand.grid( SEX=0:1, LOC=0:1)
setDT(group_codes)
group_codes[, `:=`( 
    SEX_LABEL=stan_dicts$INTtoSEX[as.character(SEX)], 
    LOC_LABEL=stan_dicts$INTtoLOC[as.character(LOC)] )]

# Locate rds files for model fits and data
dirs <- list.files(args$out.dir.exact, pattern="^run-gp", full.names = TRUE) 
tmp <-  list.files(dirs , pattern = "rds$", full.names = TRUE)
dfits <- data.table(
    F=tmp,
    RDS = fifelse( basename(tmp) %like% "standata.rds", yes ="data" , no="fit"),
    ROUND = as.integer( basename(tmp) |> stringr::str_extract("(?<=round)\\d+")),
    MODEL = basename(dirname(tmp))
)

gg_list <- list()
out <- dfits[, {
    sprintf( "Model: %s, round: %s\n", MODEL[1], ROUND[1]) |> cat()
    plot_single_posterior_fit(
        fit_rds = F[RDS == 'fit'],
        standata_rds = F[RDS == 'data'],
        model=MODEL[1],
        round=ROUND[1]
    )
}, by = c("MODEL", "ROUND")]

cmds <- lapply(
    c('prevl', 'supp-hiv', 'supp-pop'),
    aggregate_posterior_fits,
    filename_fmt = "suppfig_all_posteriors_vs_data_%s.pdf"
)

# lapply( cmds, cmd2overleaf)

# posterior predictive

gg_list_pp <- list()
out <- lapply( 
    c("run-gp-prevl","run-gp-supp-hiv","run-gp-supp-pop"),
    plot_single_posterior_pred,
    DT=dfits
)

gg_list_pp <- lapply(
    gg_list_pp, 
    function(p) 
        p + theme(legend.margin = margin(r=0,l=-10, t = 0, b=0, unit = "pt")) 
)

p_pp <- ggarrange(
    plotlist = gg_list_pp,
    ncol=1, nrow=3,
    common.legend = TRUE, legend='bottom',
    labels = "auto",
    font.label = list(size = 10, color = "black", face = "bold", family = NULL)
) + nm_reqs 
p_pp

filename <- paste0("suppfig_all_ppchecks.pdf")
cmd <- ggsave2(p = p_pp, file = filename, LALA = out.dir.figures, w = 18 , h = 24, u="cm")
# system(cmd)
# cmd2overleaf(cmd)

#####################
catn("End of script")
#####################
