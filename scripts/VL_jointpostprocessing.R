#!/bin/Rscript

self_relative_path <- "scripts/VL_jointpostprocessing.R"
########################
cat("\nStart of:", self_relative_path,"\n")
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
} |> suppressPackageStartupMessages()

################
#    PATHS     #
################

# NOTE: contributions by age/group do not make sense in the "among hiv", if we merge by dcens.
# instead, we should be merging by the N of hiv+

if( interactive() ){
    gitdir <- here::here()
}else{
    cmd <- commandArgs()
    cmd <- cmd[cmd %like% "file"]
    gitdir <- gsub(paste0(".*--file=(.*)/", self_relative_path, "$"), "\\1", cmd)
}
source(file.path(gitdir, "R/paths.R"))

file.exists(
    path.hivstatusvl.r1520,
    path.census.eligible
) |>
    all() |>
    stopifnot()

# command line options, stored in args. Then subset
opts_vec <- c(
    "viremic.viral.load",
    "detectable.viral.load",
    "out.dir.prefix",
    "out.dir.exact",
    "round",
    "jobname"
)
args <- args[names(args) %in% opts_vec]
print(args)
# testing
if(interactive()){
    args$jobname <- 'vl_1000_firstpart'
}
source(file.path(gitdir.functions, "plotting_main_figures.R"))
source(file.path(gitdir.functions, "postprocessing_helpers.R"))
source(file.path(gitdir.functions, "postprocessing_tables.R"))
source(file.path(gitdir.functions, "phsc_vl_helpers.R"))
source(file.path(gitdir.functions, "paper_statements.R"))
naturemed_reqs()

overwrite <- !interactive()
make_plots <- TRUE
make_tables <- TRUE

with(args, {
    VL_DETECTABLE <<- detectable.viral.load
    VIREMIC_VIRAL_LOAD <<- viremic.viral.load
    # need to check that we have both samples for "" and "_firstpart"
    out.dir <<- gsub("_firstpart$", "_joint", out.dir.exact)
    if(is.na(out.dir.exact)){
        out.dir <<- file.path(
            out.dir.prefix, 
            gsub('_firstpart$','_joint',jobname))
    }
})
stopifnot("out.dir must end in _joint"= out.dir %like% '_joint$')
indir.ftp <- gsub( '_joint', "_firstpart", out.dir)
indir.all <- gsub( '_joint', "", out.dir)
stopifnot(
    "did not find 2 directories necessary to run joint analysis"=all(dir.exists(c(indir.ftp, indir.all)))
)
out.dir.figures <- file.path(out.dir, "figures")
out.dir.tables <- file.path(out.dir, "tables")
dir.create(out.dir.tables) |> suppressWarnings()
dir.create(out.dir.figures) |> suppressWarnings()

####################
catn("=== MAIN ===")
####################

# load first participation rates
dfirst_prop <- get.first.participant.rates()

# get census eligible
dcens <- get.census.eligible() |>
    setnames(c("AGE_LABEL", "SEX_LABEL", "LOC_LABEL"), c("AGEYRS", "SEX", "LOC"))
dcens[, AGEGROUP := split.agegroup(AGEYRS)]

# load number of census eligible individuals (.50 too rough)
dpartrates <- readRDS(path.participation.rates) |>
    subset(select = c("ROUND", "FC", "SEX", "AGEYRS", "PARTRATE_SMOOTH.25")) |>
    setnames(c("FC", "PARTRATE_SMOOTH.25"), c("LOC", "PARTRATE"))

####################################################
catn("=== Compare model fits among FTP and ALL ===")
####################################################

dfiles_rda <- get.output.paths.ftp.and.all(regex=".rda$")
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

# particular focus on suppression among HIV positivss
if (make_plots) {
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
                LOC_LAB = community_dictionary$longest,
                ROUND_LAB = round_labs
                )
        ) +
        scale_x_continuous(expand = c(0, 0), limits = c(15, 49), breaks = seq(15, 50, by = 5)) +
        scale_y_percentage +
        scale_color_manual(values = palettes$sex, labels=sex_dictionary2) +
        scale_fill_manual(values = palettes$sex, labels=sex_dictionary2) +
        my_labs(y = "Prevalence of suppression among HIV positives", fill = "Gender", linetype="") +
        NULL
    p2 <- p + geom_ribbon(alpha = .1, aes(ymin = CL, ymax = CU, fill = SEX_LAB), color = NA)

    filenames <- paste0("fit_suppofhiv_compare_ftpvsall", c("", "CrIs"), "_r1619.pdf")

    cmd <- ggsave2(p = p, file = filenames[1], LALA = out.dir.figures, w = 16, h = 14)
    cmd <- ggsave2(p = p2, file = filenames[2], LALA = out.dir.figures, w = 16, h = 14)
}


##################################################
catn("=== Get posterior draws from rds files ===")
##################################################

dfiles_rds <- get.output.paths.ftp.and.all("round1[0-9].rds$|220729.rds$")
stopifnot(dfiles_rds[, .N, by = "F"][, all(N == 2)])
stopifnot(dfiles_rds[, .N, by = "MODEL"][, all(N == 8)])

catn("Load all previous results and plot")
# _________________________________________

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
    p1_vir <- plot.comparison.prevalence.fishinginland.oneround(DT=djoint, model = "run-gp-supp-pop", round = 19)
    p1_supp <- plot.comparison.prevalence.fishinginland.oneround(DT=djoint, model = "run-gp-supp-hiv", round = 19, ylim = 1)
    p1_hiv <- plot.comparison.prevalence.fishinginland.oneround(DT=djoint, model = "run-gp-prevl", round = 19)
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
    p1 <- plot_propofpop_of_viraemic_byagesex_stratbycommround(djoint, colorby="ROUND", cri=TRUE)
    p2 <- plot_propofpop_of_viraemic_byagesex_stratbycommround(djoint, colorby="ROUND", cri=FALSE)
    filename <- paste0("propofpop_of_viraemic_byagesex_stratbycommround", c("_cri", ""), ".pdf")
    cmd <- ggsave2(p = p1, file = filename[1], LALA = out.dir.figures, w=.w, h=.w, u="cm")
    cmd <- ggsave2(p = p2, file = filename[2], LALA = out.dir.figures, w=.w, h=.w, u="cm")

    rm(.w, .h)
}


catn("Get quantiles for population prevalences by agegroup")
# ___________________________________________________________

filename_rds <- file.path(out.dir.tables, "posterior_quantiles_agegroups.rds")

if (file.exists(filename_rds) & !overwrite) {
    djoint_agegroup <- readRDS(filename_rds)
} else {
    djoint_agegroup <- dfiles_rds[MODEL != "run-gp-supp-hiv",
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
        },
        by = c("MODEL", "ROUND")
    ]

    djoint_agegroup[is.na(AGEGROUP), AGEGROUP := "Total"]
    djoint_agegroup[is.na(SEX), SEX := "Total"]
    saveRDS(object = djoint_agegroup, filename_rds)
}

if (make_tables) {
    tmp <- paper_statements_female_prevalence(djoint_agegroup)
    tmp <- paper_statements_prevalence_viraemia(djoint_agegroup)

}

catn("Get quantiles for population prevalences aggregated over age")
# ___________________________________________________________

filename_rds <- file.path(out.dir.tables, "posterior_quantiles_ageaggregate.rds")

if (file.exists(filename_rds) & !overwrite) {
    joint_ageagrr_list <- readRDS(filename_rds)
} else {
    djoint_ageaggr <- dfiles_rds[MODEL != "run-gp-supp-hiv",
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
                        list(joint = sum(joint * ELIGIBLE_SMOOTH)),
                        by = c(dot.cols, "LOC", "SEX")
                    ]
                }
            )
        },
        by = c("MODEL", "ROUND")
    ]
    .f <- function(old, new) {
        (old - new) / old
    }

    tmp <- dcast(djoint_ageaggr,
        MODEL + .chain + .iteration + .draw + LOC + SEX ~ paste0("R", ROUND),
        value.var = "joint"
    )[
        ,
        list(
            MODEL = MODEL, LOC = LOC, SEX = SEX,
            redR17 = .f(R16, R17),
            redR18 = .f(R16, R18),
            redR19 = .f(R16, R19),
            redR1718 = .f(R17, R18),
            redR1819 = .f(R18, R19)
        )
    ] |>
        melt(
            id.vars = c("MODEL", "LOC", "SEX"),
            measure.vars = c("redR17", "redR18", "redR19", "redR1718", "redR1819")
        )

    joint_ageagrr_list <- list(
        percent_reduction = tmp[, quantile2(value), by = c("MODEL", "SEX", "LOC", "variable")],
        round_totals = djoint_ageaggr[, quantile2(joint), by = c("MODEL", "SEX", "LOC", "ROUND")]
    )
    saveRDS(object = joint_ageagrr_list, filename_rds)
    rm(djoint_ageaggr)
}

if (make_tables) {
    # age-aggregated HIV prevalence by sex, round, loc
    dcens_aggr <- dcens[, lapply(.SD, sum), .SDcols = c("ELIGIBLE_SMOOTH", "ELIGIBLE"), by = c("ROUND", "LOC", "SEX")]
    tab <- merge(
        joint_ageagrr_list$round_totals[MODEL == "run-gp-prevl"],
        dcens_aggr,
        by = c("ROUND", "LOC", "SEX")
    ) |> remove.ILIU()
    cols <- c("CL", "M", "CU")
    tab[, (cols) := lapply(.SD, function(x) x / ELIGIBLE), .SDcols = cols]
    tab[, CELL := prettify_cell(M * 100, CL * 100, CU * 100, percent = TRUE)]
    tab[ROUND == 19, {
        sprintf(
            "
            In round 19, %s of men and %s of women were estimated to be HIV positive in fishing communities,
            compared to %s and %s among men and women in inland, respectively.",
            CELL[LOC == "fishing" & SEX == "M"], CELL[LOC == "fishing" & SEX == "F"],
            CELL[LOC == "inland" & SEX == "M"], CELL[LOC == "inland" & SEX == "F"]
        ) |> cat()
    }]
    tab <- tab[, .(ROUND, LOC, SEX, CELL)]
    filename_overleaf <- file.path(out.dir.tables, "overleaf_ageaggr_hivprev.rds")
    saveRDS(object = tab, file = filename_overleaf)


    # age-aggregated # of unsuppressed by sex, and location
    tab2 <- tablify.posterior.Nunsuppressed(joint_ageagrr_list)

    filename_tex <- file.path(out.dir.tables, "table_aggregatedNunsuppressed.tex")
    write.to.tex(tab2, file = filename_tex)
    filename <- "table_aggregatedNunsuppressed.pdf"
    p <- table.to.plot(tab2)
    ggsave2(p = p, file = filename, LALA = out.dir.tables, w = 8, h = 5.5)
}

catn("Get quantiles for suppression levels in agegroups")

filename_rds <- file.path(out.dir.tables, "posterior_quantiles_suppression_agegroup.rds")

if( file.exists(filename_rds) & ! overwrite){
    dsupp_agegroup <- readRDS(filename_rds)
}else{
    dsupp_agegroup <- dfiles_rds[,{
        paths <- file.path(D,F)
        .check <- function(x) { stopifnot(length(x) == 1); return(x)}
        idx.hivprev.all <- which(FTP==FALSE & MODEL=="run-gp-prevl") |> .check()
        idx.hivprev.ftp <- which(FTP==TRUE & MODEL=="run-gp-prevl") |> .check()
        idx.supphiv.all <- which(FTP==FALSE & MODEL=="run-gp-supp-hiv") |> .check()
        idx.supphiv.ftp <- which(FTP==TRUE & MODEL=="run-gp-supp-hiv") |> .check()
        cat("Round ", unique(ROUND), "\n")
        draws_prev <- get.weighted.average.p_predict(
            readRDS(paths[idx.hivprev.all]),
            readRDS(paths[idx.hivprev.ftp]),
            round = unique(ROUND),
            expression_prereturn = {
                # find composition of prevalnce by age group
                by_cols <- c(dot.cols, "LOC", "SEX", "AGEGROUP")
                tmp <- merge(draws_all, dcens, by = c("LOC", "SEX", "AGEYRS", "ROUND"))
                tmp <- tmp[, N_HIV := joint * ELIGIBLE_SMOOTH] |> 
                    subset(select=c(dot.cols, "LOC", "SEX", "AGEGROUP", "AGEYRS", "N_HIV"))
                # tmp[, P_HIV_AGE := proportions(N_HIV), by=c(dot.cols, 'LOC', "SEX", "AGEGROUP")]
                # tmp[, P_HIV_SEX := proportions(N_HIV), by=c(dot.cols, 'LOC', "SEX")]
                # tmp[, P_HIV_LOC := proportions(N_HIV), by=c(dot.cols, "LOC")]
                return(tmp)
            }
        )
        cat("prevalence done\n")
        draws_supp <- get.weighted.average.p_predict(
            readRDS(paths[idx.supphiv.all]),
            readRDS(paths[idx.supphiv.ftp]),
            round = unique(ROUND),
            expression_prereturn = {draws_all}
        )
        cat("suppression done\n")
        dot.cols <- c(".chain", ".iteration", ".draw")
        tmp <- merge(draws_supp, draws_prev, by=c(dot.cols, "LOC", "SEX", "AGEYRS"))
        .aggr.hiv.prev <- function(DT, by_cols){
            DT[, .(S=sum(joint*proportions(N_HIV))) , by=c(dot.cols, by_cols)][, quantile2(S), by=by_cols]
        }
        list(
            .aggr.hiv.prev(tmp, by_cols=c("LOC", 'SEX', 'AGEGROUP')),
            .aggr.hiv.prev(tmp, by_cols=c("LOC", 'SEX'))[, AGEGROUP := "Total"],
            .aggr.hiv.prev(tmp, by_cols=c("LOC"))[, `:=` (SEX="Total",AGEGROUP="Total")]
        ) |> rbindlist(use.names=TRUE)
    } , by=c("ROUND")]
    saveRDS(object = dsupp_agegroup, filename_rds)
}

if(make_plots & FALSE){
    .w <- 14; .h <- 14
    p_list <- lapply(c(16:19), plot.prevalence.by.age.group, DT=dsupp_agegroup)
    filenames <- paste0("plot_suppofhiv_by_agegroup_round", c(16:19), ".pdf")
    for(i in 1:length(p_list)){
        # ggsave2(p=p_list[[i]], file=filenames[i], LALA=out.dir.figures, .w, .h)
    }
}

catn("Get quantiles for contributions by age")
# _____________________________________________

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

if ( make_plots ) {
    .w <- 10; .h <- 12

    dcontrib[, .N,by=c('LOC', 'MODEL', "SEX", "AGEYRS", "ROUND")]

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


catn("Get average age for PLHIV & viraemic population")
# _____________________________________________________

filename_rds <- file.path(out.dir.tables, "mean_ages_plhiv_viraemic.rds")

if (file.exists(filename_rds) & !overwrite) {
    dmeanage <- readRDS(filename_rds)
} else {
    dmeanage <- dfiles_rds[MODEL %in% c('run-gp-prevl', 'run-gp-supp-pop'), {
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
                tmp1 <- tmp[j=.(
                    AGEMEAN=Hmisc::wtd.mean(AGEYRS, joint * ELIGIBLE_SMOOTH),
                    AGESTD=Hmisc::wtd.var(AGEYRS, joint * ELIGIBLE_SMOOTH) |> sqrt()
                ),
                by=c(dot.cols, 'LOC', "ROUND","SEX")] |>
                    melt.data.table(
                        measure.vars = c("AGEMEAN", "AGESTD"),
                        variable.name = "TYPE",
                        value.name = "value"
                    )
                tmp1[j = quantile2(value), by=c("SEX", "LOC", "TYPE") ]
            }
        )
    },
    by = c("MODEL", "ROUND")
    ]
    saveRDS(object = dmeanage, file = filename_rds)
}

if ( make_tables){
    tmp <- rbind(
        paper_statements_meanage_population(DT=dmeanage, label='run-gp-supp-pop', type="AGEMEAN"),
        paper_statements_meanage_population(DT=dmeanage, label='run-gp-supp-pop', type="AGESTD")
    )
}

catn("Other contribution to viraemia quantiles for text")
# ________________________________________________________

filename_overleaf <- file.path(out.dir.tables, "overleaf_viraemiacontribution_custom.rds")

dcens_custom <- copy(dcens)
dcens_custom[, AGEGROUP := split.agegroup(AGEYRS, breaks = c(15, 25, 40, 50))]

if (file.exists(filename_overleaf) & !overwrite) {
    contrib_viraemia_custom <- readRDS(filename_overleaf)
} else {
    contrib_viraemia_custom <- dfiles_rds[MODEL == "run-gp-supp-pop",
        {
            dcens
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

    # print statements for paper
    paper_statements_contributions_viraemia_round(round=16)
    paper_statements_contributions_viraemia_round(round=19)

    # paper_statements_contributions_viraemia_round(round=16, agegroup="15-24")
    # paper_statements_contributions_viraemia_round(round=19, agegroup="15-24")
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
            dcens
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

tmp <- paper_statements_contributions_PLHIV_custom()

catn("Get quantiles for contributions by agegroup")
# __________________________________________________

filename_rds <- file.path(out.dir.tables, "fullpop_allcontributions_byagegroup.rds")

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
                    tmp <- tmp[, list(z = sum(joint * ELIGIBLE_SMOOTH)), by = c(dot.cols, "LOC", "SEX", "AGEGROUP") ]
                    tmp <- tmp[, list(
                        AGEGROUP = AGEGROUP,
                        SEX = SEX,
                        CONTRIBUTION = z / sum(z)
                    ), by = c(dot.cols, "LOC")]
                    tmp_totals <- rbind(
                        tmp[, .( AGEGROUP = "Total", CONTRIBUTION = sum(CONTRIBUTION)), by = c(dot.cols, "LOC", "SEX")],
                        tmp[, .( AGEGROUP = "Total", SEX="Total", CONTRIBUTION = sum(CONTRIBUTION)), by = c(dot.cols, "LOC")]
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
        AGEGROUP == "Total" & ROUND %in% c(16, 19) & LOC == 'inland' & MODEL %like% 'supp-pop', prettify_cell(100*M, 100*CL, 100*CU), by=c('SEX', "ROUND")]
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
            stopifnot(.N == 2)
            paths <- file.path(D, F)
            idx.all <- which(FTP == FALSE)
            idx.ftp <- which(FTP == TRUE)
            cat(paths[idx.all], "\n")
            get.posterior.logratios.ftp(
                readRDS(paths[idx.all]),
                readRDS(paths[idx.ftp]),
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
            stopifnot(.N == 2)
            paths <- file.path(D, F)
            idx.all <- which(FTP == FALSE)
            idx.ftp <- which(FTP == TRUE)
            cat(paths[idx.all], "\n")
            get.posterior.diff.ftp(
                readRDS(paths[idx.all]),
                readRDS(paths[idx.ftp])
            )
        },
        by = c("MODEL", "ROUND")
    ]
    saveRDS(object = d_diff, filename_rds2)
}


if (make_plots) {
    # only need to change sligtly, nice!
    .w <- 10; .h <- 12
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
            readRDS(paths[idx.all]),
            readRDS(paths[idx.ftp]),
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
                readRDS(paths[idx.all]),
                readRDS(paths[idx.ftp]),
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
                readRDS(paths[idx.all]),
                readRDS(paths[idx.ftp]),
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

    filename_overleaf <- file.path(out.dir.tables, "overleaf_")
    saveRDS(object = tab, file = filename_overleaf)
}


#####################
catn("End of script")
#####################
