#!/bin/Rscript

self_relative_path <- "src/make_paper_figures.R"

########################
cat("\nStart of: ", self_relative_path, "\n")
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
    library(patchwork)
} |> suppressPackageStartupMessages()


################
#    PATHS     #
################

if (interactive()) {
    gitdir <- here::here()
} else {
    cmd <- commandArgs()
    cmd <- cmd[cmd %like% "file"]
    gitdir <- gsub(paste0("--file=(.*)/", self_relative_path), "\\1", cmd)
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
    "jobname",
    "shared.hyper"
)
args <- args[names(args) %in% opts_vec]
# testing
if (interactive()) {
    # args$out.dir.exact <- "/home/andrea/HPC/ab1820/home/projects/2022/longvl/vl_1000_joint"
    args$out.dir.exact <- "/home/andrea/HPC/ab1820/home/projects/2022/longvl/cmdstan_alpha100sharedhyper_vl_1000"
    args$shared.hyper <- TRUE
}
print(args)
{
    source(file.path(gitdir.functions, "plotting_main_figures.R"))
    source(file.path(gitdir.functions, "plotting_functions.R"))
    source(file.path(gitdir.functions, "postprocessing_helpers.R"))
    source(file.path(gitdir.functions, "postprocessing_tables.R"))
    source(file.path(gitdir.functions, "phsc_vl_helpers.R"))
    source(file.path(gitdir.functions, "paper_statements.R"))
    naturemed_reqs()
}


overwrite <- !interactive()
make_plots <- make_tables <- make_paper_numbers <- TRUE
if (make_paper_numbers) {
    ppr_numbers <- list()
}
fetch.postprocessing.settings.from.args(args)

# Census eligible, participants, and smooth
ncen <- fread(path.census.eligible, select = c("ROUND", "TYPE", "SEX", "AGEYRS", "ELIGIBLE", "ELIGIBLE_SMOOTH")) |>
    subset(ROUND %in% 16:19)
ncen[, ROUND := as.integer(ROUND)]
setnames(ncen, c("TYPE"), c("FC"))

# HIV status and viral load among participants
cols <- c("STUDY_ID", "ROUND", "SEX", "AGEYRS", "FC", "HIV_STATUS", "VL_COPIES", "FIRST_PARTICIPATION")
dall <- get.dall(path.hivstatusvl.r1520) |>
    subset(ROUND %in% 16:19, select = cols) |>
    unique()
dall[, ROUND := as.integer(ROUND)]

# proportions
npar <- summarize.aggregates.dall()
dprop <- merge(npar, ncen)
check_more_elig_than_part <- dprop[, all(N_PART < ELIGIBLE)]
stopifnot(check_more_elig_than_part)

# paths
filename_rds_contrib <- file.path(out.dir.tables, "fullpop_allcontributions.rds")
filename_rds_prevalence <- file.path(out.dir.tables, "fullpop_posteriorquantiles_by_agesexround.rds")
filename_rds_prevalence_agegroup <- file.path(out.dir.tables, "posterior_quantiles_agegroups.rds")
filename_rds_contrib_agegroup <- file.path(out.dir.tables, "fullpop_allcontributions_byagegroup.rds")
filename_rds_supphiv_agegroup <- file.path(out.dir.tables, "posterior_quantiles_suppression_agegroup.rds")

djoint_agegroup <- readRDS(filename_rds_prevalence_agegroup)
dcontrib_agegroup <- readRDS(filename_rds_contrib_agegroup)
dsupp_agegroup <- readRDS(filename_rds_supphiv_agegroup)

######################
# Google drive stuff #
######################

if ( interactive()){
    library(googledrive)
    drive_project_base <-'longivl_paper_figures' 
    drive_outdir <- file.path(drive_project_base, basename(out.dir))
    drive_mkdir(name = drive_outdir)
    drive_find(drive_outdir)
    
    upload_to_googledrive <- function(path){
        drive_upload( media=path, path = file.path(drive_outdir, basename(path)), overwrite = TRUE )
    }
}

###########
# HELPERS #
###########

.load.model.subset <- function(file, expr = NULL) {
    expr <- enexpr(expr)
    readRDS(file) |> subset(eval(expr))
}

########################
catn("=== FIGURE 1 ===")
########################

# only make this figure if ggsn package is installed (but still load files for later)
dcontrib <- .load.model.subset(filename_rds_contrib, MODEL == "run-gp-prevl" & ROUND == 19) |> prettify_labels()
dprev <- .load.model.subset(filename_rds_prevalence, MODEL == "run-gp-prevl" & ROUND == 19) |> prettify_labels()

if (system.file(package = "ggsn") != "") {
    # Map of the Rakai communities
    fig1a.outer <- plot.uganda.map(zoom = "medium", maptype = "toner-lines", labs = TRUE)
    fig1a.inner <- plot.rakai.map(.size = 1, labs = FALSE) +
        theme(
            panel.border = element_rect(colour = "red", size = 1),
            plot.margin = margin(t = 0, b = 0, l = 0, r = 0, unit = "cm"),
            legend.key.size = unit(0.01, "cm"),
            legend.spacing.x = unit(0.01, "cm")
        )
    .delta <- .7
    fig1a <- fig1a.outer + inset_element(
        fig1a.inner,
        left = 1 - .delta, right = 1, bottom = 0, top = .delta, align_to = "panel"
    )

    # pyramid of census eligible, participants, and smooth
    fig1b <- plot.pyramid.bysexround(dprop[ROUND %in% c(19)],
        .ylab = "Number of census eligible individuals",
        NUM = NULL,
        DEN = "ELIGIBLE",
        percent_lab = FALSE
    ) +
        facet_grid(ROUND_LAB ~ FC_LAB, labeller = labeller(ROUND_LAB = round_labs2), scales = "free_x") +
        geom_hline(yintercept = 0, color = "black") +
        geom_line(aes(y = ELIGIBLE_SMOOTH * (1 - 2 * (SEX == "M")), color = SEX_LAB)) +
        my_labs(color = "Gender") +
        theme(legend.position = "none")

    # fig 1c
    fig1c <- plot_prevalenceandcontrid(dprev, dcontrib)

    fig1 <- {{
        (fig1a | fig1b) + plot_layout(widths = c(1, 1))
    } / {
        fig1c
    }} + plot_layout(heights = c(1, 3), widths = c(1, 1, 1)) &
        theme(
            legend.key.size = unit(3, "mm"),
            legend.margin = margin(1.6, 1.6, 1.6, 1.6, unit = "mm"),
        ) +
            nm_reqs + t_nomargin
    filename <- paste0("main_figure_populationcomposition.pdf")
    ggsave_nature(p = fig1, filename = filename, LALA = out.dir.figures, w = 19.5, h = 23)
    upload_to_googledrive(path=file.path(out.dir.figures, pdf2png(filename)) )
}


# Statment in section 1.
{
    tmp1 <- paper_statements_female_prevalence(djoint_agegroup)
    tmp2 <- paper_statements_female_contributions_prevalence(dcontrib_agegroup)
    sprintf(
        "while on average %s of all women had HIV in inland communities, we found that %s of all individuals with HIV in inland communities were female, and it is the latter proportion that is relevant for intervention planning in these communities. Vice versa, while on average %s of all women had HIV in fishing communities, we found that %s of all individuals with HIV in fishing communities were female.",
        tmp1[LOC == "inland" & SEX == "F"]$CELL, tmp2[LOC == "inland" & SEX == "F"]$CELL,
        tmp1[LOC == "fishing" & SEX == "F"]$CELL, tmp2[LOC == "fishing" & SEX == "F"]$CELL
    )
}

########################
catn("=== FIGURE 2 ===")
########################

djoint <- {
    file.path(out.dir.tables, "fullpop_posteriorquantiles_by_agesexround.rds") |>
        readRDS()
}
fig2a <- plot.fit.weighted.by.ftpstatus(djoint, "run-gp-supp-hiv")
tmp <- paper_statements_suppression_above_959595(djoint)

# fig2b <- {
#     file.path(out.dir.figures, "posterior_suppressionincrease_vsround16.rds") |>
#     readRDS()
# } |> plot.relative.suppression.vs.round16.ratio()
# filename <- paste0('main_figure_changesinsuppression.pdf')
# ggsave_nature(p=fig1, filename=filename, LALA=out.dir.figures, w=21, h=16)

# fig2b

fig2 <- plot.main.suppression.among.plhiv(DT = djoint, type = "point", unaids = TRUE, rev = TRUE, m = -9)
filename <- paste0("main_figure_suppression_plhiv_r1619.pdf")
cmd <- ggsave_nature(p = fig2, filename = filename, LALA = out.dir.figures, w = 17, h = 16)
if(interactive()) upload_to_googledrive(path=file.path(out.dir.figures, pdf2png(filename)) )

########################
catn("=== FIGURE 3 ===")
########################

dcens <- copy(ncen) |> setnames("FC", "LOC")
fig3 <- plot_propofpop_of_viraemic_byagesex_stratbycommround(DT=djoint, colorby='ROUND_LAB', cri=TRUE)
dlabels <- labels_propofpop_of_viraemic_byagesex_stratbycommround()
dtriangles <- point_propofpop_of_viraemic_byagesex_stratbycommround()
.f <- function(add_labels=FALSE, add_ages=FALSE){
    list(
        if(add_labels)  geom_text(data=dlabels, aes(label=CELL), color="black", x=Inf, y=Inf, hjust=1, vjust=1),
        if(add_ages)  geom_point(data=dtriangles, aes(x=M2, y=0), size=3, pch=17, alpha=.8), 
        NULL
    ) -> out 
    out[!sapply(out, is.null)]
}
fig3 <- fig3 + .f(add_ages=FALSE, add_labels=FALSE)
filename <- paste0("main_figure_profile_nonsuppressed.pdf")
cmd <- ggsave_nature(p = fig3, filename = filename, LALA = out.dir.figures, w = 17, h = 16)
if(interactive()) upload_to_googledrive(path=file.path(out.dir.figures, pdf2png(filename)) )

########################
catn(" FIGURE for KATE")
########################

{
    horizontal <- TRUE
    fig_pyr <- plot.pyramid.bysexround(dprop[ROUND == 19 & FC == "inland"],
        .ylab = "Number of census eligible individuals",
        NUM = "ELIGIBLE", DEN = "ELIGIBLE",
        percent_lab = FALSE
    ) +

        facet_grid(. ~ FC_LAB, labeller = labeller(ROUND_LAB = round_labs2, FC_LAB = community_dictionary$long), scales = "free_x") +
        geom_hline(yintercept = 0, color = "black") +
        scale_y_continuous(labels = abs, limits = c(-800, 800)) +
        theme(legend.key.size = unit(.3, "cm"), strip.text.x = element_blank()) +
        slides_reqs

    fig_kate2 <- plot_prevalenceandcontrid(
        sec_name = "Contribution of each age group to HIV prevalence",
        dprev[LOC == "inland"], dcontrib[LOC == "inland"],
        legend.key.size = unit(0.3, "cm"),
        sec_axis_scale = .1,
        slides = TRUE,
        extra_fig = fig_pyr,
        CrI = FALSE
    )
    plot_prevalence <- fig_kate2
    plot_prevalence

    # (fig1b + theme(strip.text.y=element_blank()) | fig_kate2) + plot_layout(ncol=2, widths=c(1,2.5))
    # if(!horizontal){
    #     plot_prevalence <- (
    #         (fig_pyr + theme(strip.text.y=element_blank()) + labs(tag='a')) /
    #         fig_kate2 + labs(tag='b')
    #     ) +  plot_layout(ncol=1, heights=c(1,2.4)) + slides_reqs
    # } else {
    #     plot_prevalence <- (
    #         (fig_pyr + theme(strip.text.y=element_blank()) + labs(tag='a')) |
    #         fig_kate2 + labs(tag='b')
    #     ) +  plot_layout(nrow=1, ncol=2, widths=c(1,2.4)) + slides_reqs
    # }
    rm(horizontal)
}
ggsave_nature(p = plot_prevalence, filename = "whopepfar_prevalence.pdf", LALA = out.dir.figures, w = 25, h = 15)


{
    dprev_vir <- .load.model.subset(filename_rds_prevalence, MODEL == "run-gp-supp-pop" & ROUND == 19) |> prettify_labels()
    dcontrib_vir <- .load.model.subset(filename_rds_contrib, MODEL == "run-gp-supp-pop" & ROUND == 19) |> prettify_labels()
    dprev_vir2 <- .load.model.subset(filename_rds_prevalence, MODEL == "run-gp-supp-hiv" & ROUND == 19) |> prettify_labels()
    dprev_vir2[, (c("M", "CU", "CL")) := lapply(.SD, function(x) 1 - x), .SDcols = c("M", "CL", "CU")]

    fig_kate_1c <- plot_suppandcontrib(
        dprev_vir[LOC == "inland"],
        dcontrib_vir[LOC == "inland"],
        sec_name = "Contribution of each age group to people with unsuppressed viral load",
        prevalence.label = "Prevalence of viraemia among population",
        remove.legend = TRUE,
        CrI = FALSE,
        slides = TRUE
    )
    fig_kate_1d <- plot_suppandcontrib(
        dprev_vir2[LOC == "inland"],
        dcontrib_vir[LOC == "inland"],
        sec_name = "Contribution of each age group to people with unsuppressed viral load",
        prevalence.label = "Prevalence of viraemia among HIV positive population",
        slides = TRUE,
        CrI = FALSE,
        sec_axis_scale = .069,
        UNAIDS = TRUE
    )

    # plot_suppression <- (
    #     (fig_kate_1c +  labs(tag="Considering entire population") + nm_reqs + labs(x="LALALA")) /
    #     (fig_kate_1d +  labs(tag="Considering PLHIV only (as per 95-95-95 goals)") + slides_reqs )
    # ) + plot_layout(guides="collect", ncol=1, heights=c(1,1.1)) + theme(legend.position="bottom")
    plot_suppression <- fig_kate_1d
}
ggsave_nature(p = plot_suppression, filename = "whopepfar_suppression.pdf", LALA = out.dir.figures, w = 23, h = 17)


################################################
catn("Make table with study pop characteristics")
#################################################

if (0) {
    filename_rds <- file.path(out.dir.tables, "fullpop_allcontributions.rds")
    fig1c <- readRDS(filename_rds) |>
        subset(ROUND == 19 & LOC == "inland") |>
        plot.agesex.contributions.by.roundcomm(label = "run-gp-prevl", include_baseline = TRUE)
}
cat("not completed...")


################################
catn("UNAIDS-95^3 alternatives")
################################

p <- plot.proposed.metric()
filename <- paste0("main_figure_proposedmetric.pdf")
cmd <- ggsave_nature(p = p, filename = filename, LALA = out.dir.figures, w = 17, h = 22)
if (interactive()) system(zathura2gthumb(cmd))

#####################
catn("Main table 1")
#####################

dtable1 <- make.main.table.contributions(add_asterisks_unaids = FALSE)
if(interactive()){
    write.to.googlesheets(dtable1, "Table1")
}

filename_tex <- file.path(out.dir.tables, "main_table_contributions_hiv.tex")
write.to.tex(dtable1, file = filename_tex)
p <- table.to.plot(dtable1)
cmd <- ggsave2(p = p, file = tex2pdf(filename_tex), LALA = out.dir.tables, w = 20, h = 10)
if (interactive()) system(zathura2gthumb(cmd))


dcontrib_agegroup[ROUND %in% c(16, 19) & MODEL == "run-gp-supp-pop" & AGEGROUP %in% c("15-19", "20-24"),
    sum(M),
    by = c("ROUND", "SEX", "LOC")
] |>
    dcast(SEX + LOC ~ ROUND)


####################
catn("Main table 2")
####################

filename_table <- file.path(out.dir.tables, "table_reductionHIVandUNSUPP.rds")
tab_merge <- readRDS(filename_table)
write.to.googlesheets(tab_merge, sheet = "Table2")
