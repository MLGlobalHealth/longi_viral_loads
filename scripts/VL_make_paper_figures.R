#!/bin/Rscript

self_relative_path <- "scripts/VL_make_paper_figures.R"

########################
cat("\nStart of: ", self_relative_path, "\n")
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
    library(patchwork)
    library(osmdata)
})


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

# Census-eligible, participants, and smooth
cols <- c("ROUND", "TYPE", "SEX", "AGEYRS", "ELIGIBLE", "ELIGIBLE_SMOOTH")
ncen <- fread(path.census.eligible, select=cols ) |>
    subset(ROUND %in% 16:19) |>
    setnames("TYPE", "FC")
ncen[, ROUND := as.integer(ROUND)]

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
.fp <- function(x) file.path(out.dir.tables, x)
filename_rds_contrib                 <- .fp("fullpop_allcontributions.rds")
filename_rds_prevalence              <- .fp("fullpop_posteriorquantiles_by_agesexround.rds")
filename_rds_prevalence_agegroup     <- .fp("posterior_quantiles_agegroups.rds")
filename_rds_contrib_agegroup        <- .fp("fullpop_allcontributions_byagegroup.rds")
filename_rds_supphiv_agegroup        <- .fp("posterior_quantiles_suppression_agegroup.rds")
filename_rds_supphiv_agegroup_custom <- .fp("posterior_quantiles_suppression_agegroup_custom.rds")

djoint_agegroup         <- readRDS(filename_rds_prevalence_agegroup)
dcontrib_agegroup       <- readRDS(filename_rds_contrib_agegroup)
dsupp_agegroup          <- readRDS(filename_rds_supphiv_agegroup)
dsupp_agegroup_custom   <- readRDS(filename_rds_supphiv_agegroup_custom)

######################
# Google drive stuff #
######################

if ( interactive()){
    library(googledrive)
    drive_project_base <-'longivl_paper_figures' 
    if( ! basename(out.dir) %in% drive_ls(drive_project_base)$name ){
        drive_mkdir(name = basename(out.dir),path=drive_project_base, overwrite = FALSE)
    }else{
        cat("GoogleDrive job directory already present\n")
    }
    upload_to_googledrive <- function(path){
        drive_upload( media=path, path = file.path(drive_project_base, basename(path)), overwrite = TRUE )
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
dcontrib    <- .load.model.subset(filename_rds_contrib, MODEL == "run-gp-prevl" & ROUND == 19) |> prettify_labels()
dprev       <- .load.model.subset(filename_rds_prevalence, MODEL == "run-gp-prevl" & ROUND == 19) |> prettify_labels()

if (system.file(package = "ggsn") != "" && usr == "andrea") {

    stadia_api_key <- Sys.getenv("STADIA_API_KEY")
    ggmap::register_stadiamaps(stadia_api_key)
    # package_version("ggplot2")
    library(ggmap)
    citation("ggmap")
    #  © Stadia Maps © Stamen Design © OpenMapTiles © OpenStreetMap contributor

    fig1a_inset <- plot_all_maps(sizes=list(points=1, text=3))
    # fig1a_columns <- plot_all_maps(type="columns")
    # fig1a_columns_list <- plot_all_maps(type="columns", return_list=TRUE)

    fig1b_h <- plot_paper_population_pyramid(layout="horizontal")

    # fig_map_and_hist <- patchwork::wrap_plots(
    #     plot_spacer(), 
    #     fig1a_columns,
    #     plot_spacer(), 
    #     fig1b_v
    # ) + plot_layout(widths=c(-3,7, -3, 1))

    if(0) {
        # alternative maps
        .byrow <- TRUE
        fig1b_v_list <- plot_paper_population_pyramid(layout="vertical", return_list=TRUE)
        fig_map_and_hist <- patchwork::wrap_plots(
            fig1a_columns_list[[1]] + nm_reqs,
            fig1a_columns_list[[2]] + nm_reqs + theme(legend.key.size = unit(0.05, "cm"), legend.spacing.x = unit(0.05, "cm")),
            fig1b_v_list[[1]] + nm_reqs + labs(x="Age") + if(.byrow==FALSE){labs(y="")}else{NULL},
            fig1b_v_list[[2]] + nm_reqs,
            nrow = 2, 
            byrow = .byrow,
            widths=c(1,1)
        ) + plot_annotation(tag_levels = "a") 
        filename <- "main_figure_mapandhist.pdf"
        cmd <- ggsave_nature(p = fig_map_and_hist, filename = filename, LALA = out.dir.figures, w = 18.5, h = 16)
    }
    # system(zathura2gthumb(cmd))

    # fig 1c
    fig1c <- plot.figure.main.prevalence(subtitles=TRUE, size=9)
    filename <- paste0("main_figure_hivprevalenceonly.pdf")
    cmd <- ggsave_nature(p = fig1c, filename = filename, LALA = out.dir.figures, w = 19.5, h = 16)
    # system(zathura2gthumb(cmd))
    upload_to_googledrive(path=file.path(out.dir.figures, pdf2png(filename)) )
    }
    # plot.figure.main.prevalence(subtitles=FALSE)

    # fig1 <- {
    #     fig1_top / fig1c
    # } + plot_layout(heights = c(1, 2.7), widths = c(1, 1, 1)) &
    #     theme(
    #         legend.key.size = unit(3, "mm"),
    #         legend.margin = margin(1.6, 1.6, 1.6, 1.6, unit = "mm"),
    #     ) + t_nomargin

    ls_top <- list( nm_reqs, t_nomargin, theme(legend.key.size=unit(3, 'mm'), legend.position='none'))
    ls <- list( scale_x_continuous(), scale_y_continuous(), nm_reqs, t_nomargin)
    fig1_top <- (fig1a_inset | fig1b_h) + plot_layout(widths = c(1, 1))
    fig1 <- ggarrange(
        fig1_top + ls_top,
        fig1c + ls,
        nrow=2, heights=c(1,2), align="v" )

    filename <- paste0("main_figure_populationcomposition.pdf")
    cmd <- ggsave_nature(p = fig1, filename = filename, LALA = out.dir.figures, w = 19.5, h = 22)
    system(gsub("zathura","inkscape",cmd))
    # system(zathura2gthumb(cmd))
    upload_to_googledrive(path=file.path(out.dir.figures, pdf2png(filename)) )

# }


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
.null <- paper_statements_suppression_above_959595(djoint)
.yl <- function() scale_y_continuous(limits=c(0, 1), expand=expansion(c(0,0)), labels=scales::label_percent() )

fig2a <- list(
    plot.main.suppression.among.plhiv(DT=djoint[LOC=='inland'] ,m=-5,joint=TRUE) + .yl() + labs(tag = "a"),
    plot.main.suppression.among.plhiv(DT=djoint[LOC=='fishing'],m=-5,joint=TRUE)+ .yl() + labs(tag = "b")
) |> ggarrange( plotlist=_, ncol=1, legend="bottom", common.legend = TRUE  )
fig2b <- list(
    hist_prevalence_by_age_group_custom(dsupp_agegroup_custom[ LOC == 'inland']) +
        theme(strip.text.x=element_text(color="white")) + .yl() + labs(y=NULL) + nm_reqs,
    hist_prevalence_by_age_group_custom(dsupp_agegroup_custom[ LOC == 'fishing']) + 
        theme(strip.text.x=element_text(color="white")) + .yl() + labs(y=NULL) + nm_reqs
) |> ggarrange(plotlist=_,ncol=1, legend="bottom", common.legend = TRUE)

fig2new <- ggarrange(
    fig2a, 
    fig2b + theme(legend.position = "none"), 
    ncol=2, 
    widths = c(1.4, 1))
filename <- paste0("main_figure_suppression_plhiv_r1619.pdf")
cmd <- ggsave_nature(p = fig2new,
    filename = filename,
    LALA = out.dir.figures,
    w = 19, h = 16)
# system(cmd)
if(interactive()) 
    upload_to_googledrive(path=file.path(out.dir.figures, pdf2png(filename)) )

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
        .ylab = "Number of census-eligible individuals",
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

    fig_kate_1d <- plot_suppandcontrib(
        dprev_vir2[LOC == "inland"],
        dcontrib_vir[LOC == "inland"],
        sec_name = "Age and gender profile of people exhibiting vireamia\n(black line, right y-axis)",
        prevalence.label = "Prevalence of viraemia among HIV positive population in each age band\n(orange bars, left y-axis)",
        reqs = nm_reqs,
        remove.legend = TRUE,
        CrI = TRUE,
        sec_axis_scale = .08,
        UNAIDS = FALSE
    ) |> annotate_figure(top = text_grob( community_dictionary$longest2['I'], size = 9))
    fig_kate_1f <- plot_suppandcontrib(
        dprev_vir2[LOC == "fishing"],
        dcontrib_vir[LOC == "fishing"],
        sec_name = "Age and gender profile of people exhibiting vireamia\n(black line, right y-axis)",
        prevalence.label = "Prevalence of viraemia among HIV positive population in each age band\n(orange bars, left y-axis)",
        reqs = nm_reqs,
        CrI = TRUE,
        sec_axis_scale = .08,
        UNAIDS = FALSE
    ) |> annotate_figure(top=text_grob( community_dictionary$longest2['F'] , size = 9))

    # plot_suppression <- fig_kate_1d + theme(legend.position="none")
    plot_suppression2 <- fig_kate_1d / fig_kate_1f
}

if(0){
    # used to 000e84161c933cdd1aba319a3571e36cfa3d63dc
    ggsave_nature(p = plot_suppression, filename = "whopepfar_suppression.pdf", LALA = out.dir.figures, w = 23, h = 17)
    cmd <- ggsave_nature(p = plot_suppression2, filename = "whopepfar_suppression_comms.pdf", LALA = out.dir.figures, w = 30 , h = 23)
    # system(cmd)
}

dcontrib_vir1619 <- .load.model.subset(
    filename_rds_contrib,
    MODEL == "run-gp-supp-pop" & ROUND %in%  c(16, 19)
) |> prettify_labels()
p_sex <- make_figure_3b(facet_var = "SEX_LAB")
# p_loc <- make_figure_3b(facet_var = "LOC_LAB")
# fig3_redo <- ((fig_kate_1d / fig_kate_1f)/ p_loc ) + plot_layout(ncol=1, heights=c(.9,1,.4))
fig3_redo <- ((fig_kate_1d / fig_kate_1f)/ p_sex ) + plot_layout(ncol=1, heights=c(.9,1,.4))

cmd <- ggsave_nature(
    p = fig3_redo,
    filename = "main_figure_suppression_redo.pdf",
    LALA = out.dir.figures, w = 20 , h = 22
)
system(cmd)

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

dtable1 <- make_main_table_contributions(add_asterisks_unaids = FALSE)
write.to.googlesheets(dtable1, "Table1")

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


########################
catn("Supp table 2 on mean ages")
########################

dmeanage <- file.path(out.dir.tables, "mean_ages_plhiv_viraemic.rds") |>
    readRDS()

tab <- make.supp.table.meanage()
filename_table <- file.path(out.dir.tables, "table_agemeanstd_prev_unsupp.rds")
write.to.googlesheets(tab, sheet = "SuppTable2")
