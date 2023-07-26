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
    library(patchwork)
}


################
#    PATHS     #
################

gitdir <- here::here()
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
    "indir",
    "round",
    "jobname",
    "only.firstparticipants"
)
args <- args[names(args) %in% opts_vec]
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
make_plots <- TRUE
make_tables <- TRUE

make_paper_numbers <- TRUE
if (make_paper_numbers) {
    ppr_numbers <- list()
}

with(args, {
    VL_DETECTABLE <<- detectable.viral.load
    VIREMIC_VIRAL_LOAD <<- viremic.viral.load
    # need to check that we have both samples for "" and "_firstpart"
    out.dir <<- gsub(out.dir.exact, "_firstpart$", "_joint")
    if(is.na(out.dir.exact)){
        out.dir <<- file.path(
            out.dir.prefix, 
            gsub('_firstpart$','_joint',jobname))
    }
})
stopifnot("out.dir must end in _joint"= out.dir %like% '_joint$')
indir.ftp <- gsub( '_joint', "_firstpart", out.dir)
indir.all <- gsub( '_joint', "", out.dir)
stopifnot("did not find 2 directories necessary to run joint analysis"=all(dir.exists(c(indir.ftp, indir.all))))
out.dir.figures <- file.path(out.dir, "figures")
out.dir.tables <- file.path(out.dir, "tables")


# Census eligible, participants, and smooth
ncen <- fread(path.census.eligible, select = c( "ROUND", "TYPE", "SEX", "AGEYRS", "ELIGIBLE", "ELIGIBLE_SMOOTH")) |>
    subset(ROUND %in% 16:19)
ncen[,ROUND := as.integer(ROUND)]
setnames(ncen, c('TYPE'), c('FC'))

# HIV status and viral load among participants
cols <- c("STUDY_ID", "ROUND", "SEX", "AGEYRS", "FC", 'HIV_STATUS', 'VL_COPIES', 'FIRST_PARTICIPATION')
dall <- get.dall(path.hivstatusvl.r1520) |> 
    subset(ROUND %in% 16:19, select=cols) |> unique()
dall[,ROUND := as.integer(ROUND)]

# proportions
npar <- summarize.aggregates.dall()
dprop <- merge(npar, ncen)
check_more_elig_than_part <- dprop[, all(N_PART < ELIGIBLE) ] 
stopifnot(check_more_elig_than_part)

# paths
filename_rds_contrib <- file.path(out.dir.tables, "fullpop_allcontributions.rds")
filename_rds_prevalence <- file.path(out.dir.tables, "fullpop_posteriorquantiles_by_agesexround.rds")
filename_rds_prevalence_agegroup <- file.path(out.dir.tables, "posterior_quantiles_agegroups.rds")
filename_rds_contrib_agegroup <- file.path(out.dir.tables, "fullpop_allcontributions_byagegroup.rds")

djoint_agegroup <- readRDS(filename_rds_prevalence_agegroup)
dcontrib_agegroup <- readRDS(filename_rds_contrib_agegroup)

###########
# HELPERS #
###########

.load.model.subset <- function(file, expr=NULL){
    expr <- enexpr(expr)
    readRDS(file) |> subset( eval(expr) )
}

########################
catn("=== FIGURE 1 ===")
########################

# Map of the Rakai communities
fig1a <- plot.rakai.map(.size=2) + 
    theme(panel.border = element_rect(colour = "red",size=2))

subfig1a <- plot.uganda.map(zoom="medium")


# pyramid of census eligible, participants, and smooth
fig1b <- plot.pyramid.bysexround( dprop[ROUND %in% c(19)], 
    .ylab = 'Number of participants among census eligible individuals',
    NUM="N_PART",
    DEN='ELIGIBLE',
    percent_lab = FALSE) +
    facet_grid(ROUND_LAB ~ FC_LAB, labeller = labeller(ROUND_LAB=round_labs2), scales='free_x') +
    geom_hline(yintercept=0, color='black') +
    geom_line( aes(y=ELIGIBLE_SMOOTH * (1 - 2*(SEX == 'M')), color=SEX_LAB )) +
    my_labs(color="Gender")

# fig 1c
dcontrib <- .load.model.subset(filename_rds_contrib, MODEL=="run-gp-prevl" & ROUND == 19) |> prettify_labels()
dprev <- .load.model.subset(filename_rds_prevalence, MODEL=="run-gp-prevl" & ROUND == 19) |> prettify_labels()
# p_i <- plot_2yaxis_hist_lines(dcontrib[LOC == 'inland'],  dprev[LOC == 'inland'], sec_name="")
# p_f <- plot_2yaxis_hist_lines(dcontrib[LOC == 'fishing'], dprev[LOC == 'fishing']) + labs(y="")
fig1c <- plot_prevalenceandcontrid(dprev, dcontrib)  + 
    labs(tag = "c")

# patchwork:
fig1 <- ( 
    (subfig1a + labs(tag="a") | fig1a | fig1b + theme(legend.position = "none") + labs(tag='b') ) + plot_layout(widths = c(1, 1, 1.5))
)/(
    # (p_i + labs(tag='c') + theme(legend.position="none") | p_f) & plot_layout(guides="collect")
    fig1c
)  
fig1 <- fig1 + plot_layout(heights = c(1,2)) & theme(legend.key.size = unit(0.4, "cm")) + nm_reqs  
filename <- paste0('main_figure_populationcomposition.pdf')
ggsave_nature(p=fig1, filename=filename, LALA=out.dir.figures, w=21, h=19)

# Statment in section 1.
{
    tmp1 <- paper_statements_female_prevalence(djoint_agegroup)
    tmp2 <- paper_statements_female_contributions_prevalence(dcontrib_agegroup)
    sprintf( 
        "while on average %s of all women had HIV in inland communities, we found that %s of all individuals with HIV in inland communities were female, and it is the latter proportion that is relevant for intervention planning in these communities. Vice versa, while on average %s of all women had HIV in fishing communities, we found that %s of all individuals with HIV in fishing communities were female.", 
        tmp1[LOC == "inland" & SEX=="F"]$CELL, tmp2[LOC == "inland" & SEX=="F"]$CELL, 
        tmp1[LOC == "fishing" & SEX=="F"]$CELL, tmp2[LOC == "fishing" & SEX=="F"]$CELL
    )
}

########################
catn("=== FIGURE 2 ===")
########################

# don't really like it atm 

filename_rds <- file.path(out.dir.tables, "fullpop_posteriorquantiles_by_agesexround.rds")
djoint <- readRDS(filename_rds)
fig2a <- plot.fit.weighted.by.ftpstatus(djoint, "run-gp-supp-hiv")
tmp <- paper_statements_suppression_above_959595(djoint)

filename_rds  <- file.path(out.dir.figures, "posterior_suppressionincrease_vsround16.rds")
tmp <- readRDS(filename_rds)
fig2b <- plot.relative.suppression.vs.round16.ratio(tmp) 

filename <- paste0('main_figure_changesinsuppression.pdf')
ggsave_nature(p=fig1, filename=filename, LALA=out.dir.figures, w=21, h=16)

# fig2b

########################
catn(" FIGURE for KATE")
########################

{
    horizontal <- TRUE
    fig_pyr <- plot.pyramid.bysexround( dprop[ROUND == 19 & FC=='inland'], 
        .ylab = 'Number of census eligible individuals',
        NUM="ELIGIBLE", DEN='ELIGIBLE',
        percent_lab = FALSE) +

        facet_grid(. ~ FC_LAB, labeller = labeller(ROUND_LAB=round_labs2, FC_LAB = community_dictionary$long), scales='free_x') +
        geom_hline(yintercept=0, color='black') +
        scale_y_continuous(labels=abs, limits=c(-800, 800)) + 
        theme(legend.key.size = unit(.3, "cm"), strip.text.x=element_blank()) +
        slides_reqs

    fig_kate2 <- plot_prevalenceandcontrid(
        sec_name="Contribution of each age group to HIV prevalence",
        dprev[LOC=="inland"], dcontrib[LOC=="inland"],
        legend.key.size=unit(0.3, "cm"),
        sec_axis_scale = .1,
        slides=TRUE,
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
ggsave_nature(p = plot_prevalence, filename="whopepfar_prevalence.pdf", LALA="~/Downloads", w=25, h=15)


{
    dprev_vir <- .load.model.subset(filename_rds_prevalence, MODEL=="run-gp-supp-pop" & ROUND == 19) |> prettify_labels()
    dcontrib_vir <- .load.model.subset(filename_rds_contrib, MODEL=="run-gp-supp-pop" & ROUND == 19) |> prettify_labels()
    dprev_vir2 <- .load.model.subset(filename_rds_prevalence, MODEL=="run-gp-supp-hiv" & ROUND == 19) |> prettify_labels()
    dprev_vir2[, (c("M", "CU", "CL")) := lapply(.SD, function(x) 1-x), .SDcols = c("M", "CL", "CU")]

    fig_kate_1c <- plot_suppandcontrib(
        dprev_vir[LOC=="inland"], 
        dcontrib_vir[LOC=="inland"],
        sec_name = "Contribution of each age group to people with unsuppressed viral load",
        prevalence.label = "Prevalence of viraemia among population",
        remove.legend = TRUE,
        CrI = FALSE,
        slides= TRUE
    )
    fig_kate_1d <- plot_suppandcontrib(
        dprev_vir2[LOC=="inland"], 
        dcontrib_vir[LOC=="inland"],
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
ggsave_nature(p = plot_suppression, filename="whopepfar_suppression.pdf", LALA="~/Downloads", w=23, h=17)


################################################
catn("Make table with study pop characteristics")
#################################################

filename_rds <- file.path(out.dir.tables, "fullpop_allcontributions.rds")
fig1c <- readRDS(filename_rds) |> 
    subset(ROUND==19 & LOC=='inland') |>
    plot.agesex.contributions.by.roundcomm(label = "run-gp-prevl", include_baseline =TRUE) 


