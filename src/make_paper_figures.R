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
opts_vec <- c("viremic.viral.load", "detectable.viral.load", "out.dir.prefix", "indir", "round", "jobname", "only.firstparticipants")
args <- args[names(args) %in% opts_vec]

source(file.path(gitdir.functions, "plotting_main_figures.R"))
source(file.path(gitdir.functions, "plotting_functions.R"))
source(file.path(gitdir.functions, "postprocessing_helpers.R"))
source(file.path(gitdir.functions, "postprocessing_tables.R"))
source(file.path(gitdir.functions, "phsc_vl_helpers.R"))
source(file.path(gitdir.functions, "paper_statements.R"))
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

## output directories (maybe pass as args?)
out.dir <- args$out.dir.prefix
out.dir <- file.path(out.dir, paste0("vl_", VIREMIC_VIRAL_LOAD, "_joint"))
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
subfig1a | fig1a 

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

# HIV prevalence over time (maybe round 16, 19 ?)
if(0){
    contribution <- TRUE
    fig1c <- if(contribution){
        filename_rds <- file.path(out.dir.tables, "fullpop_allcontributions.rds")
        fig1c <- readRDS(filename_rds) |> 
            subset(ROUND %in% c(16, 19)) |>
            plot.agesex.contributions.by.roundcomm(label = "run-gp-prevl", include_baseline =TRUE) 
    }else{
        filename_rds <- file.path(out.dir.tables, "fullpop_posteriorquantiles_by_agesexround.rds")
        fig1c <- readRDS(filename_rds) |>
            subset(ROUND %in% c(16, 19)) |>
            plot.fit.weighted.by.ftpstatus("run-gp-prevl",include_baseline = TRUE) 
    } + facet_grid(LOC_LAB ~ ROUND_LAB, labeller = labeller(ROUND_LAB=round_labs2)) 

    fig1 <- (((fig1a + plot_layout(guides='keep') | fig1b) + plot_layout(widths = c(1, 2)))/ fig1c) + 
        plot_layout(heights=c(1,1.2), guides='collect') + 
        plot_annotation(tag_levels = 'a') 
    fig1 <- fig1 & nm_reqs + theme(plot.tag = element_text(size=8, face='bold'), legend.key.size = unit(0.4, "cm"))
    fig1
}


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



################################################
catn("Make table with study pop characteristics")
#################################################
