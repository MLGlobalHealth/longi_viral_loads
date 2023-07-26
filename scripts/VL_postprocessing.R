#!/bin/Rscript

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

# helpers
source( file.path(gitdir.functions,'plotting_main_figures.R') )
source( file.path(gitdir.functions,'postprocessing_helpers.R') )
source( file.path(gitdir.functions,'phsc_vl_helpers.R') )
naturemed_reqs() # stores them in nm_reqs

# command line options, stored in args. Then subset
opts_vec <- c(
    'viremic.viral.load',
    'detectable.viral.load',
    'out.dir.prefix',
    'out.dir.exact',
    'indir',
    'round',
    'jobname',
    'only.firstparticipants')
args <- args[ names(args) %in% opts_vec]
print(args)
if(interactive()){
    args$jobname <- 'cmdstan_vl_1000_firstpart'
}

# Set global variables
with(args, {
    VL_DETECTABLE <<- detectable.viral.load
    VIREMIC_VIRAL_LOAD <<- viremic.viral.load
    vl.out.dir <<- out.dir.exact
    if(is.na(out.dir.exact)){
        vl.out.dir <<- file.path(out.dir.prefix, jobname)
    }
})
stopifnot('Specified vl.out.dir does not exist' = file.exists(vl.out.dir))

# other setup
#------------

make_paper_numbers <- TRUE

if(make_paper_numbers)
    ppr_numbers <- list()

# output directory with rda files
vl.out.dir.figures <- file.path(vl.out.dir, 'figures')
vl.out.dir.tables <- file.path(vl.out.dir, 'tables')
dir.create(vl.out.dir.tables) |> suppressWarnings()
dir.create(vl.out.dir.figures) |> suppressWarnings()

################
#     MAIN     #
################

rda_files <- list.files.from.output.directory('.rda', dir=vl.out.dir, args=args)

# get data
dall <- get.dall(path.hivstatusvl.r1520, make_flowchart=TRUE)


# Get census eligible and compare with participants 
# __________________________________________________

# we need to chose whether to use the smoothed version, or the actual counts (`.count`)
by_vars <- c('ROUND', 'SEX_LABEL', 'LOC_LABEL')
dcens <- get.census.eligible()
dcens[, .(N_ELIGIBLE=sum(ELIGIBLE)), by=by_vars] |> 
    dcast(ROUND + SEX_LABEL ~ LOC_LABEL, value.var='N_ELIGIBLE' ) |>
    kable(caption = 'number of census eligible by round, sex and location')

# merge with number of participant -> positives -> not suppressed
dcens <- get.participants.positives.unsuppressed(dall) |> 
    merge(dcens, y=_, by=c(by_vars, 'AGE_LABEL'), all.x=TRUE, all.y=TRUE)
last.round <- dcens[, max(ROUND)] 

# deprecated
if(0){
    cat('--- Make UNAIDS objectives plot ---\n')
    tmp <- make.unaids.plots(DT=dcens)

    cat('--- Empirical CDF for HIV positive participants with suppressed viral load---\n')
    tmp <- plot.empirical.prob.of.suppression.with.age(dcens)
}

# Summarised analyes
# __________________

cat('--- Plot Posteriors for fishing analyses ---\n')

p.list <- plot.all.gps(loc='fishing')
filename <- 'main_allanalyses_fishing.pdf'
ggsave2(p.list[['noraw']], file=filename, LALA=vl.out.dir.figures , w=20, h=24, u='cm')
filename <- 'main_allanalyses_fishing_withraw.pdf'
ggsave2(p.list[['withraw']], file=filename, LALA=vl.out.dir.figures , w=20, h=24, u='cm')

cat('--- Plot Posteriors for inland analyses ---\n')

p.list <- plot.all.gps(loc='inland')
filename <- 'main_allanalyses_inland.pdf'
ggsave2(p.list[['noraw']], file=filename, LALA=vl.out.dir.figures, w=20, h=24, u='cm')
filename <- 'main_allanalyses_inland_withraw.pdf'
ggsave2(p.list[['withraw']], file=filename, LALA=vl.out.dir.figures, w=20, h=24, u='cm')



# deprecated

if(0){ 
    # Unaids table
    # ____________
    cat('--- Make UNAIDS goals table ---\n')
    tmp <- make.table.unaids.goals()
    print(xtable::xtable(tmp), 
          include.rownames=FALSE, 
          hline.after=c(-1, seq(0, nrow(tmp), nrow(tmp)/4)),
          file=file.path(vl.out.dir, 'tables', 'main_unaids_table.tex'))
    xtable::xtable(tmp)

    # Compare Prevalence of suppression among HIV+
    # ____________________________________________
    cat('-- Comparing prevalence of suppression among HIV+ ---\n')
    tmp <- vl.vlprops.by.comm.gender.loc(dall)
    vlc <- tmp$DT
    tmp <- make.map.220810(dall, vlc, 'PVLNSofHIV_MEAN')
}

