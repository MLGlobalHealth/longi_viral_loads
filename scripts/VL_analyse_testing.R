#!/bin/Rscript

################
# DEPENDENCIES #
################

library(data.table) 
library(ggplot2)
library(Hmisc)
library(rstan)
library(here)
library(knitr)
library(patchwork)


################
#   OPTIONS    #
################

option_list <- list(
    optparse::make_option(
        "--refit",
        type = "logical",
        default = FALSE,
        help = "Flag on whether to re-run the stan models even if already exist [Defaults to FALSE]", 
        dest = 'refit'
    ),
    optparse::make_option(
        "--run-comm-analysis",
        type = "logical",
        default = FALSE,
        help = "Flag on whether to run the Stan GLM model for community-level suppression [Defaults to FALSE] ",
        dest = 'run.comm.analysis'
    ),
    optparse::make_option(
        "--viremic-viral-load",
        type = "numeric",
        default = 1000,
        help = "Duration of acute phase, in years [Defaults to 2 months]", 
        dest = 'viremic.viral.load'
    ),
    optparse::make_option(
        "--vl-detectable",
        type = "numeric",
        default = 150,
        help = "Duration of acute phase, in years [Defaults to 2 months]", 
        dest = 'vl.detectable'
    ),
    optparse::make_option(
        "--outdir",
        type = "character",
        default='/home/andrea/HPC/ab1820/home/projects/2022/longvl',
        help = "Path to output directory where to store results [Defaults to my directory]", 
        dest = 'out.dir.prefix'
    ),
    optparse::make_option(
        "--round",
        type = "numeric",
        default = c(16,17,18,19),
        help = "Path to output directory where to store results [Defaults to my directory]", 
        dest = 'round'
    )
)

args <-  optparse::parse_args(optparse::OptionParser(option_list = option_list))
print(args)

################
#    PATHS     #
################

parallelise <- FALSE

gitdir <- here()
source(file.path(gitdir,'R/paths.R'))

file.exists(
    # out.dir,
    gitdir.stan,
    path.hivstatusvl.r1520,
    path.quest_r1519_221207 
) |> all() |> stopifnot()


# load args into env variables
with(args,{
    VIREMIC_VIRAL_LOAD  <<- viremic.viral.load
    VL_DETECTABLE  <<- vl.detectable
    stopifnot( dir.exists(out.dir.prefix))
    out.dir <<- file.path(out.dir.prefix)

    vl.out.dir <<- file.path(out.dir.prefix, paste0('vl_', VIREMIC_VIRAL_LOAD) )
    out.dir.testvl <<- file.path(vl.out.dir, 'test_vs_suppression')
})

dir.create(vl.out.dir, showWarnings = FALSE)
dir.create(out.dir.testvl, showWarnings = FALSE)

################
#    HELPERS   #
################

sapply(R_scripts, source) |> invisible()

source(file.path(gitdir.functions, 'preprocessing_helpers.R'))
# source(file.path(gitdir.functions,'base_utilities.R') )
source(file.path(gitdir.functions, 'vl_testing_helpers.R'))
source(file.path(gitdir.functions,'phsc_vl_helpers.R'))

# set up parallel backend if want
if(parallelise)
    source(file.path(gitdir, 'parallel_backend.R'))

naturemed_reqs() 

################
#     MAIN     #
################

# check study round exists
stopifnot(all(args$round %in% 16:19))

# get status&suppression ; testing data; and new diagnoses
dstatus <- get.dall(path.hivstatusvl.r1520, make_flowchart=FALSE)
dtesting <- readRDS(path.processed.testing.r1519) |> process.hiv.testing.data(rounds=18:19)
dtesting[, I := .I]

newdiagnoses <- readRDS(path.processed.hivstatus.r0920) |>
    subset(ROUND %in% args$round & ROUND == ROUND_FP, select=c('STUDY_ID', 'ROUND'))
newdiagnoses[, NEW_DIAGN:=TRUE]
stopifnot(uniqueN(newdiagnoses) == nrow(newdiagnoses))

# merge and compare hiv status dataset with testing dataset
# 3 cases in which the COMM_NUM differs: just take it from dstatus
dall2 <- merge( dstatus[ROUND %in% c(18,19)], dtesting, all.x=TRUE, all.y=TRUE)
dall2[, {
    N.hivstatus <- sum(!is.na(HIV_STATUS));
    N.nohivstatus <- sum(is.na(HIV_STATUS));
    N.testing <- sum(!is.na(I));
    N.notesting <- sum(is.na(I));
    ids.notesting <<- STUDY_ID[is.na(I)]
    ids.nohivstatus <<- STUDY_ID[is.na(HIV_STATUS)]
    sprintf('%d out of %d inds (%.2f%%) in the testing data do not appear in the vl suppression data\n',
        N.nohivstatus ,N.testing, N.nohivstatus/N.testing*100) |> cat()
    sprintf('%d out of %d inds (%.2f%%) in the vl suppression data do not appear in the testing data\n',
        N.notesting , N.hivstatus, N.notesting/N.hivstatus*100)|> cat()
}]
if(0){ uniqueN(ids.notesting);uniqueN(ids.nohivstatus)}

# only include those that appear in both datasets
dall2 <- dall2[ !is.na(HIV_STATUS) & !is.na(I)]
dall2[, I := NULL]

# newdiagnoses: check if all of them are in dall2: YES
setkey(dall2, "STUDY_ID", "ROUND")
dall2[newdiagnoses][, stopifnot(.N == nrow(newdiagnoses))]
dall2 <- merge(dall2, newdiagnoses, c("STUDY_ID", "ROUND"), all.x=TRUE)
dall2[is.na(NEW_DIAGN), NEW_DIAGN := FALSE]

# can "safely" assume that NA in TEST_YEAR AGO are never tested
dall2[is.na(TEST_YEAR_AGO), table(TEST_EVER)]

####################
# PLOTS AND TABLES #
####################

# suppression vs average time to last test:
# _________________________________________
catn("Do (aware) HIV positive participants test?")
with(dall2,
    table(TEST_YEAR_AGO, HIV_STATUS) |> 
        proportions(margin=1) |> 
        addmargins(2) |> 
        round(2) |> kable(col.names=c("Hiv -", "HIV +", 'Total')) )

p_TvsS <- unique(dall2[, .(FC, COMM_NUM)]) |> 
    merge( make.testing.frequency.metric.among.negatives(dall2), by='COMM_NUM') |>
    merge( make.hivsuppression.metric.among.positives(dall2), by=c('COMM_NUM', 'SEX', 'ROUND')
) |>  
    prettify_labels() |>
    ggplot(aes(x=TEST_YA_MEAN, y=VL_SUPP_MEAN, color=FC)) +
    # geom_point() +
    geom_text(aes(label=COMM_NUM), size=2) +
    facet_grid(SEX_LAB~ROUND_LAB) +
    theme(legend.position='bottom') +
    labs(
        x='Average time since last test among HIV negatives',
        y="Prevalence of viral suppression",
        color="Community type" 
    ) +
    NULL
filename <- paste0(out.dir.testvl, 'metrictestvsmetricsupp_by_comm_round_sex.pdf')
ggsave_nature(p=p_TvsS, filename=filename, LALA="", w=18, h=15)


# even among hiv negatives, ~ 90% of people have tested,and recently
dall2[HIV_STATUS == 0, {
    table(TEST_EVER, useNA = 'ifany') |> proportions() |> round(2) |>
        kable(caption="Ever tested (among HIV -)") |> print()
    table(TEST_YEAR_AGO < 2, useNA='ifany') |> proportions() |> round(2)|> 
        kable(caption="Tested at most one year ago (among HIV -)") |> print()
    NULL
}]

p_ago <- plot.testing.answers.proportions(DT=dall2,var="TEST_YEAR_AGO", group='negatives',chi2=FALSE)
filename <- paste0(out.dir.testvl, 'hist_answerstestago_by_comm_round_sex.pdf')
ggsave_nature(p=p_ago, filename=filename, LALA="", w=21, h=16)

p_ev <- plot.testing.answers.proportions(DT=dall2,var="TEST_EVER", group='negatives',chi2=FALSE)
filename <- paste0(out.dir.testvl, 'hist_answerstestever_by_comm_round_sex.pdf')
ggsave_nature(p=p_ev, filename=filename, LALA="", w=21, h=16)

p_ago2 <- plot.testing.answers.proportions(DT=dall2,var="TEST_YEAR_AGO", group='newdiagnoses', chi2=FALSE)
filename <- paste0(out.dir.testvl, 'hist_answerstestago2_by_comm_round_sex.pdf')
ggsave_nature(p=p_ago2, filename=filename, LALA="", w=21, h=16)

p_ev2 <- plot.testing.answers.proportions(DT=dall2,var="TEST_EVER", group='newdiagnoses', chi2=FALSE)
filename <- paste0(out.dir.testvl, 'hist_answerstestever2_by_comm_round_sex.pdf')
ggsave_nature(p=p_ev, filename=filename, LALA="", w=21, h=16)

tmp <- dall2[, .(ROUND, SEX, AGEYRS, COMM_NUM, HIV_STATUS, HIV_AND_VL, TEST_LAST_YEAR=TEST_YEAR_AGO < 1)]
tmp[, TEST_LAST_YEAR := factor(TEST_LAST_YEAR, levels=c("TRUE", "FALSE"))]

for( bool in c(TRUE, FALSE))
{
    bool = TRUE
    p1 <- plot.testing.answers.proportions(
        dall2,
        var="TEST_YEAR_AGO",
        group='negatives',
        percentage=bool)
    p <- plot.testing.answers.proportions(
        tmp,
        var="TEST_LAST_YEAR",
        group='negatives',
        percentage=bool)
    p
}

p <- plot.suppression.vs.testing.answers.proportions(
    DT=dall2, 
    .group = 'negatives_and_newdiagnoses',
    type = 'evertest')
filename <- paste0(out.dir.testvl, 'scatt_pevertest_vs_psupp_comm_by_sex.pdf')
ggsave_nature(p=p, filename=filename, LALA="", w=18, h=15)

p <- plot.suppression.vs.testing.answers.proportions(
    DT=dall2,
    .group='negatives_and_newdiagnoses',
    type='test1ya')
filename <- paste0(out.dir.testvl, 'scatt_ptest1ya_vs_psupp_comm_by_sex.pdf')
ggsave_nature(p=p, filename=filename, LALA="", w=18, h=15)
