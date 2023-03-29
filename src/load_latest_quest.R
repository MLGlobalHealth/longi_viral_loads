#############
# OBJECTIVE #
#############

"
Make individual level dataset with:
- viral load 
- sero-status
- testing history (ie. ever tested positive before or not)
"

#############
# LIBRARIES #
#############

library(data.table)
library(ggplot2)
library(patchwork)
library(readxl)
library(here)
library(knitr)
# library(kableExtra)
library(optparse)

###########
# OPTIONS #
###########

option_list <- list(
    make_option(
        "--plots",
        type = "logical",
        default = FALSE,
        help = "Whether to store output plots [defaults to FALSE]",
        dest= "make_plots"
    )
)
args <- parse_args(OptionParser(option_list = option_list))


#########
# PATHS #
#########

gitdir <- here()
source(file.path(gitdir, "paths.R"))

gitdir.functions <- file.path(gitdir, 'functions')

# outdir
outdir.misc <- '/home/andrea/HPC/ab1820/home/projects/2022/longvl/misc'
dir.create( outdir.misc, showWarnings = FALSE )

# check files exist
file.exists(path.hivstatusvl.r1520) |> stopifnot()
file.exists(path.community.types) |> stopifnot()
file.exists(path.quest_r1519_221207) |> stopifnot()
# file.exists(path.negatives.r1520) |> stopifnot() # not necessary
file.exists(path.hivres.r0914) |> stopifnot()

###########
# HELPERS #
###########

source(file.path(gitdir.functions, 'base_utilities.R'))
source(file.path(gitdir.functions, 'preprocessing_helpers.R'))

########
# MAIN #
########

# Transform xlsx to csv if not done already
# _________________________________________

if(! file.exists(path.quest_r1519_221207))
{
    xlsx.path <- gsub('.csv','.xlsx',path.quest_r1519_221207)
    read_xlsx(xlsx.path) |> as.data.table() |> fwrite(path.quest_r1519_221207)
}

# Load data
# _________

catn("** study datasets from R15-R20 **")

# testing data from 15 and 19: for HIV + and -
dtesting <- load.hiv.testing.data( path.quest_r1519_221207, make_plots=TRUE, verbose=TRUE)
dtesting[, {
    cat('\n', question_dict[['RHIVEVER2']] ,'\n')
    table(round, rhivever2) |> knitr::kable() |> print() ;
    cat('\n', question_dict[['HIVPERIOD2']] ,'\n')
    table(round, hivperiod2) |> knitr::kable() |> print();
    cat('\n', question_dict[['HIVRSLT2']] ,'\n')
    table(round, hivrslt2) |> knitr::kable() |> print();
    cat('\n', question_dict[['OGHIVTS2']] ,'\n')
    table(round, oghivts1) |> knitr::kable() |> print();
}]
# dtesting[, .SD, .SDcols = ! names(dtesting) %like% '2$']

catn('get positives from R15-20 and R09-14 and then merge')
cols <- c('STUDY_ID', 'ROUND', 'HIV_STATUS')
hivstatus <- fread(path.hivstatusvl.r1520, select=cols) |> unique()
history <- get.testing.history(path.hivres.r0914)
# history[, table(LASTNEGV, useNA = 'ifany')]
# transform history to have same column names as hivstatus, and then bind
history_as_status <- subset(history, select=c('STUDY_ID', 'ROUND', 'HIV'))
history_as_status[, `:=` (HIV_STATUS=fifelse(HIV=='P', yes=1, no=0), HIV=NULL)]
hivstatus <- rbind(history_as_status, hivstatus)
hivstatus[, .N, by=c('STUDY_ID', 'ROUND')][, stopifnot(all(N == 1))]

# get first positives: defined as first positive test in the rounds of 9 to 20
setkey(hivstatus, STUDY_ID, ROUND)
hivstatus[, `:=` (
    ROUND_FP = suppressWarnings(ROUND[min(which(HIV_STATUS == 1))]),
    ANY_NEG = any(HIV_STATUS == 0)
    ) , by="STUDY_ID"]
first_positives <- hivstatus[ROUND_FP == ROUND,]
first_positives[, table(ROUND, ANY_NEG) ] |> 
    kable(format='markdown', caption='Number of first positives tests\n(grouped by any negatives before)')

catn('load communities, their type, and merge together')

dcommtype <- fread(path.community.types)[, 
    COMM_TYPE := commcode2commtype(COMM_NUM_A) ][,
    `:=` (COMM_NUM = COMM_NUM_RAW, COMM_NUM_RAW=NULL, COMM_NUM_A=NULL) ]

dcomms <- subset(data, select=c('study_id', 'round', 'comm_num')) |> unique()
names(dcomms) <- toupper(names(dcomms))
dcomms[, ROUND := round2numeric(ROUND)]
merge(dcomms, dcommtype)

# there are about 4000 first positive cases in rounds 16 to 20
first_positives_1620 <- subset(first_positives, ROUND_FP > 15.5)

####################
# Checks and plots #
####################

# check consistency with history
# in ~ 250 cases we erroneously impute 'no previous' negatives, because those were prior to round 8
# however, nothing wrong from R15 onwards
dchecks <- merge(first_positives, history, all.x=TRUE)
dchecks[!is.na(FIRSTPOS_DIAGNOSIS_VIS) & ANY_NEG == FALSE, table(LASTNEGV, ROUND_FP)] |> 
    kable(format='markdown',
        caption='"First positives" with a previous positive date before R14')

if(0)
{
    # I have the testing data for every person in hivstatus R15-19
    # conversely, I have the hivstatus for every person in dtesting R15-19
    hivstatus.1519 <- hivstatus[ROUND %between% c(15, 19), .(STUDY_ID, ROUND, STAT=TRUE)]
    dtesting.1519 <- dtesting[ROUND %between% c(15, 19), .(STUDY_ID, ROUND, TEST=TRUE)]
    merge(hivstatus.1519,dtesting.1519, by=c('STUDY_ID', 'ROUND'),all.x=TRUE, all.y=TRUE)[ ,table(STAT, TEST)]
}

# testing plots
# _____________

tmp <- copy(dtesting)
names(tmp) <- toupper(names(tmp)) 
tmp[, ROUND:= round2numeric(ROUND)]
tmp1 <- merge(first_positives_1620, tmp)
names(tmp1) <- tolower(names(tmp1)) 

dtest_firstpos <- dtesting[, cbind(.SD, ROUND=round2numeric(round)), .SDcols=!patterns('round') ] |> 
    merge(first_positives, by.x=c('study_id', 'ROUND'), by.y=c('STUDY_ID', 'ROUND'))

if(args$make_plots)
{
    # by round
    plots <- list(
        p_ever = tmp1[, table(rhivever2,round, useNA='ifany') ] |>
            .plot(q='Have you ever received your HIV results from anywhere?'),

        p_prd  = tmp1[rhivever==1, table(hivperiod2, round, useNA='ifany')] |>
            .plot(q='How long ago did you last receive your HIV results?(years)'),

        p_rslt = tmp1[rhivever==1, table(hivrslt2, round, useNA='ifany') ] |>
            .plot(q='What was the results of this last HIV test?'),

        p_where = tmp1[rhivever==1, table(oghivts2, round,useNA='ifany') ] |>
            .plot(q='From whom did you receive the test?')
    )
    p <- (plots[[1]] + plots[[2]])/(plots[[3]] + plots[[4]])
    filename = file.path(outdir.misc, 'testing_data_r1619_firstpositives.png')
    ggsave(p, filename=filename, height=11, width=14)
}
# by community and round
if(args$make_plots)
{
    p1 <- make.testing.plots.community(dtesting)
    filename = file.path(outdir.misc, 'testing_data_r1619_bycommunity.pdf')
    ggsave(p1, filename=filename, height=11, width=18)

    p1t <- make.testing.plots.community(dtesting, bytype=TRUE)
    filename = file.path(outdir.misc, 'testing_data_r1619_bycommunitytype.pdf')
    ggsave(p1t, filename=filename, height=11, width=18)

    p2 <- make.testing.plots.community(dtest_firstpos)
    filename = file.path(outdir.misc, 'testing_data_r1619_firstpos_bycommunity.pdf')
    ggsave(p2, filename=filename, height=11, width=18)

    p2t <- make.testing.plots.community(dtest_firstpos, bytype=TRUE)
    filename = file.path(outdir.misc, 'testing_data_r1619_firstpos_bycommunitytype.pdf')
    ggsave(p2t, filename=filename, height=11, width=18)
}

dtest_firstpos[hivrslt2 %like% 'Positive',
    table(oghivts2, useNA = 'ifany') |> proportions() |> kable(format='markdown') ]

# get participants saying they had a previous Rakai positive test
# they do not appear in my testing or history dataset
idx <- dtest_firstpos[hivrslt2 %like% 'Positive' & oghivts2 %like% 'Rakai',
    .(STUDY_ID=study_id, round=ROUND)]

sprintf('There are %d individuals reporting a positive Rakai test:
    - %d in the are in the `hivstatus` dataset
    ', nrow(idx))

hivstatus[STUDY_ID %in% idx$study_id]

merge(idx[, DUMMY:='DUMMY'], hivstatus, all.y=TRUE)[STUDY_ID %in% idx$STUDY_ID][, list(
    ROUND[which(DUMMY == 'DUMMY')],
    ROUND[which(HIV_STATUS ==1)[1]]
), by=STUDY_ID ]


######################
# Merge all and save #
######################

# prettify testing
names(dtesting) <- toupper(names(dtesting)) 
dtesting[, ROUND := round2numeric(ROUND)]
dtesting <- dtesting[, .SD, .SDcols = ! names(dtesting) %like% '2$'] |> unique()

# ending in: "RCCS_processed_participants_hivstatus_230328.rds")
savefile(data=hivstatus, filename=path.processed.hivstatus.r0920)
savefile(data=dtesting, filename=path.processed.testing.r1519)
