library(data.table)
library(ggplot2)
library(patchwork)
library(readxl)

#########
# PATHS #
#########



indir.deepdata <- '/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live/'
indir.deepdata.r1520 <- file.path(indir.deepdata, 'RCCS_R15_R20')
indir.deepdata.r1518 <- file.path(indir.deepdata, 'RCCS_R15_R18')
indir.deepdata.r0914 <- file.path(indir.deepdata, 'RCCS_R9_R14')

path.quest_r1520_220830 <- file.path(indir.deepdata.r1520, 'Quest_R015_R020_VOIs_August302022.csv')
path.quest_r1519_221207 <- file.path(indir.deepdata.r1520, 'quest_R15_R19_VoIs_Dec072022.csv')

indir.deepsequence_analyses <- '/home/andrea/HPC/project/ratmann_xiaoyue_jrssc2022_analyses/live/PANGEA2_RCCS1519_UVRI'
path.community.types <- file.path(indir.deepsequence_analyses, 'community_names.csv')


env <- new.env()
load(path.community.types, envir=env)
ls(env)
env$comgps

# outdir
outdir.misc <- '/home/andrea/HPC/ab1820/home/projects/2022/longvl/misc'
dir.create( outdir.misc, showWarnings = FALSE )

###########
# HELPERS #
###########

catn <- function(x) cat('\n----', x, '----\n')

round2numeric <- function(x) 
{
    x <- gsub('S$','.5',x)
    as.numeric(gsub( '[A-z]', '', x))
}

empty2NA <- function(DT)
{
    toNA <- function(vec) 
    {
        vec[vec==""] <- NA
        return(vec)
    }
    DT[, lapply(.SD, toNA), ]
}

question_dict <- list(
    rhivever2= 'Have you ever received your HIV results from anywhere?',
    RHIVEVER2= 'Have you ever received your HIV results from anywhere?',
    hivperiod2= 'How long ago did you last receive your HIV results?(years)',
    HIVPERIOD2= 'How long ago did you last receive your HIV results?(years)',
    hivrslt2= 'What was the results of this last HIV test?',
    HIVRSLT2= 'What was the results of this last HIV test?',
    oghivts2= 'From whom did you receive the test?',
    OGHIVTS2= 'From whom did you receive the test?'
)

palettes <- list( 
    comm2 = c("#9C964A", "#85D4E3","#FF9933")
)

########
# MAIN #
########

# Transform xlsx to csv if not done already
if(! file.exists(path.quest_r1519_221207))
{
    xlsx.path <- gsub('.csv','.xlsx',path.quest_r1519_221207)
    read_xlsx(xlsx.path) |> as.data.table() |> fwrite(path.quest_r1519_221207)
}


# load hiv testing related stuff
load.hiv.testing.data <- function(path, print_statements=FALSE, make_plots=FALSE, verbose = FALSE)
{

    if(verbose==TRUE)
    {
        tmp <- path |> fread()
        hiv_cols <- grep( 'hiv', names(tmp), value=TRUE )
        hiv_tmp <- subset(tmp, select=c('round', hiv_cols))
        cat('Ever tested for HIV?\n')
        hiv_tmp[, round(100*table(round, rhivever)/.N, 2) ,]
        cat('When last?\n')
        hiv_tmp[, round(100*table(round, hivperiod)/.N, 2) ,]
        cat('Which result?\n')
        hiv_tmp[, round(100*table(round, hivrslt)/.N, 2) ,]
        cat('Which oghivts1?\n')
        hiv_tmp[, round(100*table(round, oghivts1)/.N, 2) ,]
    }

    cols_informing_hivtesting <- c('rhivever', 'hivperiod', 'hivrslt', 'oghivts1')
    cols_identifying_person_round <- c('study_id', 'round')

    # load data
    cols <- c(cols_informing_hivtesting, cols_identifying_person_round)
    dtesting <- fread(path, select=cols)

    # understand labels
    dictionaries <- list(
        rhivever = c(
            '1'='1.Yes',
            '2'='2.No',
            '8'='8.Not Applicable',
            '9'='9.NR'
        ),
        hivperiod= c(
            '1'= '1.(<1)',
            '2'= '2.(1-2)',
            '3'= '3.(3-4)',
            '4'= '4.(>4)',
            '8'= 'NA?'),
        hivrslt = c(
            '1'='1.Negative',
            '2'='2.Positive',
            '3'='3.Indeterminate',
            '7'="7.Don't know",
            '8'='8. ?',
            '9'='9.NR'),
        oghivts1 = c(
            '1'='1. Rakai Project',
            '2'='2. RAIN',
            '3'='3. Other NGO (specify)',
            '4'='4. Gov. Doctor/nurse',
            '5'='5. Priv. Doctor/nurse/HW',
            '6'='6. Other (specify)',
            '7'='7.97.98.998. ?',
            '88'='88. No additional response',
            '97'='7.97.98.998. ?',
            '98'='7.97.98.998. ?',
            '998'='7.97.98.998. ?'
        )
    )

    # check no doubles
    idcols <- c('study_id', 'round')
    dtesting[, .N, by=idcols][, stopifnot(all(N==1))]

    # rename with labels
    dtesting[, rhivever2  := dictionaries$rhivever[as.character(rhivever)] ]
    dtesting[, hivperiod2 := dictionaries$hivperiod[as.character(hivperiod)] ]
    dtesting[, hivrslt2   := dictionaries$hivrslt[as.character(hivrslt)] ]
    dtesting[, oghivts2   := dictionaries$oghivts1[as.character(oghivts1)] ]

    # translate labels:
    dtesting[, table(rhivever, rhivever2, useNA='ifany')]
    dtesting[, table(hivperiod, hivperiod2, useNA='ifany')]
    dtesting[, table(hivrslt, hivrslt2, useNA='ifany')]
    dtesting[, table(oghivts1, oghivts2, useNA='ifany')]


    if(print_statements)
    {
        .cat <- function(x) cat('\n',x, '\n')
        .pr <- function(x) { table(x, round, useNA='ifany') |> knitr::kable() |> print()}

        .cat('Have you ever received your HIV results from anywhere')
        dtesting[, table(rhivever2,round, useNA='ifany') |> print() ]
        .cat('How long ago did you last receive your HIV results?(years)')
        dtesting[, table(hivperiod2, round, useNA='ifany') |> print() ]
        .cat('What was the results of this last HIV test?')
        dtesting[, table(hivrslt2, round, useNA='ifany') |> print() ]
        .cat('From whom did you receive the test?')
        dtesting[, table(oghivts2, round,useNA='ifany') |> print() ]
    }

    if(make_plots)
    {
        .plot <- function(x, q)
        {
            DT <- as.data.table(x) 
            # DT[, round := gsub("R0|S", "", round)  ]
            question <- names(DT)[1]
            p <- ggplot(DT, aes_string(x="round", fill=question)) +
                geom_col(aes(y=N)) + 
                theme_bw() +
                theme(legend.position='bottom') +
                scale_y_continuous(expand=expansion(mult=c(0,.1)))+
                labs(x='Rounds', y='Number of answers', title=q, fill='') +
                NULL
            force(p)
            p
        }
 
        plots <- list(
            p_ever = dtesting[, table(rhivever2,round, useNA='ifany') ] |>
                .plot(q='Have you ever received your HIV results from anywhere?'),

            p_prd  = dtesting[rhivever==1, table(hivperiod2, round, useNA='ifany')] |>
                .plot(q='How long ago did you last receive your HIV results?(years)'),

            p_rslt = dtesting[rhivever==1, table(hivrslt2, round, useNA='ifany') ] |>
                .plot(q='What was the results of this last HIV test?'),

            p_where = dtesting[rhivever==1, table(oghivts2, round,useNA='ifany') ] |>
                .plot(q='From whom did you receive the test?')
        )

        p <- (plots[[1]] + plots[[2]])/(plots[[3]] + plots[[4]])
        filename = file.path(outdir.misc, 'testing_data_r1519.png')
        ggsave(p, filename=filename, height=11, width=14)
        plots <<- plots
        p_testing <<- p
    }

    # check output before returning
    dtesting[, .N, by=idcols][, stopifnot(all(N==1))]

    return(dtesting)
}

dtesting <- load.hiv.testing.data( path.quest_r1519_221207, make_plots=TRUE)


data <- fread(path.quest_r1519_221207)
cols <- data |> names()
#   [1] "study_id" "round"    "comm_num" "curr_id"  "intdate"  "birthmo" 
#   [7] "birthyr"  "birthdat" "ageyrs"   "agecomp"  "sex"      "occup1"  
#  [13] "occup2"   "pregnow"  "evermarr" "currmarr" "polymar"  "currrltn"
#  [19] "pregtest" "eversex"  "ag1stsex" "sexyear"  "sexp1yr"  "sexp1out"
#  [25] "sexpever" "rltn1"    "rltn2"    "rltn3"    "rltn4"    "sxm1"    
#  [31] "rnyrcon1" "rnyrcon2" "rnyrcon3" "rnyrcon4" "rltongo1" "rltongo2"
#  [37] "rltongo3" "rltongo4" "hivquest" "rhivrlst" "oghivts1" "circum"  
#  [43] "rldyslt1" "rlwkslt1" "rlmoslt1" "rlyrslt1" "rltnage1" "rltnyrs1"
#  [49] "rltnlst1" "rltnhh1"  "rldyslt2" "rlwkslt2" "rlmoslt2" "rlyrslt2"
#  [55] "rltnage2" "rltnyrs2" "rltnlst2" "rltnhh2"  "rldyslt3" "rlwkslt3"
#  [61] "rlmoslt3" "rlyrslt3" "rltnage3" "rltnyrs3" "rltnlst3" "rltnhh3" 
#  [67] "rldyslt4" "rlwkslt4" "rlmoslt4" "rlyrslt4" "rltnage4" "rltnyrs4"
#  [73] "rltnlst4" "rltnhh4"  "days1"    "weeks1"   "months1"  "years1"  
#  [79] "days2"    "weeks2"   "months2"  "years2"   "days3"    "weeks3"  
#  [85] "months3"  "years3"   "days4"    "weeks4"   "months4"  "years4"  
#  [91] "rltncm1"  "rltncm2"  "rltncm3"  "rltncm4"  "occup11"  "occup21" 
#  [97] "occup12"  "occup22"  "occup13"  "occup23" 
#  [1] "occup14"   "occup24"   "id"        "sexp1yrb" 
#  [6] "cuarvmed"  "cuseptmed" "rhivever"  "pregnow2"  "arvmed"   
# [11] "rltncom1"  "rltncom2"  "rltncom3"  "rltncom4"  "hivrslt"  
# [16] "hivperiod" "hivcare"   "religion"  "educate"   "educyrs"  

rltn_data <- data |> subset(select=cols[cols %like% 'rltn'])

cols <- c('study_id',  'sex',  'birthdat', 'ageyrs', 'round', 'comm_num')
data |> subset(select=cols)



if(0)   # RAPID TESTS FROM ROUND 19 ONLY...
{
    cols <- c('study_id', 'round', 'rapid')
    tmp <- file.path( indir.deepdata.r1520, "flowR19_VOIs.dta") |>
        haven::read_dta() |> as.data.table() |> subset(select=cols)
    tmp[!is.na(study_id), table(round, rapid)]
}

catn("** study datasets from R15-R20 **")

            
catn('get positives?')

hivstatus <- file.path(indir.deepdata.r1520, 'all_participants_hivstatus_vl_220729.csv') |>
    fread(select=c('STUDY_ID', 'ROUND', 'HIV_STATUS')) |> unique()


catn('get negatives')

negatives <- file.path(indir.deepdata.r1520, 'R016_R020_Data_for_HIVnegatives.csv') |>
    fread(select=c('study_id', 'round')) |> unique()
names(negatives) <- toupper(names(negatives)) 
negatives[, ROUND:= round2numeric(ROUND)]
check_negs_all_in_hivstatus <- merge(negatives, hivstatus)[, sum(HIV_STATUS==0) == nrow(negatives)]
stopifnot(check_negs_all_in_hivstatus)


catn("** study datasets from R09-R14 **")

list.files( indir.deepdata.r0914, pattern="")
history <- file.path(indir.deepdata.r0914, "HIV_R09_R14.csv") |>
    fread(select=c('study_id', 'round', 'hiv', 'lastnegv', 'firstpos_diagnosis_vis')) |> 
    empty2NA() |> 
    subset(!is.na(hiv))

history[, table(firstpos_diagnosis_vis, useNA = 'ifany')]
cols <- c('round', 'firstpos_diagnosis_vis')
history[,  (cols) := lapply(.SD, round2numeric), .SDcols=cols]
names(history) <- toupper(names(history))

# what is the definition of firstpos diagnosis visit? 
idx <- history[ HIV == 'P', unique(STUDY_ID) ]
history[STUDY_ID %in% idx]
idx <- history[, any(HIV == 'P') & any(HIV=='N'), by=STUDY_ID][V1==TRUE, STUDY_ID]
history[STUDY_ID %in% idx, sum(!is.na(FIRSTPOS_DIAGNOSIS_VIS))]

catn("** merge **")

history_as_status <- subset(history, select=c('STUDY_ID', 'ROUND', 'HIV'))
history_as_status[, `:=` (HIV_STATUS=fifelse(HIV=='P', yes=1, no=0), HIV=NULL)]

hivstatus2 <- rbind(history_as_status, hivstatus)

# get first positives
setkey(hivstatus2, STUDY_ID, ROUND)
hivstatus2[, `:=` (
    ROUND_FP = ROUND[min(which(HIV_STATUS == 1))],
    ANY_NEG = any(HIV_STATUS == 0)
    ) , by="STUDY_ID"]
first_positives <- hivstatus2[ROUND_FP == ROUND,]
cat("The number of first positive test in the RCCS by round is:\n")
first_positives[, table(ROUND, ANY_NEG) ] |> knitr::kable(col.names = )

# check consistency with history
# in ~ 250 cases we erroneously impute 'no previous' negatives, because those were prior to round 8
dchecks <- merge(first_positives, history, all.x=TRUE)
dchecks[!is.na(FIRSTPOS_DIAGNOSIS_VIS) & ANY_NEG == FALSE, table(LASTNEGV, ROUND_FP)]
dchecks[FIRSTPOS_DIAGNOSIS_VIS != ROUND_FP, table(ROUND_FP,FIRSTPOS_DIAGNOSIS_VIS) ]

# there are about 4000 first positive cases in rounds 16 to 20
first_positives_1620 <- subset(first_positives, ROUND_FP > 15.5)
first_positives_1620

# let's check their answers
tmp <- copy(dtesting)
names(tmp) <- toupper(names(tmp)) 
tmp[, ROUND:= round2numeric(ROUND)]
tmp1 <- merge(first_positives_1620, tmp)
names(tmp1) <- tolower(names(tmp1)) 

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


# can we get it by community? 

catn("Load community types")
dcommtype <- fread(path.community.types)[, 
    COMM_TYPE := fcase(
        COMM_NUM_A %like% '^f', 'fishing',
        COMM_NUM_A %like% '^a', 'agrarian',
        COMM_NUM_A %like% '^t', 'trading') ][,
    `:=` (COMM_NUM = COMM_NUM_RAW, COMM_NUM_RAW=NULL, COMM_NUM_A=NULL) ]


catn('load communities')
dcomms <- subset(data, select=c('study_id', 'round', 'comm_num')) |> unique()
names(dcomms) <- toupper(names(dcomms))
dcomms[, ROUND := round2numeric(ROUND)]
merge(dcomms, dcommtype)


make.plots.community <- function(dtest, bytype=FALSE)
{
    # dtest <- copy(dtesting); bytype=TRUE
    names(dtest) <- toupper(names(dtest)) 
    if( !is.numeric(dtest$ROUND))
    {
        dtest[, ROUND := round2numeric(ROUND)]
    }
    dtest <- merge(dcomms, dtest, by=c('STUDY_ID', 'ROUND'))
    dtest <- merge(dtest, dcommtype, by=c('COMM_NUM'))

    aggregation_var <- fifelse(bytype==TRUE, yes="COMM_TYPE", no="COMM_NUM")

    dtest_aggr_comm <- melt(dtest,
        id.vars=c(aggregation_var, 'ROUND'),
        measure.vars= c('RHIVEVER2', 'HIVPERIOD2', 'HIVRSLT2', 'OGHIVTS2'),
        value.name = 'VALUE', variable.name = 'VARIABLE')
    dtest_aggr_comm[, N := .N , c(aggregation_var,'ROUND', 'VARIABLE', 'VALUE')]
    dtest_aggr_comm <- unique(dtest_aggr_comm)
    dtest_aggr_comm[, TOT := sum(N), by=c(aggregation_var, 'ROUND', 'VARIABLE')]
    dtest_aggr_comm[, PROP := N / TOT, by=c(aggregation_var, 'ROUND', 'VARIABLE')]

    plot.question.results.by.community <- function(var, DT)
    {
        var <- as.character(var)
        cat(question_dict[[var]], '\n')
        dtest_aggr_comm2 <- DT |> subset(VARIABLE == var & ROUND > 15.5)
        dtest_aggr_comm2[, AGGVAR:=lapply( .SD, as.factor), .SDcols=aggregation_var]

        # plot settings
        leg_nrow = fifelse(bytype==TRUE, yes=1, no=2)


        p <- dtest_aggr_comm2 |> 
            ggplot(aes(x=ROUND, color=AGGVAR, y=PROP)) +
            geom_line() +
            facet_grid(~ VALUE) +
            theme_bw() +
            theme(legend.position='bottom') + {
                if(bytype) scale_color_manual(values=palettes$comm2) 
            } +
            scale_y_continuous(labels=scales::percent) + 
            labs(title = question_dict[[var]], x='RCCS Round', y='Proportion of answers', color='community') +
            guides(color=guide_legend(nrow=leg_nrow,byrow=TRUE))
            

        force(p)
        return(p)
    }

    plots_comm <- lapply(unique(dtest_aggr_comm$VARIABLE), plot.question.results.by.community, DT=dtest_aggr_comm)
    p_aggregated <- ((plots_comm[[1]] + plots_comm[[2]])/(plots_comm[[3]] + plots_comm[[4]]) ) & theme(legend.position ='bottom') 
    p_aggregated <- p_aggregated + plot_layout(guides = "collect")
    p_aggregated
}

p1 <- make.plots.community(dtesting)
filename = file.path(outdir.misc, 'testing_data_r1619_bycommunity.pdf')
ggsave(p1, filename=filename, height=11, width=18)

p1t <- make.plots.community(dtesting, bytype=TRUE)
filename = file.path(outdir.misc, 'testing_data_r1619_bycommunitytype.pdf')
ggsave(p1t, filename=filename, height=11, width=18)

dtest_firstpos <- dtesting[, cbind(.SD, ROUND=round2numeric(round)), .SDcols=!patterns('round') ] |> 
    merge(first_positives, by.x=c('study_id', 'ROUND'), by.y=c('STUDY_ID', 'ROUND'))
p2 <- make.plots.community(dtest_firstpos)
filename = file.path(outdir.misc, 'testing_data_r1619_firstpos_bycommunity.pdf')
ggsave(p2, filename=filename, height=11, width=18)

p2t <- make.plots.community(dtest_firstpos, bytype=TRUE)
filename = file.path(outdir.misc, 'testing_data_r1619_firstpos_bycommunitytype.pdf')
ggsave(p2t, filename=filename, height=11, width=18)

dtest_firstpos[hivrslt2 %like% 'Positive', table(oghivts2)/.N] |>
    knitr::kable()

# get participants saying they had a previous Rakai positive test
idx <- dtest_firstpos[hivrslt2 %like% 'Positive' & oghivts2 %like% 'Rakai', .(study_id, round=ROUND)]
names(idx) <- toupper(names(idx))

history[study_id %in% idx$study_id]

dtesting[study_id %in% idx$study_id]

merge(idx[, DUMMY:='HELLO'], hivstatus, all.y=TRUE)[STUDY_ID %in% idx$STUDY_ID][, list(
    ROUND[which(DUMMY == 'HELLO')],
    ROUND[which(HIV_STATUS ==1)[1]]
), by=STUDY_ID ]
hivstatus[STUDY_ID %in% idx$study_id]
