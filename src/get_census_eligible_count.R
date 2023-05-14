{
    library(data.table)
    library(dplyr)
    library(ggplot2)
    library(scales)
    library(lubridate)
    library(haven)
    library(here)
} |> suppressPackageStartupMessages()

gitdir <- here::here()
source(file.path(gitdir, "R/paths.R"))
source(file.path(gitdir.R, 'base_utilities.R'))
source(file.path(gitdir.R, 'base_plots.R'))
source(file.path(gitdir.R, 'dictionaries.R'))
source(file.path(gitdir.functions, 'plotting_functions.R'))


outdir <- '/home/andrea/HPC/ab1820/home/projects/2022/longvl'


# load files
dcomm <- get.dcomm(split.inland = FALSE)
dcomm[, `:=` (
    COMM_IDX = get.dcomm.idx(COMM_NUM),
    COMM_NUM_A = NULL, 
    COMM_NUM_RAW = NULL
)]
flow.1518 <- fread(path.flow.r1518)
flow.19 <- as.data.table(read_dta(path.flow.r19))

#
# combine flow across rounds
# 

cols_flow <- c('comm_num', 'locate1', 'locate2', 'resident', 'ageyrs', 'sex', 'round', 'curr_id')

flow <- lapply(
    list(flow.1518, flow.19), 
    function(DT){
        tmp <- select(DT, cols_flow)
        tmp[, curr_id := format(curr_id, scientific=FALSE)]
}) |> rbindlist() |> unique()

# check no id is repeated by round
is_all_unique <- flow[, .N , by=c('curr_id','round')][, all(N==1) ]
stopifnot(is_all_unique)

#
# Find census eligible count
#

# find  community
flow <- merge(flow, dcomm, by.x = 'comm_num', by.y = 'COMM_NUM')

# Code for ineligibility
flow[, reason_ineligible := NA_character_]

flow[locate1==10 & locate2==8, reason_ineligible := "Out_migrated"]
flow[locate1==2 & locate2==10, reason_ineligible := "Out_migrated"]
flow[locate1==13 & locate2==8, reason_ineligible := "Out_migrated"]
flow[locate1==3 & locate2==10, reason_ineligible := "Out_migrated"]
flow[locate1==5 & locate2==10, reason_ineligible := "Out_migrated"]
flow[locate1==3 & locate2==13, reason_ineligible := "Out_migrated"]
flow[locate1==6 & locate2==13, reason_ineligible := "Out_migrated"]
flow[locate1==6 & locate2==10, reason_ineligible := "Out_migrated"]

flow[locate1==7 & locate2==8, reason_ineligible := "Already_seen"]
flow[locate1==2 & locate2==7, reason_ineligible := "Already_seen"]
flow[locate1==6 & locate2==7, reason_ineligible := "Already_seen"]
flow[locate1==3 & locate2==7, reason_ineligible := "Already_seen"]
flow[locate1==5 & locate2==7, reason_ineligible := "Already_seen"]
flow[locate1==7 & locate2==7, reason_ineligible := "Already_seen"]
flow[locate1==17 & locate2==8, reason_ineligible := "Already_seen"]
flow[locate1==7 & locate2==88, reason_ineligible := "Already_seen"]
flow[locate1==7 & locate2==2, reason_ineligible := "Already_seen"]

flow[locate1==11 & locate2==8, reason_ineligible := "Dead"]
flow[locate2==11, reason_ineligible := "Dead"]

flow[resident==0, reason_ineligible := "not_resident"]

flow[ageyrs<15 | ageyrs > 49, reason_ineligible := "Not_within_eligible_age_range"]

flow[is.na(reason_ineligible), reason_ineligible := 'none']

# SET ROUND 15S IN INLAND AS 15

flow[, PARTICIPATED_TO_ROUND_RO15 := any(round == 'R015'), by= 'curr_id']
flow[round == 'R015S' & TYPE == 'inland' & PARTICIPATED_TO_ROUND_RO15 == F, round := 'R015']
flow <- flow[!(round == 'R015S' & TYPE == 'inland' & PARTICIPATED_TO_ROUND_RO15 == T)]

# find count eligible (Out migrated count as half, assuming could have migrated after?)
by_cols <- c('reason_ineligible', 'round', 'TYPE', 'ageyrs', 'sex')
cols_reasons <- unique(flow$reason_ineligible)
re <- flow[, list(count = .N), by = by_cols] |> 
    dcast.data.table(round + TYPE + ageyrs + sex ~ reason_ineligible,
        value.var = 'count') |> 
    setnafill(fill=0, cols=cols_reasons)
re[, ELIGIBLE := round(none + Out_migrated / 2)]
re <- re[ELIGIBLE != 0]

# additional variable
colnames(re) <- toupper(colnames(re))
re[ROUND %like% '^R' ,ROUND := substring(ROUND, 3) ]

# find index sex and comm
re <- re[order(ROUND, SEX, TYPE, AGEYRS)]
re[, SEX_INDEX := ifelse(SEX == 'M', 1, 0)]
re[, COMM_INDEX := ifelse(TYPE == 'fishing', 1, 0)]

# find smooth count with loess smooth
rounds <- unique(re$ROUND)
ageyrs_topredict <- re[, sort(unique(AGEYRS))]

ncen <- re[, {

    .f <- function(x)
        loess(ELIGIBLE ~ AGEYRS, span = x) |> predict(ageyrs_topredict)

    list(
        AGEYRS = ageyrs_topredict,
        ELIGIBLE = ELIGIBLE, 
        ELIGIBLE_SMOOTH.25 = .f(0.25),
        ELIGIBLE_SMOOTH.50 = .f(0.50),
        ELIGIBLE_SMOOTH.75 = .f(0.75)
    )
    
}, by=c('ROUND', 'TYPE', 'SEX')]
ncen <- merge(ncen, re, by=c('ROUND', 'TYPE', 'SEX', 'AGEYRS', 'ELIGIBLE'))

p_ncen <- plot.n.census.eligible.smooth()
filename <- "Smooth_census_eligible_count_all_round.png"
ggsave2(p_ncen, file=filename, LALA=outdir)

# choose smoothing 50
ncen[, ELIGIBLE_SMOOTH := ELIGIBLE_SMOOTH.50]
ncen <- select(ncen, -c('ELIGIBLE_SMOOTH.25', 'ELIGIBLE_SMOOTH.50', 'ELIGIBLE_SMOOTH.75'))



# save
filename <- file.path(indir.repository, 'data', 'census_eligible_individuals_table_230514.csv')
fwrite(ncen, filename , row.names = F)

# make a table for paper writing.
.replace.na <- function(x, sub){
    if(!is.numeric(x))
        x[is.na(x)] <- sub 
    return(x)
}

# get numbers
ncen_table <- subset(ncen, ! ROUND  %like% '15') |> 
    cube(j=list(N=sum(ELIGIBLE)), ,by=c('ROUND', 'TYPE', 'SEX')) |> 
    lapply(.replace.na, sub='Total') |> 
    as.data.table() |> 
    dcast.data.table(ROUND + TYPE ~ SEX, value.var='N')

# get average over rounds
.f <- function(x)
    as.integer(as.integer(x)/4)
cols <- c('M', 'F', 'Total')
ncen_table[ ROUND == 'Total', ROUND := "Average" ]
ncen_table[ ROUND == 'Average', (cols) := lapply(.SD, .f ), .SDcols = cols]

# get average over N of communities
cols <- c('M', 'F', 'Total')
new_cols <- paste(cols, 'bycomm', sep='-')
dcomm_N <- dcomm[, .(N_COMM=uniqueN(COMM_IDX)), by='TYPE']
ncen_table <- merge( ncen_table, dcomm_N, by='TYPE' )
ncen_table[, (new_cols) := lapply(.SD,  function(x) as.integer(x/N_COMM)), .SDcols=cols]

filename <- file.path(indir.repository, 'data', 'census_eligible_individuals_table_230514.rds')
fwrite(ncen, filename , row.names = F)


# NOTE: we may also want to study participation rates here
