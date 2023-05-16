# cuarvmed = current arvmed
# everarvmed = ever used arv

require(data.table)
require(ggplot2)
require(readxl)

# TODO:
# There seem to be a couple of individuals whose firstpositive visit is successive to
# a few 0 hiv_vl measurements. So:
# - remove those measurements
# - study other individuals for which first HIV_VL measuremnts are 0 in dvl

# TODO: load the new dataset  from Kate and compare it with what we have


######### 
# Paths # 
#########

gitdir <- here::here()
source(file.path(gitdir, 'R/paths.R'))

outdir <- '/home/andrea/HPC/ab1820/home/projects/2022/longvl'
outdir.figures <- file.path(outdir, 'figures')
outdir.tables <- file.path(outdir, 'tables')

###########
# HELPERS #
###########

source(file.path(gitdir.functions, 'preprocessing_helpers.R'))
source(file.path(gitdir.functions, 'plotting_functions.R'))

########
# MAIN #
########

FINAL_COLS <- c('STUDY_ID', 'ROUND', 'SEX', 'AGEYRS' , 'COMM', 'COMM_NUM',  'CURR_ID', 'HIV_STATUS', 'HIV_VL', 'HIVDATE', 'ARVMED')
FINAL_COLS_LOWER <- tolower(FINAL_COLS)

# Set option
THRESHOLD <- 50
fix.incoherent.arv.reporting <- 1
fix.death.dates <- 1
save.images <- FALSE

# explore dataset
dvl <- read.csv(path.viral.loads1, fill=TRUE, comment.char="") |> 
    as.data.table() |> 
    subset(nchar(study_id) < 8) |> 
    empty2naDT() |> 
    remove.columns.with.unique.entry()

dvl[, round:=round2numeric(round)]
setkey(dvl, study_id, hivdate)
fix_visit_dates_before_1990(dvl)
dvl <- community.keys2type(DT=dvl, by.DT='comm_num')

make_vl_samplesize_table(dvl, rounds = 16:19)

# conf age always == ageyrs, except in 2 cases
dvl[, conf_age := NULL]

# Store separate info in dsero, darv, dcd4 and dintdates, dicd
dvl <- make_relational_database(dvl)
sapply(c('dsero','darv','dcd4','dintdates','ddeath', 'dlocate', 'dbirth'), exists) |> all() |> stopifnot()

# load Allpcr dataset and substitute missing hiv_vl 
tmp <- fill_na_vls_with_allpcr_data(file=path.viral.loads2)
dvl <- tmp$dvl

dvl |> 
    subset(select=intersect(names(dvl), FINAL_COLS_LOWER) ) |> 
    subset(round %between% c(16, 19)) |>
    dplyr::pull(hiv_vl) |> is.na() |>
    table() |> proportions()

study_low_level_viremia(dvl)

# extract and plot 
dvl[round %in% 16:19, 
    as.data.table(table(HAS_VL=!is.na(hiv_vl))),
    by='round']

plot_rounds_vl_collection(DT=dvl, atleast=1)


# save dvl in R15_R20 folder:
filename <- file.path(indir.deepdata.r1520, "viral_loads_r15r20_processed_230502.csv")

if( ! file.exists(filename) )
{
    sprintf('Saving file %s ...', filename) |> cat()
    fwrite(dvl, filename)
}else{
    sprintf('File %s already exists...', filename) |> catn()
}


# Process darv
#_____________

melt(darv, id.vars="study_id",) -> tmp

tmp[, round := round2numeric( gsub('^.*med', '', variable) )]
tmp[, variable := gsub('^(.*med).*$', '\\1', variable) ]

dcast(tmp, 
      study_id + round ~ variable
) -> darv

cat('Check inconsisten ARV reporting, eg: cuarvmed but not arvmed\n')
darv[ cuarvmed == 1 & is.na(arvmed),
    sprintf(' -%s participants reported current, but not ever ARVMED', .N) |> cat()]

cat('If a participant reported cuarvmed or arvmed at a given round, set arvmed=1 for succesive rounds\n')
darv[ , arvmed := {
        idx <- ! is.na(arvmed) | ! is.na(cuarvmed)
        idx <- which(idx)[1]
        if ( ! is.na(idx) )
                arvmed[idx:length(arvmed)] <- 1L
        arvmed
} , by='study_id']

#########################
# Negative Participants #
#########################


dneg <- fread(path.negatives.r1520)

# Fix date, round and get community type
dneg[, `:=` (
    conf_age = NULL,  # again, almost equal to ageyrs
    int_date = as.Date(int_date, '%d/%m/%Y'),
    hivdate = as.Date(hivdate, '%d%b%Y'),
    round = round2numeric(round))] 
dneg <- community.keys2type(dneg)

# tables
by_cols <- c('round',  'comm', 'sex')
tbl <- dneg[, .N, by=by_cols]
setkeyv(tbl, by_cols) 
# tbl

# remove columns not in dvl and merge
dneg[, `:=` (hiv_vl=0, hiv=NULL, HIV_STATUS=0, int_date=NULL, hh_num=NULL, resident=NULL, mobility=NULL)]
dvl[,  `:=` (HIV_STATUS=1, int_date=NULL)]


# stopifnot(names(dvl) %in% names(dneg))
stopifnot(names(dneg) %in% names(dvl))
idx <- names(dneg) %in% names(dvl)
cols <- names(dneg)[idx]
dvl[, hivdate := as.Date(hivdate)]

dall <- rbind(dneg, dvl[, ..cols]) |> suppressWarnings()
dall <- merge(dall, darv[, .(study_id, round, arvmed)], 
      by=c('study_id', 'round'), all.x=TRUE)

setnames(dall, names(dall), toupper(names(dall)) )
setcolorder(dall, FINAL_COLS)

##############################
# Add first participant status
##############################


firstparticipants <- haven::read_dta(path.participation) |> 
        as.data.table() |> 
        subset(select=c('study_id', 'round', 'qst_1strnd'))
names(firstparticipants) <- toupper(names(firstparticipants)) 

firstparticipants <- firstparticipants[, .(
    STUDY_ID, 
    ROUND=round2numeric(ROUND),
    FIRST_PARTICIPATION = as.integer(ROUND == QST_1STRND))]
firstparticipants <- subset(firstparticipants, ROUND %between% c(16, 19))

dall <- merge(dall, firstparticipants,
    all.x=TRUE, all.y=TRUE,
    by=c('STUDY_ID', 'ROUND'))


#####################
# Prettify and save #
#####################

# write to file
filename <- file.path(indir.deepdata.r1520,'all_participants_hivstatus_vl_230515.csv')
if( ! file.exists(filename) )
{
    sprintf('Saving file %s...', filename) |> catn()
    fwrite(dall, filename)
}else{
    sprintf('File %s already exists...', filename) |> catn()
}


if(0)   # double check participants
{
    missing2 <- subset(dall,
        ROUND %between% c(16, 19) &  (is.na(HIV_STATUS) | is.na(FIRST_PARTICIPATION)),
        select=c('STUDY_ID', 'ROUND', 'HIV_STATUS', 'FIRST_PARTICIPATION'))
    missing2[, `:=` (
        MISSING = fifelse(is.na(HIV_STATUS), yes="HIV status", no="First questionaire round"),
        HIV_STATUS = NULL, FIRST_PARTICIPATION = NULL
        )]
    filename <- file.path(outdir.tables, 'studyidround_missingfrom_firstparticipantsORvlandnegativedata.rds')
    saveRDS(missing2, filename)
}
