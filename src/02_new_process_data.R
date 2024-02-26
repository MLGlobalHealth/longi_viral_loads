library(data.table)
library(consort)
library(lubridate)
library(consort)    # For CONSORT scheme
library(grid)
read_dta <- haven::read_dta

gitdir <- here::here()
source(file.path(gitdir, "R/paths.R"))

objective_cols <- fread(path.hivstatusvl.r1520) |> names()

###########
# HELPERS #
###########

source(file.path(gitdir.functions, 'preprocessing_helpers.R'))
source(file.path(gitdir.functions, 'plotting_functions.R'))
source(file.path(gitdir.functions, "paper_statements.R"))

.summarise <- function(DT){
    dtable <- copy(DT)
    if ( ! "COMM" %in% names(dtable) ){
        dtable[, COMM := "inland"]
    }
    groupingsets(dtable, 
        j = comma(.N),
        sets = list(c("ROUND", "COMM"), c("COMM","SEX", "ROUND")),
        by = c("COMM","ROUND", "SEX")
    )  |> dcast( ROUND + COMM ~ SEX) |> 
        subset(! is.na(COMM)) |> 
        setkey(COMM, ROUND) |> 
        setnames( "NA", "Total") |>
        print()
}

########
# MAIN #
########

file.exists(
    path.nm.participants,
    path.participation,
    path.negatives.r1520,
    path.flow.r1518,
    path.flow.r19,
    path.viral.loads1,
    path.viral.loads2
) |> all() |> stopifnot()


naturemed_participants <- readRDS( path.nm.participants )
naturemed_participants[, `:=` (
    ART=NULL,
    PARTICIPATED_TO_ROUND_RO15=NULL,
    ROUND=round2numeric(ROUND),
    AGEYRS=as.integer(AGEYRS)
)]

# Let's make it as minimal as possible then: test results, viral loads.

# I need to get their sex, HIV status, age, ARVMED?

## Previous processing...
by_cols <- c("STUDY_ID", "ROUND")
allp <- process_path_participation(path=path.participation)
allp.r19 <- subset(allp, ROUND == 19)
flow <-  process_flow_datasets(path1=path.flow.r1518, path2=path.flow.r19)
flow.r19 <- subset(flow, ROUND == 19)
allp.r19 <- merge(allp.r19, flow.r19[, .SD, .SDcols=c(by_cols, "SEX", "BIRTHDAT")], by=by_cols, all.x=TRUE)
allp.r19[, `:=` (AGEYRS = floor(interval(BIRTHDAT, INT_DATE)/years(1)))]
allp.r19[, `:=` (BIRTHDAT=NULL, INT_DATE=NULL)]
dfp <- allp[, .(STUDY_ID, ROUND, FIRST_PARTICIPATION)] |> unique()

# Add R19 individuals to the "Naturemed participants" dataset 
naturemed_participants_ext <- rbind(
    naturemed_participants,
    allp.r19[, .SD, .SDcols=names(naturemed_participants)]
)
naturemed_participants_ext[, AGEYRS := as.integer(AGEYRS)]

# Get HIV tests and VL for all rounds
to_rm <- c("CURR_ID", "REGION", "SEX", "HIVDATE", "INT_DATE", "COMM_NUM")
dneg <- process_hiv_negatives(path=path.negatives.r1520)
dvl  <- process_hiv_vls_for_hivpositives(path1=path.viral.loads1, path2=path.viral.loads2) 
dneg[, (to_rm) := NULL]
dvl[, (to_rm) := NULL] 
tests <- rbind(dneg, dvl, fill=TRUE)

# merge to participants
check <- merge(naturemed_participants_ext, tests, by=c("STUDY_ID", "ROUND"), all.x=TRUE)
check[, table(COMM.y == COMM.x, useNA = "ifany")]
check[, table(AGEYRS.y == AGEYRS.x, useNA = "ifany")]
check[ abs(AGEYRS.y - AGEYRS.x)> 3, STUDY_ID] -> ids
str(check)
check[j=`:=` (
    AGEYRS = fcoalesce(AGEYRS.x, AGEYRS.y), 
    COMM = fcoalesce(COMM.x, COMM.y),
    AGEYRS.x=NULL, AGEYRS.y=NULL,
    COMM.x=NULL, COMM.y=NULL
)]
check[ is.na(SEX), SEX := "U"]

# Merge with first participation
check <- merge(check, dfp, by=c("STUDY_ID", "ROUND"), all.x=TRUE)
setkey(check, STUDY_ID, ROUND)
check[ STUDY_ID == "K122497", FIRST_PARTICIPATION := c(TRUE, FALSE)]
check[, {
    fp_ids <- which(FIRST_PARTICIPATION == TRUE)
    l <- length(fp_ids) 
    if( l > 1){ FALSE } else if (l == 0){ TRUE}else{
        fp_ids == 1
    }
}, STUDY_ID][V1 == FALSE, stopifnot(.N==0)]

.summarise(subset(check, ROUND %between% c(16, 19) ))
.summarise(subset(naturemed_participants, ROUND %between% c(16, 19) ))

cat("
##                            
##  Make supplementary tables 
##                            
")

by_cols <- c("COMM", "SEX", "ROUND")
# Number of census-eligible is calculated in get_census_eligible_count.R
# accounting for migration, death, and other factors
ncen <- fread(path.census.eligible)
ncen[, `:=` (ROUND = round2numeric(ROUND))]
setnames(ncen, "TYPE", "COMM")
ncen_tot <- ncen[, .(N_EL=sum(ELIGIBLE)), by=by_cols]

check_r1619 <- subset(check, ROUND %between% c(16, 19) & SEX != "U")
tots_r1619 <- check_r1619[, .(N_PAR=uniqueN(STUDY_ID)), by=c("COMM","SEX", "ROUND")]

dtable <- merge(ncen_tot, tots_r1619, by=c("COMM", "SEX", "ROUND"))
dtable[, `:=` (PERC= round(100*N_PAR/N_EL, 1))]
make.table.eligible.participants(dtable, cols=c("N_EL", "N_PAR"), splitrow=TRUE) 

# check_r1619[, table(is.na(HIV_STATUS)), by=by_cols]
# check_r1619[is.na(HIV_STATUS), hist(AGEYRS)]
make.table.first.time.participants(check_r1619)

# Add dates and then finish on Inkscape
options(txt_gp = gpar(cex = 1, fontsize=6, lw=.5, col = "black", lex=1))
g <- build_grid(make_consort_diagram(check_r1619))
cmd <- ggsave2("consort_diagram.pdf", 
    p = g, 
    LALA=file.path(OUTDIR, "figures"), 
    w=21, h=15, u="cm")
system(cmd)

# Save
filename <- file.path(indir.deepdata.r1520,'all_participants_hivstatus_vl_240214.csv')
if( ! file.exists(filename) )
{
    sprintf('Saving file %s...', filename) |> catn()
    fwrite(check_r1619, filename)
}else{
    sprintf('File %s already exists...', filename) |> catn()
}
