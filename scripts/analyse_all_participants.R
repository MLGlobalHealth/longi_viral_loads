# AIMS:
# - adapt older code from Oli to multi-round, longitudinal settings

################
# DEPENDENCIES #
################

library(data.table)
library(ggplot2)
library(Hmisc)
library(rstan)
# For parallelisation across cores:
library(foreach)
library(doParallel)


################
#    PATHS     #
################

usr <- Sys.info()[['user']]
if(usr == 'andrea')
{
        indir.repository <-'~/git/longi_viral_loads'
        indir.deepsequence.data <- '~/Documents/Box/ratmann_pangea_deepsequencedata'

}else{
        indir.repository <-'~/git/longi_viral_loads'
        indir.deepsequence.data <- '~/rds/general/projects/LALALADEEPDATA'
}

out.dir <- file.path(indir.repository,'results', '220729_oli')
path.stan <- file.path(indir.repository, 'stan')
path.tests <- file.path(indir.deepsequence.data, 
                        'RCCS_R15_R20',
                        "all_participants_hivstatus_vl_220729.csv")

file.exists(
        out.dir,
        path.stan,
        path.tests
) |> all() |> stopifnot()

################
#    HELPERS   #
################

source( file.path(indir.repository,'functions/base_utilities.R') )
# source( file.path(indir.repository,'functions/preprocessing_helpers.R') )
source( file.path(indir.repository,'scripts/phyloscan.viral.load.project.R'))

# set up parallel backend
n.cores <- min(4, parallel::detectCores() - 1 )
my.cluster <- parallel::makeCluster(
        n.cores,
        type='FORK',
        outfile='.parallel_log.txt')
doParallel::registerDoParallel(cl = my.cluster)
print(my.cluster)


################
#     MAIN     #
################

VL_DETECTABLE = 400
VIREMIC_VIRAL_LOAD = 500

# dall[ HIV_STATUS == 1, mean(!is.na(HIV_VL)), by=ROUND]

# Load data: exclude round 20 as incomplete
dall <- fread(path.tests)
dall <- dall[ROUND >= 16 & ROUND <= 19]

# rename variables according to Oli's old script + remove 1 unknown sex
setnames(dall,
         c('HIV_VL', 'COMM'),
         c('VL_COPIES', 'FC') )
dall[, HIV_AND_VL := ifelse( HIV_STATUS == 1 & !is.na(VL_COPIES), 1, 0)]
dall <- dall[! SEX=='']

# tmp$p for plot and tmp$DT for 'vlc' datatable
tmp <- vl.vlprops.by.comm.gender.loc(dall, write.csv=FALSE)


# Estimate prevalences
#_____________________

if(0) # already done and takes time!
{
        # Run GP to estimate prevalence by rounds.
        vl.prevalence.by.gender.loc.age.gp(dall)
        vl.prevalence.by.gender.loc.age.icar(dall)
}


# Estimate mean viral load
# ________________________
if(0) # TORUN:
{
        vl.meanviralload.by.gender.loc.age.icar(dall)
}

# Estimate suppressed pop
# _______________________
if(0)
{
        vl.suppofinfected.by.gender.loc.age.icar(dall)
        vl.suppofinfected.by.gender.loc.age.gp(dall)
}

if(0)
{
        vl.suppofpop.by.gender.loc.age.gp(dall)
        vl.suppofpop.by.gender.loc.age.icar(dall)
}



