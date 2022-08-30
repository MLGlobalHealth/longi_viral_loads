# AIMS:
# - adapt older code from Oli to multi-round, longitudinal settings
# TODO: check why we do not have any ARVMED == 2
# TODO: discuss: we are removing individuals with missing VLs. 

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
VIREMIC_VIRAL_LOAD = 1000

vl.out.dir <- file.path(out.dir, paste0('vl_', VIREMIC_VIRAL_LOAD) )
stopifnot(dir.exists(vl.out.dir))

dall <- get.dall(path.tests)

if(0) # Study ARVMED
{
        darv <- dall[HIV_STATUS == 1]
        
        cat('Assume NA ARV means no ARV usage...\n')
        tmp <- darv[ HIV_AND_VL == 1]
        tmp[is.na(ARVMED) | ARVMED != 1, ARVMED := 0]
        by_cols <- c('ROUND', 'FC', 'SEX', 'ARVMED')
        cols <- c('M', 'CL', 'CU')
        tmp[, Y := as.integer(VL_COPIES <= VIREMIC_VIRAL_LOAD)]
        tmp <- tmp[, binconf( sum(Y) , .N, return.df=T), by=by_cols]
        names(tmp) <- c(by_cols, cols)
        setkeyv(tmp, by_cols)
        supp.prop.by.arv <- copy(tmp)
        
        ggplot(supp.prop.by.arv, aes(x=FC, colour=SEX)) + 
                geom_point(aes(y=M, pch=as.factor(ARVMED)), position=position_dodge(width=.5) ) +
                geom_linerange(aes(ymin=CL, ymax=CU, linetype=as.factor(ARVMED)), position=position_dodge(width=.5) ) +
                facet_grid( ~ ROUND) + 
                scale_y_continuous(labels=scales:::percent, limits=c(0,1), expand=c(0,0)) +
                scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
                theme(legend.position='bottom') + 
                labs(x='Community type', y='Proportion of suppressed measurements',
                     linetype='Ever reported ARV', pch='Ever reported ARV',
                     title='Suppression by ARV reporting...')


        cat('Suppression among participants reporting ARV \n')
        darv[ARVMED == 1 & is.na(VL_COPIES), STUDY_ID ] -> idx
        tmp <- darv[STUDY_ID %in% idx, any(VL_COPIES==0, na.rm=T) , by='STUDY_ID']
        tmp[, cat('Out of', .N, 'HIV + participants with NA VL measurements', sum(V1), 
                  'were measured suppressed at least once\n')]

        tmp <- darv[HIV_AND_VL == 1 & ARVMED == 1]
        by_cols <- c('ROUND', 'FC', 'SEX')
        cols <- c('M', 'CL', 'CU')
        tmp[, Y := as.integer(VL_COPIES <= VIREMIC_VIRAL_LOAD)]
        tmp <- tmp[, binconf( sum(Y) , .N, return.df=T), by=by_cols]
        names(tmp) <- c(by_cols, cols)
        setkeyv(tmp, by_cols)
        supp.prop.among.report <- copy(tmp)

        # Sex plays a bigger role than community.
        ggplot(supp.prop.among.report, aes(x=FC, colour=SEX)) + 
                geom_point(aes(y=M), position=position_dodge(width=.5) ) +
                geom_linerange(aes(ymin=CL, ymax=CU), position=position_dodge(width=.5) ) +
                facet_grid( ~ ROUND) + 
                scale_y_continuous(labels=scales:::percent, limits=c(.5,1), expand=c(0,0)) +
                scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
                theme(legend.position='bottom') + 
                labs(x='Community type', y='Proportion of suppressed measurements', title='Among people reporting ever ARV...')


}


if(0)
{
        # Proportion of participants with HIV+ and na viral loads
        .f <- function(x) paste0(round(x*100, 2), '%')
        dall[, .f(mean(is.na(VL_COPIES) & HIV_STATUS==1)), by='ROUND']
}

if(0)
{
        # tmp$p for plot and tmp$DT for 'vlc' datatable
        tmp <- vl.vlprops.by.comm.gender.loc(dall, write.csv=FALSE)
        tmp$p
}


# Estimate HIV prevalence
#________________________

if(0) # already done and takes time!
{
        # Run GP to estimate prevalence by rounds.
        vl.prevalence.by.gender.loc.age.gp(dall, refit=TRUE)
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

# Among HIV positive
if(0)
{
        vl.suppofinfected.by.gender.loc.age.gp(dall, refit=FALSE)
        vl.suppofinfected.by.gender.loc.age.icar(dall)
}

# Among Entire population
if(0)
{
        # TORUN?
        vl.suppofpop.by.gender.loc.age.gp(dall, refit=TRUE)
        vl.suppofpop.by.gender.loc.age.icar(dall)
}


# GET POSTERIORS ON SUPP AMONG POP
# ________________________________

.f <- function(file)
{
        tmp <- new.env()
        round <- as.integer(gsub('^.*?round([0-9][0-9]).*?$', '\\1',file))
        load(file, envir=tmp)
        ls(tmp)
        tmp <- tmp$nspop.by.age
        tmp[, ROUND:=round]
        tmp
}
dsupp <- list.files(vl.out.dir, '220729f_suppAmongPop.*rda', full.names=T)
dsupp <- lapply(dsupp, .f)
dsupp <- rbindlist(dsupp)


# CAN WE GET INCIDENCE NOW?
# ________________________

dinc <- file.path(indir.repository, 'data', 'RCCS_1518_incidence.csv')
dinc <- fread(dinc)

dinc[, .(INCIDENCE*PY, NEWINF)]

tmp <- grep( 'ROUND|COMM|AGEYRS|SEX|INCIDEN|PREVALEN',names(dinc), value=TRUE)
dinc <- dinc[, ..tmp]

p <- ggplot(dinc, aes(x=AGEYRS, colour=ROUND)) + 
        geom_line(aes(y=INCIDENCE, group=ROUND)) + 
        facet_grid(COMM ~ SEX) + 
        viridis::scale_color_viridis() +
        theme_bw() 

# COMPARE TO AGE-SPECIFIC
dinc[, unique(AGEYRS)]
dinc
