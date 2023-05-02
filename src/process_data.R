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

###########
# HELPERS #
###########

source(file.path(gitdir.R, 'base_utilities.R'))
source(file.path(gitdir.functions, 'preprocessing_helpers.R'))
source(file.path(gitdir.functions, 'plotting_functions.R'))

.dict.class <- function(x)
{
    lbls <- c("durably_suppressed", "durably_viremic", "newly_suppressed", "intermittently_viremic", "newly_viremic")
    nms <- gsub('_',' ',lbls)
    names(nms) <- lbls
    unname(nms[x])
}

########
# MAIN #
########

# Set option
THRESHOLD <- 50
fix.incoherent.arv.reporting <- 1
fix.death.dates <- 1
save.images <- FALSE

# explore dataset
dvl <- read.csv(path.viral.loads1, fill=TRUE, comment.char="") |> 
    as.data.table()
setkey(dvl, study_id, hivdate)

cat('Excluding R20:')
make_vl_samplesize_table(dvl, excludeR20 = TRUE)
cat('Including R20:')
make_vl_samplesize_table(dvl, excludeR20 = FALSE)

# Empty Entries to NA
cols <- dvl[, sapply(.SD, is.character)]
cols <- names(cols)[cols]
dvl[, (cols):=lapply(.SD, .empty2na), .SDcols=cols]

# Remove uninformative columns (ie: only one entry)
idx <- dvl[, which(lapply(.SD, uniqueN) == 1)]
cols <- names(idx)
cat('Removing constant columns, as uninformative:\n')
dvl[, lapply(.SD,unique) , .SDcols=cols]
dvl[, (cols):=NULL]

# Translate round from string to numeric
dvl[, round:=round2numeric(round)]

# Fix unrealistic visit dates
fix_visit_dates_before_1990(dvl)

# Translate community keys to type
dvl <- community.keys2type(DT=dvl)

# repeated_names <- gsub('R0[0-9][0-9].*$','', names(dvl))
# repeated_names <- table(repeated_names) > 1

# 'unique columns' are: age, community, (CD4HTSS15), community number, current id, death data and cause, 1st post/last neg date and visit, hiv and hiv_vl and date, 
# 'repeated columns concern' arvmed, current arvmed, cd4 count and date, interview date
# some columns are repeated for Rounds 11-20 (consider also round R015S)!

# ocmpare ageyrs vs conf_age
tmp1 <- unique(dvl[,.(study_id, ageyrs,conf_age)]) 
tmp1[ageyrs != conf_age, 
     cat('There are', .N, 'individuals with ageyrs != conf_age. (REMOVING conf_age). \n')]
dvl[, conf_age := NULL]

# study curr_id
tmp <- unique(dvl[, .(study_id, round, curr_id)])
idx <- tmp[,uniqueN(curr_id), by=c('study_id', 'round')][V1 > 1, {cat( uniqueN(study_id), '(study_id,round)s have multiple curr_id\n'); list(study_id, round)}]
# tmp[study_id %in% idx]

# Store separate info in dsero, darv, dcd4 and dintdates, dicd
dvl <- make_relational_database(dvl)
tmp <- sapply(c('dsero','darv','dcd4','dintdates','ddeath', 'dlocate', 'dbirth'), exists)
stopifnot(all(tmp))

# load Allpcr dataset and substitute missing hiv_vl 
tmp <- fill_na_vls_with_allpcr_data(file=path.viral.loads2)
cat('The new dataset provided viral load measurements in rounds:\n')
tmp$tbl
dvl <- tmp$dvl

study_low_level_viremia(dvl)

# extract and plot 
cat('How many participants had measurements in both rounds 15 and 15S?\n')
dvl[round %in% c(15, 15.5), sum(!is.na(hiv_vl)) > 1, by='study_id'][, table(V1)]

if(save.images)
{
    filename <- file.path(out.dir, 'nm_vls_per_round.png')
    plot_rounds_vl_collection(DT=dvl, rnd=15, filename=filename)
    plot_rounds_vl_collection(DT=dvl, rnd=15, 
                              community='fishing', atleast=4)
    filename <- file.path(out.dir, 'nm_3vls_per_round.png')
    plot_rounds_vl_collection(DT=dvl, rnd=15, atleast=3, filename=filename)
    filename <- file.path(out.dir, 'nm_4vls_per_round.png')
    plot_rounds_vl_collection(DT=dvl, rnd=15, atleast=4, filename=filename)
    plot_rounds_vl_collection(DT=dvl, rnd=15, atleast=5)
}

if(save.images) # STUDY NA VL MEASUREMENTS
{
    # Loooooads of NA VL measurements!
    # most of them are in the first 3 rounds! 
    dvl[, mean(is.na(hiv_vl))]
    
    # 3 plots to do, assuming:
    # - NA correspond to an actual measurement whose result is unknown
    # - NA correspond to a non-measurement 
    # - one with NA alone


    cols <- c('study_id', 'round', 'hiv_vl')
    tmp <- dvl[, ..cols]
    tmp[, round := round2factor(round)]


    if(0) # New request
    {
        tmp1 <- tmp[, DUMMY:=as.integer(gsub('[A-Z]', '', round))]
        tmp1 <- tmp1[, .(study_id, DUMMY)]

        tmp1 <- dcast(tmp1, study_id ~ DUMMY)
        colnames(tmp1)[-1] <- paste0('R',colnames(tmp1)[-1])

        .f <- function(x) round(mean(x != 0)*100, 2 )
        rbind(
                tmp1[R15 > 0, lapply(.SD, .f)],
                tmp1[R16 > 0, lapply(.SD, .f)],
                tmp1[R17 > 0, lapply(.SD, .f)],
                tmp1[R18 > 0, lapply(.SD, .f)],
                tmp1[R19 > 0, lapply(.SD, .f)],
                tmp1[R20 > 0, lapply(.SD, .f)]
        ) -> tmp1
        setnames(tmp1, 'study_id', 'round')
        tmp1[, round := paste0('R', 15:20)]
        knitr::kable(tmp1)
        # For all study_ids with DUMMY == 15, how many have 
    }

    tmp[, count_with_na := 1]
    tmp[, count_without_na := dplyr::if_else(is.na(hiv_vl), 0 , 1)]

    g <- ggplot(data=tmp, aes(x=round, fill=!is.na(hiv_vl))) + 
            geom_histogram(stat='count', position='stack') + 
            theme_bw() + 
            theme(legend.position='bottom') + 
            labs(title='Missing Viral Load measurements',
                 subtitle='Most of the reported VL measurement in rounds 15, 16 are NA',
                 y='Number of VL measurements',
                 x='Rounds',
                 fill='Non missing VL measurement')
    filename=file.path(out.dir, 'na_vlmeasurements_byround.png')
    ggsave(filename, g, w=8, h=6)

    # By number of successive visits
    setkey(tmp, study_id, round)
    .f <- function(x){ sum(x) - cumsum(x) }
    tmp[, N_successive_visits := .f(count_with_na) , by=study_id]
    tmp[, N_successive_visits := as.factor(N_successive_visits)]
    tmp[, first_measurement :=  cumsum(count_with_na) == 1,by='study_id']
    yrange <- tmp[round=='R15', sum(count_with_na)]

    g <- ggplot(data=tmp, aes(x=first_measurement, fill=N_successive_visits )) + 
        geom_histogram(stat='count', position='stack') + 
        theme_bw() + 
        facet_grid(~round) + 
        viridis::scale_fill_viridis(discrete=TRUE) + 
        ylim(0, yrange) + 
        theme(panel.margin = grid::unit(0, "lines"),
              legend.position='bottom',
              legend.box.just='center',
              legend.title.align=0.5,
              legend.box='horizontal') + 
        guides(fill=guide_legend(nrow=1, byrow=TRUE,
                                 title.position='left',
                                 label.position='bottom',
                                 label.hjust=0.5,
                                 title.hjust=0.5)) + 
        labs(title='Measurements follow up',
             subtitle='New participants are more likely to be lost to follow up.',
             fill='Successively collected VL measurements per participant',
             x='First measurement in R15-R20?',
             y='Number of VL measurements') 
    g

    filename=file.path(out.dir, 'vlmeasurements_followup.png')
    ggsave(filename, g, w=8, h=6)

    # what if we actually use the rounds of diagnosis instead
    all(tmp[, unique(study_id)] %in% dsero[, unique(study_id)])
    cols <- c('study_id', 'firstposv')
    tmp1 <- dsero[, ..cols]
    if(!is.numeric(tmp1$firstposv)){
        tmp1[, firstposv := round2numeric(firstposv)]
        tmp1[, firstposv := round2factor(firstposv)]
        tmp1[ grepl('NA', firstposv), firstposv:=NA]
    }
    tmp <- merge(tmp, tmp1, by='study_id')
    cat('The proportion of participants with unknown first positive visit is:\n')
    cat(tmp[, is.na(firstposv)[1], by='study_id'][, round(mean(V1),2)], '\n')
    idx <- tmp[as.character(round)==as.character(firstposv) & first_measurement == FALSE, study_id]
    tmp[study_id %in% idx]
    tmp[as.character(round)!=as.character(firstposv), first_measurement := FALSE]

    tmp[, is.na(firstposv)[1], by='study_id'][, mean(V1)]

    g <- ggplot(data=tmp, aes(x=first_measurement, fill=N_successive_visits )) + 
        geom_histogram(stat='count', position='stack') + 
        theme_bw() + 
        facet_grid(~round) + 
        viridis::scale_fill_viridis(discrete=TRUE) + 
        ylim(0, yrange) + 
        theme(panel.margin = grid::unit(0, "lines"),
              legend.position='bottom',
              legend.box.just='center',
              legend.title.align=0.5,
              legend.box='horizontal') + 
        guides(fill=guide_legend(nrow=1, byrow=TRUE,
                                 title.position='left',
                                 label.position='bottom',
                                 label.hjust=0.5,
                                 title.hjust=0.5)) + 
        labs(title='Measurements follow up (adjusted by `firstposvd`)',
             subtitle='New participants are more likely to be lost to follow up.',
             fill='Successively collected VL measurements per participant',
             x='First measurement in R15-R20?',
             y='Number of VL measurements') 
    g
    
    filename=file.path(out.dir, 'vlmeasurements_followup_adjusted.png')
    ggsave(filename, g, w=8, h=6)
}

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
darv[ cuarvmed == 1 & is.na(arvmed), cat(' -',.N, 'participants reported current, but not ever ARVMED ')]

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

# remove conf_age
dneg[ ageyrs != conf_age]
dneg[ , conf_age := NULL]
# Fix date, round and get community type
dneg[, int_date := as.Date(int_date, '%d/%m/%Y')]
dneg[, hivdate := as.Date(hivdate, '%d%b%Y')]
dneg[, round := round2numeric(round)]
dneg <- community.keys2type(dneg)

# tables
by_cols <- c('round',  'comm', 'sex')
tbl <- dneg[, .N, by=by_cols]
setkeyv(tbl, by_cols) 
# tbl

# remove columns not in dvl and merge
dneg[, `:=` (hiv_vl=0, HIV_STATUS=0, int_date=NULL)]
dvl[,  `:=` (HIV_STATUS=1, int_date=NULL)]


# stopifnot(names(dvl) %in% names(dneg))
stopifnot(names(dneg) %in% names(dvl))
idx <- names(dneg) %in% names(dvl)
cols <- names(dneg)[idx]
dvl[, hivdate := as.Date(hivdate)]

dall <- rbind(dneg, dvl[, ..cols]) |> suppressWarnings()

dall <- merge(dall, darv[, .(study_id, round, arvmed)], 
      by=c('study_id', 'round'), all.x=TRUE)
# study arvmed
# tmp <- dall[HIV_STATUS == 1, mean( !is.na(arvmed)) , by=c('round', 'sex')]
# knitr::kable(dcast(tmp,round ~ sex))

setnames(dall, names(dall), toupper(names(dall)) )
cols <- c('STUDY_ID', 'ROUND', 'SEX', 'AGEYRS' , 'COMM', 'COMM_NUM',  'CURR_ID', 'HIV_STATUS', 'HIV_VL', 'HIVDATE', 'ARVMED')
setcolorder(dall, cols)


# write to file
filename <- file.path(indir.deepdata.r1520,'all_participants_hivstatus_vl_230502.csv')
if( ! file.exists(filename) )
{
    sprintf('Saving file %s...', filename) |> catn()
    fwrite(dall, filename)
}else{
    sprintf('File %s already exists...', filename) |> catn()
}

files <- list.files( indir.deepdata.r1520, pattern="viral", full.names =TRUE)
lapply(files, fread) |> 
    lapply(function(DT) DT[, .(uniqueN(study_id)), by='round'])

if(0)
{
    if(0)
    {
        # participants moving from inland to fishing communities.
        idx <- dvl[, uniqueN(comm), by='study_id'][V1 > 1, unique(study_id)] 
        dvl[study_id %in% idx]
    }

    # define trajectory types
    #________________________

    date_utt <- 2013 # according to grabowsky2021

    # Maybe this not needed as reported above
    # Round 15 only has VLs for fishing communities
    dvl_15 <- count_vls_by_round(DT=dvl, rnd=15, subset_n=4)
    dvl_16 <- count_vls_by_round(DT=dvl, rnd=16, subset_n=4)


    # Define trajectories
    cat('Trajectories for round 15 and larger')
    dclass_15 <- define_trajectories(dvl_15)
    plot_classification(DT=dvl_15, dclass_15)

    cat('Trajectories for round 16 and larger')
    dclass_16 <- define_trajectories(dvl_16)
    plot_classification(DT=dvl_16, dclass_16)

    # process outputs of relational database
    # __

    # DCD4 <- copy(dcd4)
    dcd4 <- process_dcd4(dcd4)

    stopifnot(dcd4[!is.na(cd4), all(!is.na(cd4date))] & dcd4[!is.na(cd4date), all(!is.na(cd4))])
    tmp <- dcd4[, any(!is.na(cd4)) ,by='study_id']
    cat(sum(tmp$V1), 'out of ', tmp[, .N], 'study_ids (', round(tmp[, 100*sum(V1)/.N], 2)  ,'%)have at least a CD4 count measurement.\n')

    tmp <- dcd4[, .(N=sum(!is.na(cd4))) ,by='study_id']
    cat('Number of study_ids with with N CD4 counts.\n')
    knitr::kable(tmp[, table(N)])

    # sex:

    # Classifications stratified by age, sex, community
    if(save.images)
    {
        # R15 :
        plot_classes_by_sex_age(dclass_15, dvl_15,
                                file.path(out.dir, 'classes_by_agesex_15.png'))
        plot_classes_by_comm_age(dclass_15, dvl_15,
                                 file.path(out.dir, 'classes_by_agecomm_15.png'))
        # R16 :
        plot_classes_by_comm_age(dclass_16, dvl_16,
                                 file.path(out.dir, 'classes_by_agecomm_16.png'))
        plot_classes_by_sex_age(dclass_16, dvl_16,
                                file.path(out.dir, 'classes_by_agesex_16.png'))
    }

    # get arv

    darv <- process_darv(darv)
    cat('The proportion of participant who did not report ARV on their first visit is:\n',
        darv[, firstarv[1], by=study_id][, mean(is.na(V1))], '\n')

    dvl_tmp <- merge(darv[, .(study_id, round, arvmed, cuarvmed)], dvl, by=c('study_id', 'round'), all.y=T)

    # joint cd4 and vl
    if(save.images)
    {
        # higher CD4 counts in F vs M consistent with 
        # https://www.researchgate.net/publication/253336206_Population-Based_CD4_Counts_in_a_Rural_Area_in_South_Africa_with_High_HIV_Prevalence_and_High_Antiretroviral_Treatment_Coverage/figures?lo=1

        plot_scatter_cd4vl_bygroup(dcd4, dvl, group='sex', 
                                   filename=file.path(out.dir, 'scatter_cd4vl_by_sex.png'))
        plot_scatter_cd4vl_bygroup(dcd4, dvl, group='comm', 
                                   filename=file.path(out.dir, 'scatter_cd4vl_by_comm.png'))
        plot_scatter_cd4vl_bygroup(dcd4, dvl_tmp, group='arvmed', 
                                   filename=file.path(out.dir, 'scatter_cd4vl_by_arvmed.png'))
    }
    rm(dvl_tmp)

    # get deaths
    ddeath <- process_deaths(ddeath)


    # compare demographics
    cols <- c('study_id', 'region', 'comm', 'sex')
    tmp <- dclass_15[, .(study_id,class)]

    # Participants with defined first arv use are underepresented in the durably viremic class.
    # Makes sense
    # However NOTE THAT THE MERGING MAY NOT BE COMPLETE!!!!
    tmp1 <- unique(darv[,.(study_id, firstarv)])
    tmp <- merge(tmp, tmp1, by='study_id')
    tmp[, list(dvl,
               N_firstarv=sum(!is.na(firstarv)),
               N_tot=.N,
               mean=round(mean(!is.na(firstarv)),3)),
      by=class]


    # save image
    rm(tmp, tmp1, tmp2)
    cat('process_data.R completed, saving environment...\n')
    filename <- file.path(out.dir, 'preprocessed_data.RData')
    save.image(file=filename)
    cat('done\n')


    tmp <- dclass_16[, .(study_id, class = as.factor(class))]
    tmp1 <- dvl_16[, lapply(.SD, function(x) x[length(x)]), by='study_id']
    tmp1 <- tmp1[, .(study_id, comm, sex, ageyrs)]
    tmp <- merge(tmp, tmp1)

    # to specify baseline:
    tmp[, class := relevel(class, ref='durably_suppressed')]
    test <- nnet::multinom(data = tmp, class ~ sex + comm)

    summary(test)

    # compare predictions
    cbind(
          tmp[, .(study_id, class)], 
          fitted(test)
          ) -> tmp1
}
