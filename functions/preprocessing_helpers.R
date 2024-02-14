table_f <- function(x, print=FALSE)
{
    # TODO: for some reason it does not print the result
    # when print=FALSE
    tmp <- as.data.table(table(x))
    tmp[, `F(%)`:= round(100*N/sum(N), 2)]
    print(knitr::kable(tmp))
    return(tmp)
}

participant_comm2mode <- function(DT, quietly=FALSE)
{
    if(!quietly)
    {
            cat(
            'Some participants travelled in fishing and inland communities.\n',
            'We only want to attribute one classification to each of these. \n',
            'So, set community to mode of where measurements were taken inland/fishing,\n',
            'with preference for fishing\n')
    }

    getmode <- function(x){
            x <- table(x)
            names(which.max(x))
    }
    tmp <- DT[,.(mode=getmode(comm)), by='study_id']

    # substitute by data.table keys
    setkey(tmp, study_id)
    setkey(DT,study_id)
    DT[tmp, comm:=mode]

    # return
    stopifnot(DT[, uniqueN(comm), by='study_id'][V1 > 1, .N == 0])
    DT
}

count2suppstatus <- function(x, character=FALSE)
{
        y <- (x < THRESHOLD)
        if(character) return(dplyr::if_else(y,'suppressed','viremic'))
        y
}

community.keys2type <- function(DT, by.DT='comm_num')
{
    # read in community keys
    # DT <- copy(dvl); by.DT <- 'comm_num'
    if('comm' %in% names(DT)) 
        stop("Comm already exists")

    dkeys <- fread(path.community.types)
    names(dkeys) <- tolower(colnames(dkeys))
    dkeys[, comm := ifelse(grepl('^f',comm_num_a), 'fishing', 'inland')]
    dkeys[, comm_num_a := NULL]

    # merge
    idx <- which(names(DT) == by.DT)
    tmp <- merge(DT, dkeys, all.x=TRUE, 
        by.x=by.DT , by.y='comm_num_raw', sort=FALSE)
    if(tmp[, any(is.na(comm))]) warning('At least one comm_num not found')

    # reorder as original, with comm preceding comm_num
    rnk <- names(DT) 
    rnk <- c(rnk[1:(idx-1)], 'comm', rnk[idx:length(rnk)])
    tmp <- tmp[, .SD, .SDcols=rnk]
    tmp
}

make_vl_samplesize_table <- function(DT, rounds)
{
    cols <- c('study_id', 'sex', 'round','hivdate', 'hiv_vl')
    tmp <- DT[, .SD, .SDcols=cols]
    tmp <- tmp[!is.na(hiv_vl),]
    tmp <- tmp[round %in% rounds]

    tmp1 <- tmp[, .(VL=length(hiv_vl)) ,by=c('study_id', 'sex')]
    tmp1 <- tmp1[, .N, by=c('VL', 'sex')]
    tmp1 <- dcast(tmp1, VL~sex)
    tmp1[, Total:=M+F]
    tmp1[, Total2:=Total*VL]
    knitr::kable(tmp1, caption='Number of vl measurements per person')
}

fill_na_vls_with_allpcr_data <- function(file=file.viral.loads2, DT=dvl)
{ 
    cat('===\n Studying ', file, '...\n===\n\n')

    # Warnings about 6 dates prior to 1900, but shouldn t be an issue
    dvl2 <- suppressWarnings(readxl::read_xlsx(file))  |> 
        as.data.table()
    setkey(dvl2, study_id, round)
    setcolorder(dvl2, c('study_id', 'round'))

    dvl2[, round := gsub('PLA', '', round)] 
    dvl2[, `:=` (studyid = NULL, round = round2numeric(round))]

    if(0)
    {
        cols <- c('study_id', 'round')
        idx <- DT[!is.na(round), .(study_id, round, data='first')]
        idx2 <- DT[!is.na(round), .(study_id, round, data='second')]
        tmp <- merge(idx, idx2, all.x=TRUE, all.y=TRUE, by=cols)
        tmp[is.na(data.y)]
        
        merge(
            DT[, .(study_id, round, hiv_vl)],
            dvl2[, .(study_id, round, COPIES = fcoalesce(new_copies, copies))],
            by=cols
        ) -> tmp

        tmp[!is.na(hiv_vl) & hiv_vl == 0, table(COPIES)]
        tmp[!is.na(hiv_vl) & COPIES == 'BD', table(hiv_vl)]
        tmp[!is.na(hiv_vl), table(hiv_vl < 500 & hiv_vl > 0)]
        
        diff <- tmp[!is.na(hiv_vl) & hiv_vl != copies]
        diff[, COPIES := fcoalesce(new_copies, copies)]
        diff[hiv_vl == 0 & copies != 'BD' & !is.na(as.integer(copies)), ]
    }
    
    # What to do with non-numeric entries?]
    dvl2[, guessed := FALSE]
    dvl2[copies %like% '[A-Z]|<|>', table(copies)] |> knitr::kable()

    if(0)
    {
        dvl2[copies %like% '[A-Z]|<|>', table(new_copies)]
        dvl2[copies %like% '<|>', table(new_copies)]
        dvl2[copies %like% 'BD', table(new_copies)]

        dvl3 <- copy(dvl2)
        dvl3[copies %like% 'BD' & !(new_copies %like% 'BD' | is.na(new_copies)), 
             copies2 := new_copies]

        cols <- c('study_id', 'round')
        merge(
            DT[, .SD, .SDcols=c(cols, 'hiv_vl')],
            dvl2[, .SD, .SDcols=c(cols, 'copies', 'new_copies')],
            by=cols
        ) -> dcomp

        dcomp[ copies == 'BD',  table(hiv_vl, new_copies)] |> knitr::kable()


        dcomp[copies=='BD' & !is.na(copies2),]
        dvl2[copies == 'BD']

    }

    dvl2[copies %like% '[A-Z]|<|>' & ! new_copies %like% '[A-Z]|<|>' & !is.na(new_copies)]
    dvl2[copies %like% '[A-Z]|<|>', guessed := TRUE]

    cat(' - For "<" ranges, keep the upper bound...\n')
    cols <- c('copies', 'new_copies')
    dvl2[, (cols) := lapply(.SD, function(x) gsub(',', '', x) ) , .SDcols=cols]
    dvl2[, (cols) := lapply(.SD, function(x) gsub('<|< ', '', x) ) , .SDcols=cols]

    # It would seem that new-copies should be kept if it differs from copies!
    # BD is generally equivalent with 0 hiv_vl (~70% of times, else NA)
    copiestmp <- merge(DT[, .(study_id, round, hiv_vl, hivdate)], dvl2, by=c('study_id', 'round'))
    # tmp[hiv_vl != new_copies & hiv_vl != 0]
    # tmp[new_copies %like% "[A-Z]", .(hiv_vl, new_copies)][, mean(!is.na(hiv_vl))]

    # 
    cat(' - setting BD, ND to 0...\n')
    cols <- c('copies', 'new_copies')
    dvl2[, (cols) := lapply(.SD, function(x) gsub('^BD|^ND|^Not Dete.*?$', '0',x)), .SDcols=cols]
    dvl2[, (cols) := lapply(.SD, function(x) gsub('[A-Z]', '0',x)), .SDcols=cols]
    dvl2[, hiv_vl := copies]
    dvl2[copies != new_copies, hiv_vl := new_copies]
    dvl2[, hiv_vl := as.numeric(hiv_vl)]

    # merge?
    tmp <- DT[is.na(hiv_vl), .(study_id, round)]
    tmp <- merge(tmp, dvl2, all.x=TRUE)
    tmp[, cat("Out of ", .N, " NA viral loads in the first dataset, ", round(100*mean(!is.na(hiv_vl)),2), '% are reported in the second one', '\n' )]

    cat('Substituting....\n')
    setnames(tmp, 'hiv_vl', 'hiv_vl2')
    tmp[!is.na(hiv_vl2), round]
    capt <- "The new dataset provided viral load measurements in rounds"
    tbl <- knitr::kable(tmp[!is.na(hiv_vl2), table(round)], caption=capt)

    DT <- merge(DT,
          tmp[!is.na(hiv_vl2), .(study_id, round, hiv_vl2)],
          all.x=TRUE, by=c('study_id', 'round'))
    DT[!is.na(hiv_vl2), hiv_vl := hiv_vl2]
    DT[, hiv_vl2:=NULL]

    return(list(dvl=DT, tbl=tbl))
}

make_relational_database <- function(DT)
{
        # cols contains the names of the columns in DT that we want to 
        # remove and store in other data.tables. Additional columns of
        # interest which we do not want to delete are additionaly 
        # specified in the .SDcols argument

        if(!exists('dsero'))
        {
                # TODO: check there is an individual D022911 with double lastnegvd
                cols <- grep('^last|^first', names(DT), value=TRUE)
                dsero <<- unique(DT[, .SD, .SDcols=c('study_id', cols)])
                DT <- DT[, .SD, .SDcols=!cols]
        }

        if(!exists('darv'))
        {
                cols <- grep('^arv|^cuarv', names(DT), value=TRUE)
                darv <<- unique(DT[, .SD, .SDcols=c('study_id', cols)])
                DT <- DT[, .SD, .SDcols=!cols]
        }
        
        if(!exists('dcd4'))
        {
                cols <- grep('^cd4', names(DT), value=TRUE)
                dcd4 <<- unique(DT[, .SD, .SDcols=c('study_id', cols)])
                DT <- DT[, .SD, .SDcols=!cols]
        }

        if(!exists('dintdates'))
        {
                cols <- grep('^intdate', names(DT), value=TRUE) 
                dintdates <<- unique(DT[, .SD, .SDcols=c('study_id', 'hivdate', 'round', cols)])
                DT <- DT[, .SD, .SDcols=!cols]
        }

        if(!exists('dlocate'))
        {
                cols <- grep('^loc|^resident|^mobility|^hh_num', names(DT), value=TRUE) 
                cols2 <- grep('^comm|study_id|region|round', names(DT), value=TRUE) 
                dlocate <<- unique(DT[, .SD, .SDcols=c(cols2, cols)])
                DT <- DT[, .SD, .SDcols=!cols]
        }

        if(!exists('ddeath'))
        {
                cols <- grep('^icd|^fina_|death', names(DT), value=TRUE) 
                ddeath <<- unique(DT[, .SD, .SDcols=c('study_id', cols)])
                DT <- DT[, .SD, .SDcols=!cols]
        }

        if(!exists('dbirth'))
        {
                cols <- c('study_id', 'ageyrs', 'sex', 'hivdate')

                dbirth <- unique(DT[, .SD, .SDcols=cols])
                tmp <- dbirth[, uniqueN(sex) == 1, by='study_id']
                stopifnot(tmp[, all(V1)])
                cat('Median date of birth estimation\n')
                dbirth[, agedays := round(ageyrs*365.25), ]
                dbirth[, birthdate2 := hivdate - agedays]
                dbirth[, birthdate := median(birthdate2), by='study_id']
                cols <- c('study_id', 'sex', 'birthdate')
                dbirth <- dbirth[, unique(.SD), .SDcols=cols]
                stopifnot(dbirth[, uniqueN(study_id) == .N])

                dbirth <<- dbirth
        }

        return( unique(DT) )
}

study_low_level_viremia <- function(DT)
{
        tmp <- copy(DT)

        # define viremic classes
        tmp[, hiv_vl_type := ifelse(hiv_vl > 200, 'low_level_viremia', 'non_viremic')]
        tmp[hiv_vl > 1000, hiv_vl_type := 'viremic']
        tmp[is.na(hiv_vl_type), hiv_vl_type := 'unknown']
        tmp[, round := round2factor(round)]
        tmp[, hiv_vl_type := ordered(hiv_vl_type,
                                     levels=c('non_viremic', 'low_level_viremia', 'viremic', 'unknown'))]

        ggplot(tmp, aes(round, fill=hiv_vl_type)) +
                geom_histogram(position='stack', stat='count') + 
                theme_bw() + 
                theme(legend.position='bottom') + 
                labs(y='', x='Sample collection round', 
                     title='Viremia types',
                     subtitle='The proportion of low-viremic participants is decreasing with rounds'
                )

        .f <- function(remove_unknown=TRUE, remove_viremic=TRUE)
        {
                tmp1 <- tmp[, .N , by=c('round', 'hiv_vl_type')]
                tmp1 <- dcast(tmp1, round ~ hiv_vl_type, value.var='N')
                cols <- setdiff(names(tmp1), 'round')
                if(remove_unknown) cols <- setdiff(cols, 'unknown')
                if(remove_viremic) cols <- setdiff(cols, 'viremic')
                tmp1[, (cols):=lapply(.SD, na2zero) , .SDcols=cols]
                setcolorder(tmp1, c('round','non_viremic', 'low_level_viremia', 'viremic', 'unknown'))
                tmp1$TOTAL <- rowSums(tmp1[, ..cols])
                tmp2 <- tmp1[, lapply(.SD, function(x) scales::label_percent(accuracy=0.01)(x/TOTAL)) , .SDcols=cols]
                tmp2[, `:=` (round = tmp1$round, N=tmp1$TOTAL)]
                setcolorder(tmp2, 'round')
                knitr::kable(tmp2)
        }

        cat('Among suppressed participants, the prevalence of low level viremia is decreasing')
        .f()

        cat('Among participants with non-NA measurements, the prevalence of low level viremia is decreasing')
        .f(remove_unknown=TRUE, remove_viremic=FALSE)
        
}

process_darv <- function(DT=darv)
{
        darv2 <- melt(DT, id.vars='study_id')
        darv2[, `:=` (
                    round=gsub('.*?arvmed', '',variable),
                    variable=gsub('R0.*?$', '', variable)
                    )]
        darv2  <- dcast(darv2, study_id + round ~ variable )
        darv2[, round:=round2numeric(round) ]

        .f <- function(x)
        {
                y <- rep(TRUE, length(x))
                y[is.na(x)] <- FALSE
                y
        }
        cols <- c('arvmed', 'cuarvmed')
        darv2[, (cols):=lapply(.SD,.f) , .SDcols=cols]

        cat('Sanity check: if a participant reported ever taking ARV in a given round, he should be reporting the same in the following rounds.\n')
        tmp <- darv2[, list(test=!is.unsorted(arvmed)) , by='study_id']
        cat('- The number of participants reporting ARV use at least once and reporting coherent successive ARV is:\n')
        print(knitr::kable(tmp[, table(test)]))

        if(fix.incoherent.arv.reporting)
        {
                cat('- Fixing those entries.\n')
                idx <- tmp[test==FALSE, unique(study_id)]

                tmp <- darv2[study_id %in% idx, {
                        z <- which(arvmed == TRUE)[1]
                        arvmed[z:.N] <- TRUE
                        list(round=round, arvmed2=arvmed, cuarvmed2=cuarvmed)
                }, by=study_id]
                tmp[arvmed2==TRUE, cuarvmed2 := TRUE]
                setkey(tmp, study_id, round)
                setkey(darv2, study_id, round)
                darv2[tmp, `:=` (arvmed=arvmed2, cuarvmed=cuarvmed2)]
                stopifnot(darv2[, !is.unsorted(arvmed),by='study_id'][, all(V1)])
        }

        darv2[, firstarv := NA_real_]
        idx <- darv2[,any(arvmed),by='study_id'][V1==TRUE, study_id]
        darv2[study_id %in% idx, firstarv:=round[which(arvmed == TRUE)[1]], by='study_id']
        darv2
}

process_deaths <- function(DT=ddeath)
{
        DT[, cat(sum(!is.na(death_date)),'/',.N,' participants have an associated date of death\n')]
        if(fix.death.dates)
        {
                DT[!is.na(death_date) & grepl('998', death_date), death_date := '01/01/1998' ]
                DT[!is.na(death_date) , death_date := gsub('97/|98/', '01/', death_date) ]
        }
        DT[, death_date := as.Date(death_date, format='%d/%m/%Y')]
        DT[, hivdeath:=NA]
        DT[ !is.na(death_date) & grepl('HIV', fina_cause), hivdeath:=TRUE ]
        DT[!is.na(fina_cause) & !(hivdeath == TRUE), hivdeath:=FALSE]
        DT
}

count_vls_by_round <- function(DT=dvl, rnd, subset_n=NULL, quietly=TRUE)
{
        # counts the number of individuals in rounds >= rnd
        # with at least n VL counts (stratified by comm type)
        tmp <- DT[round >= rnd & !is.na(hiv_vl)]

        tmp <- participant_comm2mode(tmp, quietly=quietly) 
        tmp <- tmp[, .N ,by=c('study_id', 'comm')]
        
        if(!is.null(subset_n))
        {
                idx <- tmp[N>=subset_n, study_id] 
                # TODO: maybe consider subsetting also by round >= rnd.
                out <- DT[study_id %in% idx & !is.na(hiv_vl)]
        }

        .f <- function(n, x){sum(x>=n)}

        tmp <- tmp[, lapply(1:5, .f, x=N), by='comm']
        hd <- colSums(tmp[,-1])
        hd <- data.table(comm='Total', t(hd))
        tmp <- rbind(hd, tmp)
        
        names(tmp) <- c('Community', paste0(' >= ', 1:5))
        cat('Participants in Rounds ', rnd, ' or after with at least N viral load counts')
        print(knitr::kable(tmp))
        
        if(!is.null(subset_n)) out
}

fix_visit_dates_before_1990 <- function(DT)
{
    # DT <- copy(dvl)
    DT[, hivdate := as.Date(hivdate)]
    idx <- DT[hivdate <= '1990-01-01', {
            cat(.N, 'observation dated before 1990. Redate based on age\n'); study_id}
    ]
    tmp <- DT[study_id %in% idx, .(study_id, ageyrs, hivdate)]
    tmp <- tmp[, {
            z <- ageyrs[which(hivdate==min(hivdate))]
            tmp <- hivdate + 365*(z - ageyrs);
            list(hivdate2=mean(tmp[-which.min(hivdate)]), ageyrs=z)
    }, by='study_id']

    setkey(DT, study_id, ageyrs)
    setkey(tmp, study_id, ageyrs)

    DT[tmp, hivdate:=hivdate2]
}

followup_viremic_suppressed_by_type <- function(DT)
{
        tmp <- DT[grepl('^R15',round), ]
        # if 2 measurements in R15, R15.5, only consider the second one.
        idx <- tmp[,.N, by='study_id'][N>1, study_id]
        tmp <- tmp[!(study_id %in% idx & round=='R15-R15.5')]
        cat('The', tmp[, .N], 'visitpairs with first measurement in R15 are classified as:')
        tmp[, table_f(type, print=TRUE)]

        # for each class, see whether there will be a future viremic measurement
        tmp <- tmp[, .(study_id, hivdate_2), by='type']

        tmp1 <-  DT[, .(study_id, hivdate_2, hiv_vl_2)]
        setnames(tmp1, c('hivdate_2', 'hiv_vl_2'), c('followup_date','followup_vl'))

        tmp <- merge(tmp, tmp1, by='study_id')
        tmp1 <- tmp[followup_date > hivdate_2, 
            list(
                 type,
                 any_viremic=any(followup_vl > log10(THRESHOLD)+0.01) , 
                 any_suppressed=any(followup_vl <= log10(THRESHOLD)+0.01)
            ),
            by='study_id']

        tmp1 <- tmp1[, list(
                            any_viremic=sum(any_viremic),
                            any_suppressed=sum(any_suppressed),
                            total=.N), by='type']
        cat('\nThe number of participants successively testing viremic or suppressed at least once are:')
        knitr::kable(tmp1)
}

make_visit_pair <- function(DT)
{
        tmp <- DT[!is.na(hiv_vl), .(study_id, round, hivdate, hiv_vl)] 
        tmp[, hiv_suppressed:=TRUE]
        tmp[hiv_vl > THRESHOLD, hiv_suppressed:=FALSE]

        cols <- grep('^hiv|^round', names(tmp) ,value=TRUE) 
        newcols <- paste0(cols, '_2')
        tmp[, (newcols):=lapply(.SD, shift, n=-1L),by='study_id', .SDcols=cols]
        tmp <- tmp[!is.na(hivdate_2)]

        cols <- grep('round', names(tmp), value=TRUE)
        tmp[, (cols):=lapply(.SD, function(x){paste0('R', x)} ) , .SDcols=cols]
        tmp[, `:=`(round=paste(round, round_2, sep='-'), round_2=NULL)]

        .f <- function(x,y)
        {
                dict <- c('durably_suppressed',
                          'newly_suppressed',
                          'viral_rebound', 
                          'persistently_viremic')
                names(dict) <- as.character(c(1.5, 1, .5, 0))

                out <- as.character(x/2 + y)
                out <- as.factor(unname(dict[out]))
                return(out)
        }
        tmp[, type:=.f(hiv_suppressed, hiv_suppressed_2)]
        tmp[,  `:=` (hiv_suppressed=NULL, hiv_suppressed_2=NULL)]
        cols <- grep('hiv_vl', names(tmp),value=TRUE)
        tmp[, (cols):=lapply(.SD,function(x) log10(x+0.01)) , .SDcols=cols]

        setcolorder(tmp, c("study_id","hivdate","hivdate_2","hiv_vl","hiv_vl_2","type"))
        tmp[cat('There are ', uniqueN(study_id), ' participants accounting for',.N,' visit-pairs \n')]

        tmp1 <- tmp[, .N, by=c('round', 'type')]
        tmp1 <- dcast(tmp1, round ~ type, value.var='N')
        cols <- names(tmp1)
        tmp1[, (cols):=lapply(.SD, function(x){ x[is.na(x)] <- 0; x }) , .SDcols=cols]
        tmp1[, TOTAL:= durably_suppressed + newly_suppressed + persistently_viremic + viral_rebound]
        print(knitr::kable(tmp1))

        tmp
}

process_dcd4 <- function(DT=dcd4)
{
        # What does HTS15 mean?
        stopifnot(DT[, uniqueN(study_id) == .N])
        names(DT) <- gsub('(R|HTS)','-\\1\\2',names(DT))
        names(DT)

        idx <- c('cd4-', 'cd4date-')
        .reformat <- function(label)
        {
                cols <- grep(label, names(DT), value=TRUE)

                if(label == 'cd4-')   DT[, (cols):=lapply(.SD, as.numeric), .SDcols=cols]

                cols <- c('study_id', cols)
                tmp <- melt(DT[, ..cols], id.vars='study_id')
                tmp[, `:=` (
                            round=gsub('^.*?-', '',variable),
                            variable=gsub('-.*?$', '', variable)
                            ) ]
                setnames(tmp, 'value', tmp$variable[1])
                tmp[, variable := NULL]

                tmp
        }

        tmp <- lapply(idx, .reformat)
        tmp <- do.call('merge', tmp)

        cat('Only', tmp[round=='HTS15', sum(!is.na(cd4))],' CD4 counts with HTS15 round,removing...\n')
        tmp <- tmp[round !='HTS15']
        tmp[, round:=round2numeric(round)]
        # TODO: maybe remove NA measurements?
        tmp
}

make_joseph_table <- function()
{

        tmp <- dvl[round >= 15 & round <= 19 & !is.na(hiv_vl)]
        cols <- c('study_id', 'sex')
        tmp <- tmp[, .(VL=.N), by=cols]
        tmp <- tmp[, .N, by=c('sex', 'VL')]
        tmp <- dcast(tmp, VL~sex)

        tmp[, Total:=F+M]
        tmp
}

define_trajectories <- function(DT=dvl_15)
{
    .plot <- function(DT)
    {
            tmp <- melt(DT, id.vars=c('study_id', 'class'))
            tmp <- tmp[!is.na(value), ]
            .f <- function(x) as.integer(gsub('^V','',x))
            tmp[, `:=` (round= .f(variable), variable=NULL) ]
            
            tmp[ , v2 := 10*(1-value)+runif(.N, min=-2, max=2)]


            .f <- function(idx)
            {
                    ggplot(data=tmp[class==idx], 
                           aes(color=class, y=v2, x=round)) + 
                            geom_line(aes(group=study_id))
            }
            idx <- tmp[, unique(class)]
            plots <- lapply(idx, .f) 
            plots[[5]]

    }

    # Get wide table with round results (TRUE for suppressed and FALSE for viremic)
    cols <- c('study_id', 'hivdate', 'hiv_vl')
    tmp <- DT[, .SD,.SDcols=cols]

    setkey(tmp, study_id, hivdate)
    setorder(tmp, study_id, -hivdate)
    # stopifnot(tmp[, !is.unsorted(-hivdate), by='study_id'][, all(V1)])

    tmp[, `:=` (
        hiv_suppressed =  count2suppstatus(hiv_vl),
        visit_order = paste0('V',-seq_along(hivdate))
    ), by=study_id] 
    tmp[, `:=` (hiv_vl=NULL, hivdate=NULL)]
    tmp <- dcast(tmp, study_id ~ visit_order, value.var='hiv_suppressed')
    tmp

    # Define classes based on viremic-non-viremic measurements
    tmp[, class := NA_character_]
    tmp[`V-1` == TRUE & `V-2` == TRUE, class := 'durably_suppressed']
    tmp[`V-1` == TRUE & `V-2` == FALSE, class := 'newly_suppressed']
    tmp[`V-1` == FALSE & `V-2` == FALSE & `V-3` == FALSE, class := 'durably_viremic']
    tmp[`V-1` == FALSE & (
         `V-2` == TRUE &  (`V-3` == FALSE | `V-4` == FALSE | `V-5` == FALSE | `V-6`==FALSE) | 
         `V-2` == FALSE & `V-3` == TRUE & ( `V-4` == FALSE | `V-5` == FALSE | `V-6`==FALSE) 
        ), 
        class:='intermittently_viremic'
        ]
    .naortrue <- function(x) return(is.na(x) | x == TRUE)
    tmp[is.na(class) & `V-1` == FALSE & 
        ( .naortrue(`V-3`) & .naortrue(`V-4`) & .naortrue(`V-5`) & .naortrue(`V-6`) ),
        class:='newly_viremic']

    tmp1 <- tmp[ class=='durably_suppressed',
                all(`V-3`, `V-4`, `V-5`, `V-6`, na.rm=TRUE), by='study_id']
    tmp1[, cat('-',tmp1[, sum(!V1)],' out of ', tmp1[, .N], 'durably suppressed pariticipants',
          'had previous viremic measurements\n')]

    tmp1 <- tmp[ class=='newly_suppressed',]
    tmp2 <- tmp1[,all(!`V-3`, !`V-4`, !`V-5`, !`V-6`, na.rm=TRUE), by='study_id']
    tmp2[, cat('- For', tmp1[, sum(V1)], 'out of', tmp1[, .N], 'newly suppressed, the last measurement corresponded to the first suppressed results.\n' )]

    tmp1 <- tmp[ class=='durably_viremic',
                any(`V-4`, `V-5`, `V-6`, na.rm=TRUE), by='study_id']
    tmp1[, cat('-',tmp1[, sum(V1)],' out of ', tmp1[, .N],
               'durably viremic pariticipants previously had non-viremic measurements\n')]

    tmp1 <- tmp[ class=='intermittently_viremic', ]
    tmp1 <- tmp1[, mean(c(`V-6`,`V-5`, `V-4`, `V-3`, `V-2`, `V-1`), na.rm=T), by='study_id']
    tmp1[, cat( '- Out of ', .N, 'intermittently viremic individuals:\n -',
               sum(V1 > .5), 'historically had more viremic measurements \n -',
               sum(V1 < .5), 'historically had more non-viremic measurements\n\n')]

    # summarise
    cat('The ', tmp[, .N], 'trajectories with 4 or more VL measurements were classified as:\n')
    print( tmp[, knitr::kable(table(class))] )

    stopifnot(tmp[is.na(class) , .N == 0])
    return(tmp)
    # return(list(tmp, .plot(tmp)))
}


# From load_latest_quest.R


get.testing.history <- function(path)
{
    cols <- c('study_id', 'round', 'hiv', 'lastnegv', 'firstpos_diagnosis_vis')
    history <-  fread(path, select=cols)|>
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
    return(history)
}

get.negatives <- function(path)
{
    negatives <- fread(path,select=c('study_id', 'round')) |> unique()
    names(negatives) <- toupper(names(negatives)) 
    negatives[, ROUND:= round2numeric(ROUND)]

    if(exists('hivstatus'))
    {
        cat('\n checking negatives consistency with hivstatus\n')
        check_negs_all_in_hivstatus <- merge(negatives, hivstatus)[, sum(HIV_STATUS==0) == nrow(negatives)]
        stopifnot(check_negs_all_in_hivstatus)
        cat('All negatives can be found in `hivstatus`\n')
    }
    negatives
}

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
    cols_identifying_person_round <- c('study_id', 'round', 'comm_num')

    # load data
    cols <- c(cols_informing_hivtesting, cols_identifying_person_round)
    dtesting <- fread(path, select=cols)

    # understand labels
    # check no doubles
    idcols <- c('study_id', 'round')
    dtesting[, .N, by=idcols][, stopifnot(all(N==1))]

    # rename with labels
    dtesting[, rhivever2  := answers_dict$rhivever[as.character(rhivever)] ]
    dtesting[, hivperiod2 := answers_dict$hivperiod[as.character(hivperiod)] ]
    dtesting[, hivrslt2   := answers_dict$hivrslt[as.character(hivrslt)] ]
    dtesting[, oghivts2   := answers_dict$oghivts1[as.character(oghivts1)] ]

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
        .plot <<- function(x, q)
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

    setcolorder(dtesting, cols_identifying_person_round)

    return(dtesting)
}

make.testing.plots.community <- function(dtest, bytype=FALSE)
{
    # dtest <- copy(dtesting); bytype=TRUE
    names(dtest) <- toupper(names(dtest)) 
    if( !is.numeric(dtest$ROUND))
    {
        dtest[, ROUND := round2numeric(ROUND)]
    }

    # get comm numbers if necessary
    if(! 'COMM_NUM' %in% names(dtest) ) 
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

check_repeated_visits <- function(x, check=TRUE){
    
    dup <- x[, .N, by=c("STUDY_ID", "ROUND")][ N > 1, .(STUDY_ID, ROUND)]
    if ( ! check ) return(dup)
    stopifnot(
        "Rep entries" = nrow(dup) == 0
    )
}

check_unique_values_per_id <- function(DT, id="STUDY_ID", var){
    DT[, uniqueN(.SD), by=id, .SDcols=var][, all(V1 == 1)] |> 
        stopifnot()
}

check_increasing_dates <- function(DT, var){
    setkey(DT, STUDY_ID, ROUND)
    DT[, ! is.unsorted(.SD) , by=c("STUDY_ID"), .SDcols=var][, all(V1)] |> 
        stopifnot()
}

visit_diffs <- function(DT1, DT2, labs){
    dt1 <- subset(DT1, select=c("STUDY_ID", "ROUND"));
    dt2 <- subset(DT2, select=c("STUDY_ID", "ROUND"));
    dt1[, `:=` (FIRST = TRUE)]
    dt2[, `:=` (SEC = TRUE)]
    m <- merge(dt1, dt2, by=c("STUDY_ID", "ROUND"), all=TRUE)
    m[is.na(SEC), sprintf("There are %i visits in %s who are not in %s\n", .N, labs[1], labs[2])] |> cat()
    m[is.na(FIRST), sprintf("There are %i visits in %s who are not in %s\n", .N, labs[2], labs[1])] |> cat()

    out <- m[is.na(SEC) | is.na(FIRST)]
    setnames(out, c("FIRST", "SEC"), labs)
    return(out)
}


process_path_participation <- function(path=path.participation){

    # There are 3 ids in Melodie's list that are not here 
    #   (G115933, H123385, K122497)

    # First get all participants from the participantion file
    cols_of_int <- c("STUDY_ID", "ROUND", "INT_DATE", "COMM", "COMM_NUM", "REGION", "FIRST_PARTICIPATION", "CURR_ID" )
    allparticipants <- read_dta(path) |> setDT() 
    names(allparticipants) <- toupper(names(allparticipants))
    allparticipants[, `:=` (
        FIRST_PARTICIPATION = QST_1STRND == ROUND,
        ROUND = round2numeric(ROUND),
        COMM = comm_num2comm_type(COMM_NUM),
        INT_DATE = as.Date(INT_DATE, format="%d/%m/%Y")
    )]
    allparticipants <- subset(allparticipants, select=cols_of_int)
    check_repeated_visits(allparticipants)
    check_increasing_dates(allparticipants, "INT_DATE")
    return(allparticipants)
}

process_flow_datasets <- function(path1=path.flow.r1518, path2=path.flow.r19){ 

    cols_of_int <- c("study_id", "round", "curr_id", "sex", "comm_num", "birthdat" )
    flow <- rbind(
        {
            flow1 <- fread(path1) |> 
                subset(select=cols_of_int)
            flow1[, birthdat := as.Date(birthdat, format="%d-%b-%y")]
            flow1[ birthdat > "2030-01-01" , birthdat := as.Date(gsub("^20", "19", birthdat))]
        } ,
        {
            read_dta(path2) |> 
                as.data.table() |> subset(select=cols_of_int)
        }
    ) |> empty2NA()
    names(flow) <- toupper(names(flow))
    flow[, `:=` (
        ROUND = round2numeric(ROUND),
        COMM = comm_num2comm_type(COMM_NUM)
    ) ]
    flow <- flow[ ROUND %between% c(16, 19)]

    # CHECK SEX CONSISTENCY
    dsex <- flow[! is.na(STUDY_ID), .(STUDY_ID, SEX)] |> unique()
    # dsex |> check_unique_values_per_id(var = "SEX")
    ids_sex <- dsex[, uniqueN(SEX), by= "STUDY_ID"][V1 > 1,{
        sprintf(
            "There are %i individuals with double sex. Keep the last\n",
            .N) |> cat()
        STUDY_ID
    }]
    flow[STUDY_ID %in% ids_sex, SEX := SEX[.N]  , by="STUDY_ID"]
    # CHECK BDAY CONSISTENCY
    dbirth <- flow[! is.na(STUDY_ID), .(STUDY_ID, BIRTHDAT)] |> unique()
    # dbirth |> check_unique_values_per_id(var = "BIRTHDAT")
    idx_birth <- dbirth[! is.na(BIRTHDAT), uniqueN(BIRTHDAT), by= "STUDY_ID"][V1 > 1,{
        sprintf(
            "There are %i individuals with double birthdate. Keep the last\n",
            .N) |> cat()
        STUDY_ID
    }]
    flow[STUDY_ID %in% idx_birth, BIRTHDAT := if( diff(range(BIRTHDAT)) == 1){ 
        BIRTHDAT[1]}else{BIRTHDAT}, by=STUDY_ID]
    return(flow)
}

process_hiv_negatives <- function(path=path.negatives.r1520){
    cols_of_int <- c("STUDY_ID", "SEX", "ROUND", "INT_DATE", "HIVDATE", "COMM", "COMM_NUM", "REGION", "AGEYRS", "HIV_STATUS")
    dneg <- fread(path)
    dneg[, `:=` (
        conf_age = NULL,  # again, almost equal to ageyrs
        int_date = as.Date(int_date, '%d/%m/%Y'),
        hivdate = as.Date(hivdate, '%d%b%Y'),
        round = round2numeric(round),
        comm = comm_num2comm_type(comm_num),
        hiv_status = fifelse(hiv == "P", yes=1, no=0)
    ) ] 
    names(dneg) <- toupper(names(dneg))
    dneg <- subset(dneg, select=cols_of_int)
    dneg[, table(INT_DATE == HIVDATE, useNA = "always")]
    cat("Swapping AGEYRS for F131718, ROUND 20\n")
    dneg[ STUDY_ID == "F131718" & ROUND == 20 , AGEYRS := 20 ]
    dneg <- unique(dneg)
    check_repeated_visits(dneg)
    dneg
}

process_hiv_vls_for_hivpositives <- function(path1=path.viral.loads1, path2=path.viral.loads2){
    cols_of_int <- c("STUDY_ID", "ROUND", "INT_DATE", "HIVDATE", "COMM", "COMM_NUM", "REGION", "AGEYRS", "HIV_STATUS", "HIV_VL")
    # For some reason fread doesn't work below here:
    dvl <- read.csv(path1, fill=TRUE, comment.char="") |> 
        as.data.table() |> 
        subset(nchar(study_id) < 8) |> 
        empty2naDT() |> 
        remove.columns.with.unique.entry()
    # dvl[, round:=round2numeric(round)]
    setkey(dvl, study_id, hivdate)
    fix_visit_dates_before_1990(dvl)
    dvl <- community.keys2type(DT=dvl, by.DT='comm_num')
    dvl[, `:=` ( conf_age = NULL)]
    # Store separate info in dsero, darv, dcd4 and dintdates, dicd

    # rm(dsero,darv,dcd4,dintdates,ddeath, dlocate, dbirth)
    dvl <- make_relational_database(dvl)
    # sapply(c('dsero','darv','dcd4','dintdates','ddeath', 'dlocate', 'dbirth'), exists) |> all() |> stopifnot()
    dvl[, round := round2numeric(round)]

    tmp <- fill_na_vls_with_allpcr_data(file=path2, DT=copy(dvl))
    dvl <- tmp$dvl
    names(dvl) <- toupper(names(dvl))
    dvl$HIV_STATUS <- 1L
    return(dvl)
}
