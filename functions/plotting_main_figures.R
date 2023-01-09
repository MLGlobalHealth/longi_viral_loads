.p <- function(x)
    paste0(round(100*x, 2), '%')

dfacets <- list(
                sex = setNames(c('Male', 'Female'), c('M', 'F')),
                comm = setNames(c('Fishing', 'Inland'), c('fishing', 'inland') )
)

get.census.eligible <- function()
{
    # Recall that COUNT and TOTAL_COUNT do not agree with HIV_N and N in our vla.
    # but also recall I remove (very few) HIV+ without VLs.
    vla <- .preprocess.ds.oli(dall)
    vla <- .preprocess.make.vla(vla, c('N', 'HIV_N'))

    .load.dcens <- function(file)
    {
        dcens <- fread(file)
        setnames(dcens,
                 c('COMM', 'AGEYRS', 'SEX'),
                 paste0(c('LOC', 'AGE', 'SEX'), '_LABEL'))
        dcens <- dcens[ ! ROUND %like% '15']
        dcens[, ROUND := as.integer(ROUND)]
        cols <- grep('LABEL|ELIGIBLE$|ROUND', names(dcens), value=TRUE)
        dcens <- dcens[, ..cols]
        dcens <- merge(vla, dcens,
                        by=c('ROUND', 'LOC_LABEL', 'SEX_LABEL','AGE_LABEL'))
        dcens[, `:=` (SEX=NULL, AGE=NULL, LOC=NULL)]
        dcens
    }
    dcens <- .load.dcens(path.census.eligible)
    dcens
}

plot.all.gps <- function( loc='fishing')
{
    # rda_files <- list.files(outdir, pattern='.rda', full.names=T)

    # Plotting helpers
    # ________________

    # Cleanly format round facets

    drounds <- data.table(
                          ROUND=16:19,
                          LABS=paste0('Round ', 16:19, ':\n'),
                          START=c('07/2013', '02/2015', '10/2016','06/2018'),
                          END=c('01/2015', '09/2016', '05/2018','05/2019')
    )
    drounds[, LABS:=paste0(LABS, START, ' to ', END)]
    round_labs <- drounds$LABS
    names(round_labs) <- drounds$ROUND

    .load.draw.dfiles <- function(lab)
    {
        lab2 <- fcase(
              lab %like% 'prevalence', "prev.hiv.by.age",
              lab %like% 'suppAmongInfected', "nsinf.by.age",
              lab %like% 'suppAmongPop', "nspop.by.age",
              default=NA
        )
        stopifnot(! is.na(lab2) )

        dfiles <- grep( lab, rda_files, value=TRUE)

        # dfit contains the posterior, draw the empirical data
        dfit <- list()
        draw <- list()
        for (file in dfiles)
        {
            cat('Loading file:', basename(file), '...\n')
            tmp_env <- new.env()
            load(file, envir=tmp_env)
            rlang::env_has(tmp_env, 'DT')
            dfit[[file]] <- tmp_env[[lab2]]
            draw[[file]] <- tmp_env[["DT"]]
        }
        rm(tmp_env)
        dfit <- rbindlist(dfit, idcol='ROUND')
        draw <- rbindlist(draw)
        .f <- function(x) as.integer(gsub('^.*?round([0-9]{2}).*?$', '\\1', x))
        dfit[, ROUND := .f(ROUND)]

        stopifnot("ROUND" %in% names(draw) & "ROUND" %in% names(dfit))
        list(draw, dfit)
    }

    .plot.hiv.prevalence <- function(DFIT, DRAW)
    {
        ppDT <- copy(DRAW)
        cols <- c('M', 'CU', 'CL')
        ppDT[, (cols) := binconf(HIV_N, N, return.df=T), ]

        tmp <- DFIT[ROUND == 16, .(M, AGE_LABEL, SEX_LABEL)]
        tmp1 <- ppDT[, max(M, na.rm=T) + 0.01]

        ggplot(DFIT, aes(x=AGE_LABEL, colour=SEX_LABEL)) + 		
                geom_ribbon(aes( ymin=CL, ymax=CU, fill=SEX_LABEL), alpha=0.2, col=NA) +
                geom_line(aes( y=M)) +
                geom_line(data=tmp,aes(y=M), linetype='dotted') +
                scale_x_continuous( expand=c(0,0) ) + 
                scale_y_continuous(labels=scales:::percent, expand=c(0,0), limits=c(0,tmp1)) +
                scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
                scale_fill_manual(values=c('M'='royalblue3','F'='deeppink2')) +
                facet_wrap(~ ROUND, ncol=4, labeller=labeller(ROUND=round_labs)) +
                geom_point(data=ppDT, aes(y=M, pch=SEX_LABEL, size=N)) +
                scale_size(range = c(0, 3)) +
                theme_bw() +
                theme(legend.position='bottom') + 
                labs(x='\nage at visit (years)', 
                     y='HIV prevalence\n(95% credible interval)', 
                     pch='gender', fill='gender', color='gender',
                     size='population size'
                )
    }

    .plot.supp.hiv <- function(DFIT, DRAW)
    {
        # DFIT = dfit[LOC_LABEL==loc]; DRAW =  draw[LOC_LABEL==loc])
        ppDT <- copy(DRAW)
        cols <- c('M', 'CU', 'CL')
        ppDT[, (cols) := binconf(HIV_N - VLNS_N, HIV_N, return.df=T)]

        tmp <- DFIT[ROUND == 16, .(M, AGE_LABEL, SEX_LABEL)]
        tmp1 <- ppDT[, max(M, na.rm=T) + 0.01]

        ggplot(DFIT, aes(x=AGE_LABEL, colour=SEX_LABEL)) + 
                geom_ribbon(aes( ymin=CL, ymax=CU, fill=SEX_LABEL), alpha=0.2, col=NA) +
                geom_line(aes( y=M)) +
                geom_line(data=tmp,aes(y=M), linetype='dotted') +
                geom_hline(yintercept=c(0.9^2, 0.95^2), linetype='dashed') +
                geom_richtext(size=2.5, fill=NA, label.color=NA, color='black',
                              aes(15, 0.9^2, 
                                  label = 'UNAIDS 0.90<sup>2</sup>', 
                                  vjust = 'bottom', hjust='left')) + 
                geom_richtext(size=2.5, fill=NA, label.color=NA, color='black',
                              aes(15, 0.95^2, 
                                  label = 'UNAIDS 0.95<sup>2</sup>', 
                                  vjust = 'bottom', hjust='left')) + 
                scale_x_continuous( expand=c(0,0) ) + 
                scale_y_continuous(labels=scales:::percent, expand=c(0,0), limits=c(0,1)) +
                scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
                scale_fill_manual(values=c('M'='royalblue3','F'='deeppink2')) +
                facet_wrap(~ ROUND, ncol=4, labeller=labeller(ROUND=round_labs)) +
                geom_point(data=ppDT, aes(y=M, pch=SEX_LABEL, size=N)) +
                scale_size(range = c(0, 3)) +
                theme_bw() +
                theme(legend.position='bottom') + 
                labs(x='\nage at visit (years)', 
                     y='HIV+ individuals with suppressed viral load\n(95% credible interval)\n', 
                     colour='gender', fill='gender',
                     pch='gender', 
                     size='population size'
                )
    }

    .plot.unsupp.pop <- function(DFIT, DRAW)
    {
        ppDT <- copy(DRAW)
        cols <- c('M', 'CU', 'CL')
        ppDT[, (cols) := binconf( VLNS_N, N, return.df=T)]

        tmp <- DFIT[ROUND == 16, .(M, AGE_LABEL, SEX_LABEL)]
        tmp1 <- ppDT[, max(M, na.rm=T) + 0.01]

        p <- ggplot(DFIT, aes(x=AGE_LABEL,colour=SEX_LABEL)) + 		
            geom_ribbon(aes(ymin=CL, ymax=CU, fill=SEX_LABEL), colour=NA, alpha=0.2) +			
            geom_line(aes(x=AGE_LABEL, y=M)) +
            geom_line(data=tmp,aes(y=M), linetype='dotted') +
            geom_point(data=ppDT, aes(y=M, pch=SEX_LABEL, size=N)) +
            scale_x_continuous( expand=c(0,0) ) + 
            scale_y_continuous(labels=scales:::percent, limits=c(0,tmp1), expand=c(0,0)) +
            scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
            scale_fill_manual(values=c('M'='royalblue3','F'='deeppink2')) +
            scale_size(range = c(0, 3)) +
            facet_wrap(~ ROUND, ncol=4, labeller=labeller(ROUND=round_labs)) +
            theme_bw() +
            theme(legend.position='bottom') +
            labs(x='\nage at visit (years)', 
                 y='population with unsuppressed viral load\n(95% credible interval)\n', 
                 pch='gender', colour='gender', fill='gender',
                 size="population size")
            p
    }

    # Load HIV+ prev data
    # ___________________

    tmp <- .load.draw.dfiles( "prevalence" )
    draw <- tmp[[1]]
    dfit <- tmp[[2]]
    p1 <- .plot.hiv.prevalence(dfit[LOC_LABEL==loc], draw[LOC_LABEL==loc])

    tmp <- .load.draw.dfiles( 'suppAmongInfected')
    draw <- tmp[[1]]
    dfit <- tmp[[2]]
    p2 <- .plot.supp.hiv(dfit[LOC_LABEL==loc], draw[LOC_LABEL==loc])

    tmp <- .load.draw.dfiles( 'suppAmongPop')
    draw <- tmp[[1]]
    dfit <- tmp[[2]]
    p3 <- .plot.unsupp.pop(dfit[LOC_LABEL==loc], draw[LOC_LABEL==loc])

    names(draw)
    draw[ ROUND==19 & SEX_LABEL=="M" & AGE_LABEL %between% c(30, 35) & LOC_LABEL == loc, .(AGE_LABEL, N, HIV_N, VLNS_N, F=VLNS_N/N)]

    # Merge together

    mod <-  theme(strip.background = element_blank(), strip.text.x = element_blank()) + rremove("xlab")
    p <- ggarrange(p1 + rremove("xlab"),
              p2 + mod,
              p3 + mod, 
              align="v",
              ncol=1,
              common.legend = TRUE,
              legend='bottom'
    )
    p
}

make.table.unaids.goals <- function(age_group_width=5)
{
    # Helpers
    .p <- scales::label_percent(accuracy=.01)

    .load.dfit <- function(lab)
    {

        lab2 <- fcase(
                      lab %like% 'prevalence', "prev.hiv.by.age",
                      lab %like% 'suppAmongInfected', "nsinf.by.age",
                      lab %like% 'suppAmongPop', "nspop.by.age",
                      default=NA
        )
        stopifnot(! is.na(lab2))
        dfiles <- grep( lab, dfiles, value=TRUE)

        # dfit contains the posterior
        dfit <- list()
        for (file in dfiles)
        {
            cat('Loading file:', basename(file), '...\n')
            tmp_env <- new.env()
            load(file, envir=tmp_env)
            rlang::env_has(tmp_env, 'DT')
            dfit[[file]] <- tmp_env[[lab2]]
            # draw[[file]] <- tmp_env[["DT"]]
        }
        rm(tmp_env)
        dfit <- rbindlist(dfit, idcol='ROUND')
        .f <- function(x) as.integer(gsub('^.*?round([0-9]{2}).*?$', '\\1', x))
        dfit[, ROUND := .f(ROUND)]

        stopifnot("ROUND" %in% names(dfit))

        # Fix age groups
        dfit <- dfit[AGE_LABEL %% 1 == .5]
        dfit[ , AGEYRS := AGE_LABEL - .5 ]
        dfit[ , AGE_LABEL := NULL]

        dfit
    }

    .aggregate_by_agebins <- function(DT, width=age_group_width)
    {
        # DT <- copy(tmp); width=5
        dage <- data.table(AGEYRS = sort(unique(DT$AGEYRS)))

        tmp1 <-  range(dage$AGEYRS)
        tmp1 <- seq(from=tmp1[1], to=tmp1[2], by=width)
        dage[, DUMMY := sum(AGEYRS >= tmp1), by=AGEYRS]
        dage[, GROUP := paste0(min(AGEYRS), '-', max(AGEYRS)),by=DUMMY]
        dage[, DUMMY:= NULL]

        tmp1 <- merge(dage, DT)
        tmp1 <- tmp1[,{
            list(
                 ELIGIBLE=sum(ELIGIBLE),
                 HIV_N = sum(ELIGIBLE*PREV),
                 SUPP_N = sum(ELIGIBLE*PREV*SUPP)
            )
        } 
        , by=c('ROUND',  'LOC_LABEL', 'SEX_LABEL','GROUP')]

        tmp1[, SUPP_P:=SUPP_N/HIV_N]
        tmp1[, UNAIDS90_P := .90^2 ]
        tmp1[, UNAIDS95_P := .95^2 ]
        tmp1[, UNAIDS90_N := HIV_N * .90^2 ]
        tmp1[, UNAIDS95_N := HIV_N * .95^2 ]

        cols <- grep('_N$', names(tmp1), value=TRUE)
        tmp1[, (cols):=lapply(.SD, round) , .SDcols=cols]
        cols <- grep('_P$', names(tmp1), value=TRUE)
        tmp1[, (cols):=lapply(.SD, .p) , .SDcols=cols]

        tmp1[, UNAIDS90_DIFF := pmax(UNAIDS90_N - SUPP_N,0)]
        tmp1[, UNAIDS95_DIFF := pmax(UNAIDS95_N - SUPP_N,0)]
        setkey(tmp1, ROUND, LOC_LABEL, SEX_LABEL, GROUP)


        tmp2 <- tmp1[, list(
                            ROUND, 
                            LOC_LABEL,
                            SEX_LABEL, 
                            GROUP,
                            ELIGIBLE,
                            HIV_N,
                            SUPP=paste0(SUPP_N, ' (',SUPP_P,')'),
                            UNAIDS90=paste0(UNAIDS90_N, ' (',UNAIDS90_P,')'),
                            UNAIDS90_DIFF,
                            UNAIDS95=paste0(UNAIDS95_N, ' (',UNAIDS95_P,')'),
                            UNAIDS95_DIFF
                            )]

        p <- ggplot(data=tmp2, aes(x=GROUP)) + 
            geom_col(aes(y=UNAIDS95_DIFF, fill='95-95-95')) + 
            geom_col(aes(y=UNAIDS90_DIFF, fill='90-90-90')) + 
            scale_x_discrete(expand=c(0,0)) + 
            facet_grid(LOC_LABEL ~ SEX_LABEL, scales = 'free_y', labeller=labeller(SEX_LABEL=dfacets$sex, LOC_LABEL=dfacets$comm)) +
            theme_bw() + 
            theme(legend.position='bottom') + 
            labs(x='Age', y='Number needed to suppress', fill='UNAIDS fast-track') 

        filename <- 'main_unaidsNneeded.pdf'
        ggsave2(p, file=filename, w=18, h=12, u='cm')

        p <- ggplot(data=tmp2, aes(x=GROUP)) + 
            geom_col(aes(y=UNAIDS95_DIFF/HIV_N, fill='95-95-95')) + 
            geom_col(aes(y=UNAIDS90_DIFF/HIV_N, fill='90-90-90')) + 
            scale_x_discrete(expand=c(0,0)) + 
            facet_grid(LOC_LABEL ~ SEX_LABEL, scales = 'free_y', labeller=labeller(SEX_LABEL=dfacets$sex, LOC_LABEL=dfacets$comm)) +
            theme_bw() + 
            theme(legend.position='bottom') + 
            labs(x='Age', y='Proportion needed to suppress', fill='UNAIDS fast-track') 
        filename <- 'main_unaidsPneeded.pdf'
        ggsave2(p, file=filename, w=18, h=12, u='cm')

        tmp2
    }


    # Load HIV+ and suppAmongPop gps
    # ______________________________

    cat('Focus on last round only for this table...\n')
    dfiles <- grep(last.round,rda_files, value=TRUE)

    dprev <- .load.dfit('prevalence')
    dsuppinf <- .load.dfit('suppAmongInfected')

    # Merge with census eligible 
    # __________________________
    # to get participants ratio

    dcens1 <- copy(dcens)

    dcens1 <- dcens1[ROUND %between% c(16, 19),]
    last.round <- dcens1[, max(ROUND)]
    dcens1 <- dcens1[ROUND == last.round]
    setnames(dcens1, 'AGE_LABEL', 'AGEYRS')

    by_cols <- c('ROUND', 'LOC_LABEL', 'SEX_LABEL','AGEYRS')
    cols <- c(by_cols, 'M')

    tmp <- merge(dcens1, dprev[, ..cols], by=by_cols)

    tmp

    # tmp[, HIV_N2 := round(ELIGIBLE * M)]
    setnames(tmp, 'M', 'PREV')

    tmp <- merge(tmp, dsuppinf[, ..cols], by=by_cols)

    .f <- function(a,b, c)
    {
        y <- a*b 
        list(round(y), paste0(' (',  .p(y/a), ')'))
    }
    cols <- c('SUPP_N', 'SUPP_P')
    tmp[, (cols) := .f(ELIGIBLE*PREV , M)]
    setnames(tmp, 'M', 'SUPP')

    cols <- c('UNAIDS90_N', 'UNAIDS90_P')
    tmp[,  (cols):=.f(ELIGIBLE*PREV, .9^2)]
    tmp[, UNAIDS90_DELTA := pmax(0, UNAIDS90_N - SUPP_N)]
    cols <- c('UNAIDS95_N', 'UNAIDS95_P')
    tmp[,  (cols):=.f(ELIGIBLE*PREV, .95^2)]
    tmp[, UNAIDS95_DELTA := pmax(0, UNAIDS95_N - SUPP_N)]

    tab <- .aggregate_by_agebins(tmp)
    return(tab)
}

make.unaids.plots <- function(DT)
{
    .load.dfit <- function(lab)
    {
        # lab <- 'prevalence'

        # decide which chain to load
        lab2 <- fcase(
                      lab %like% 'prevalence', "prev.hiv.by.age",
                      lab %like% 'suppAmongInfected', "nsinf.by.age",
                      lab %like% 'suppAmongPop', "nspop.by.age",
                      default=NA
        )
        stopifnot(! is.na(lab2))

        # get rdas filenames containing the fit
        dfiles <- grep( lab, rda_files, value=TRUE)

        # dfit will store contains the posterior
        dfit <- list()
        for (file in dfiles)
        {
            cat('Loading file:', basename(file), '...\n')
            tmp_env <- new.env()
            load(file, envir=tmp_env)
            dfit[[file]] <- tmp_env[[lab2]]
            rm(tmp_env)
        }

        # Join the estimates by round 
        dfit <- rbindlist(dfit, idcol='ROUND')
        .f <- function(x) as.integer(gsub('^.*?round([0-9]{2}).*?$', '\\1', x))
        dfit[, ROUND := .f(ROUND)]
        stopifnot("ROUND" %in% names(dfit))

        # Fix age groups: do not extract "half year" ages as 24.5
        dfit <- dfit[AGE_LABEL %% 1 == .5]
        dfit[ , AGEYRS := AGE_LABEL - .5 ]
        dfit[ , AGE_LABEL := NULL]

        dfit
    }

    .aggregate_by_agebins <- function(DT, width=age_group_width)
    {
        # DT <- copy(tmp); width=5
        dage <- data.table(AGE_LABEL = sort(unique(DT$AGE_LABEL)))

        tmp1 <-  range(dage$AGE_LABEL)
        tmp1 <- seq(from=tmp1[1], to=tmp1[2], by=width)
        dage[, DUMMY := sum(AGE_LABEL >= tmp1), by=AGE_LABEL]
        dage[, GROUP := paste0(min(AGE_LABEL), '-', max(AGE_LABEL)),by=DUMMY]
        dage[, DUMMY:= NULL]

        tmp1 <- merge(dage, DT)
        cols <- c('ROUND', 'LOC_LABEL', 'SEX_LABEL','GROUP',
                  'VIR_STATUS', 'PARTICIPATION_STATUS')

        tmp2 <- tmp1[,list(
                             ELIGIBLE=sum(ELIGIBLE),
                             COUNT=sum(count)
                     ), by=setdiff(cols, 'GROUP')]
        tmp2[, GROUP :='overall']

        tmp1 <- tmp1[,list(
                              ELIGIBLE=sum(ELIGIBLE),
                              COUNT = sum(count)
                     ), by=cols]
        tmp1 <- rbind(tmp1, tmp2)
        setkeyv(tmp1, cols)
        
        # should we also exclude participation status??
        cols <- setdiff(cols, c('VIR_STATUS', 'PARTICIPATION_STATUS'))
        cols <- union(cols, 'ELIGIBLE')

        tmp1[VIR_STATUS != 'NEG', 
             .( HIV_N = sum(COUNT),
                UNAIDS_SUPP_90=sum(COUNT)*.9^3,
                UNAIDS_SUPP_95=sum(COUNT)*.95^3
             ), by=cols] -> tmp2

        tmp1[VIR_STATUS == 'SUPP',{
                z_in <- COUNT[PARTICIPATION_STATUS == 'IN']
                z_out <- COUNT[PARTICIPATION_STATUS == 'OUT']
                list(
                     POST_SUPP_HOMOGENEOUS = z_in + z_out,
                     POST_SUPP_HALF = z_in + z_out/2,
                     POST_SUPP_ZERO = z_in 
                )
        },   by=cols] -> tmp1

        dplot <- merge(tmp1, tmp2, by=cols)
        dplot
    }

    plot.unaids.goals <- function(DT)
    {
        #dplot1 <- melt(DT, id.vars=c('ROUND', 'LOC_LABEL', 'SEX_LABEL', 'GROUP'))
        #DPLOT <- copy(dplot1)

        DT <- copy(out)
        dplot1 <- out[GROUP == 'overall', {
                z <- value[variable=='HIV_N']
                list(variable=variable, 
                     value=.p(value/z)
                )},by=c('ROUND', 'LOC_LABEL', 'SEX_LABEL', 'GROUP')] 
        DT <- DT[GROUP !='overall']

        # Make labels
        dplot1 <- dplot1[variable %like% 'POST_SUPP',] 
        .s <- function(z, y=NULL)
        {
                z <- tolower(gsub('POST_SUPP_', '', z))
                if( length(y) & ! is.null(y))
                        z <- paste0(z, ': ', y)
                z

        }
        dplot1[, LAB:=.s(variable, value)]
        dplot1[, LAB := paste(c('UNAIDS: 73% 86%',LAB), collapse='\n'),
               by=c('ROUND', 'LOC_LABEL', 'SEX_LABEL', 'GROUP')]
        dplot1[, `:=` (HJUST='right', VJUST='top')]
        dplot1[, `:=` (GROUP=NULL)] 

        # dplot[variable %like% 'POST_SUPP', variable:=.s( variable) ]
        .labfill <- function(x)
        {
            x <- gsub('POST_SUPP_', '',,x) |> tolower()
            fcase(x == 'homogeneous', '',
                  x == 'half', '',
                  x == 'zero', '')
        }

        p <- ggplot(DT, aes(x=GROUP, y=value, group=1, fill=variable)) + 
                geom_col(data=DT[ variable %like% 'HOMO|HALF|ZERO'],
                         position='identity',width=1) +
                geom_col(data=DT[ variable %like% 'UNAIDS' ],
                         color='darkred', linetype='dashed', fill=NA, width=1,
                         position='identity', size=1) +
                geom_col(data=DT[ variable=='HIV_N'],
                         color='black', fill=NA, width=1, size=1) +
                geom_text(data=dplot1,
                          aes(label=LAB, hjust=HJUST, vjust=VJUST),
                          x=Inf, y=Inf ) + #, vjust=VJUST, hjust=HJUST), x=100, y=100, color='red')
                facet_grid(LOC_LABEL ~ SEX_LABEL,
                           scales = 'free_y',
                           labeller=labeller(SEX_LABEL=dfacets$sex, LOC_LABEL=dfacets$comm)) +
                # scale_fill_discrete() +
                viridis::scale_colour_viridis(discrete=TRUE, option='A', begin=.4,end = .8) +
                viridis::scale_fill_viridis(discrete=TRUE, option='A',begin=.4, end = .8, 
                                            labels= c('homogeneous among census eligible',
                                              'half-prevalence of suppression out-of-study',
                                              'no out-of-study-suppression')) +
                scale_y_continuous(expand=expansion(mult = c(0, .15))) +
                theme(legend.position='bottom', legend.direction = "vertical") +
                labs(x='Age', y='HIV positive among census eligible population',
                     fill='Suppression assumptions',
                     title='Estimated suppression among HIV positive',
                     subtitle='(Under different assumptions on non-participants)')  +
                nm_reqs  

        filename <- paste0('main_unaids_goals_N_by_sex_loc_assumptions_round',max.round,'.pdf')
        ggsave2(p, file=filename, LALA=file.path(vl.out.dir, 'figures'), w=18, h=20, u='cm')

        return(DT)
    }

    # Extract Census for last round and load posteriors
    # DT <- copy(dcens)
    # DT <- copy(dcens)
    max.round <- DT[, max(ROUND)]
    dunaids <- DT[ROUND==max.round]
    dprev <- .load.dfit('prevalence')
    dsupp_pop <- .load.dfit('suppAmongPop')

    # Get the posterior estimate for the number of
    # HIV cases & of viraemic (unsuppressed) individuals
    .f <- function(DT, colnameN=NA, colnameM=NA)
    {
        stopifnot(! is.na(colnameN) & ! is.na(colnameM))
        tmp1 <- DT[, .(ROUND, LOC_LABEL, SEX_LABEL, AGE_LABEL=AGEYRS, M)]
        merge(
              dunaids, tmp1,
              by=c('ROUND', 'LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL')
        ) -> dunaids
        dunaids[, DUMMY:=M*ELIGIBLE]
        setnames(dunaids, 'DUMMY', colnameN)
        setnames(dunaids, 'M', colnameM)
        return(dunaids)
    }
    dunaids <- .f(dprev, colnameN='HIV_N_M', colnameM='PREVALENCE_M')
    dunaids <- .f(dsupp_pop, colnameN='HIV_UNSUPP_M', colnameM='UNSUPP_PROP_M')

    # Consistency: there cannot be more unsuppressed than positive
    tmp <- dunaids[HIV_N_M < HIV_UNSUPP_M, .N]
    if(tmp > 0)
    {
        cat('Forcing posterior unsuppressed to be lower than HIV+ cases in',
            tmp, 'cases ...\n')
        dunaids[HIV_N_M < HIV_UNSUPP_M, HIV_UNSUPP_M := HIV_N_M]
    }

    # Summarise HIV status among participants and non-participants
    tmp <- dunaids[, .(ROUND, LOC_LABEL, SEX_LABEL, AGE_LABEL,
                       ELIGIBLE, N, HIV_N_M, HIV_UNSUPP_M)]
    tmp[, HIV_NEG_M := N - HIV_N_M]
    tmp[, HIV_SUPP_M := HIV_N_M - HIV_UNSUPP_M]
    tmp[, HIV_N_M := NULL]
    # Need to go wide to long
    cols <- grep('_M$', names(tmp), value=TRUE)
    tmp <- melt(tmp,  measure.vars=cols,
                value.name='PARTICIPATION_STATUS_IN', 
                variable.name='VIR_STATUS')
    tmp[, VIR_STATUS := gsub('^HIV_|_M', '', VIR_STATUS)]
    tmp[, VIR_STATUS := ordered(VIR_STATUS, c('UNSUPP', 'SUPP', 'NEG'))]
    tmp[, PARTICIPATION_STATUS_OUT := PARTICIPATION_STATUS_IN*(ELIGIBLE-N)/N]

    # Estimate distribution of viramic status among out-of-study individuals
    # TODO: (?) If we want to change HIV homogeneity assumptions, here:
    cols <- grep('PARTICIPATION_STATUS', names(tmp), value=TRUE)
    tmp <- melt(tmp, measure.vars=cols,
         value.name='count',
         variable.name='PARTICIPATION_STATUS')

    # prepare for plots 
    by_cols <- c('ROUND', 'LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL')
    tmp[, PARTICIPATION_STATUS := gsub('PARTICIPATION_STATUS_',
                                       '', PARTICIPATION_STATUS)]
    setkeyv(tmp, c(by_cols, 'PARTICIPATION_STATUS', 'VIR_STATUS' ))

    dplot <- copy(tmp)
    dplot[, cumcount:=cumsum(count), by=by_cols]
    # dplot[, PARTICIPATION_STATUS := ordered(PARTICIPATION_STATUS, c('OUT', 'IN'))]

    # ggplot label for participation rate:
    dlab <- dplot[, .(LAB=paste0('Part. rate:  \n', .p(sum(N)/sum(ELIGIBLE))),
                    VJUST='top', 
                    HJUST='right'
                    ), by=c('ROUND', 'SEX_LABEL', 'LOC_LABEL')]

    # relabelling 
    .rl1 <- function(x)
    {
        lvls <- c('HIV negative', 'suppressed viral load', 'unsuppressed viral load')
        out <- fcase(x %like% 'UNSUPP', 'unsuppressed viral load',
              x %like% '^SUPP', 'suppressed viral load',
              x %like% 'NEG', 'HIV negative')
        out <- ordered(out, levels=rev(lvls))
    }
    
    .rl2 <- function(x)
    {
        lvls <- c('out-of-study', 'in-study')
        out <- fcase(x %like% 'IN', 'in-study',
              x %like% 'OUT', 'out-of-study')
        out <- ordered(out, levels=lvls)
    }
    
    p <- ggplot(dplot[ROUND==max.round], aes(x=AGE_LABEL)) +
            geom_col(aes(y=count, alpha=.rl2(PARTICIPATION_STATUS), fill=.rl1(VIR_STATUS)), colour='grey') +
            geom_step(aes(y=N, x=AGE_LABEL - .5), position='dodge', alpha=.6, linetype='dashed') + 
            geom_step(aes(y=ELIGIBLE, x=AGE_LABEL - .5), position='dodge', alpha=.6) + 
            geom_text(data=dlab, aes(label=LAB, vjust=VJUST, hjust=HJUST), x=Inf, y=Inf)+
            facet_grid(LOC_LABEL ~ SEX_LABEL, scales = 'free_y', labeller=labeller(SEX_LABEL=dfacets$sex, LOC_LABEL=dfacets$comm)) +
            scale_y_continuous(expand=expansion(mult=c(0,.05))) +
            scale_x_continuous(expand=c(0,0)) +
            scale_alpha_discrete(range = c(.6, 1)) +
            scale_fill_manual(values=palettes$supp_status) + 
            theme(legend.position='bottom', legend.box='vertical') +
            theme(legend.key.height= unit(.25, 'cm'),
                  legend.key.width= unit(.25, 'cm'), 
                  legend.margin= unit(.4, 'cm')
                  ) + 
            guides(alpha = guide_legend(override.aes = list(colour = NULL))) + 
            guides(fill = guide_legend(override.aes = list(colour = NULL))) + 
            labs(x='Age', y='Census Eligible Population',
                 fill='HIV status', color='HIV status',
                 alpha='Participation status',
                 title='Participants among Census Eligibles Population',
                 subtitle='(assuming homogeneous population)') +
          nm_reqs
    filename <- paste0('main_hivstatus_by_participation_loc_sex_round',max.round,'.pdf')
    ggsave_nature(p, filename=filename,LALA=file.path(vl.out.dir, 'figures') , w=15, h=17)

    out <- .aggregate_by_agebins(dplot, width=5)
    out <- melt(out, id.vars=c('ROUND', 'LOC_LABEL', 'SEX_LABEL', 'GROUP')) 

    if(make_paper_numbers)
    {
        out[LOC_LABEL %like% '^f' & SEX_LABEL %like% '^F' & GROUP=='overall'] 

        # get proportion of suppressed among HIV positive, assuming...
        tmp <- out[ variable %like% 'HIV_N|HOMOGENEOUS' & GROUP =='overall']
        goal909090bysex <- tmp[, {
            z1 <- (variable == 'HIV_N')
            value[!z1]/value[z1]
        }, by=c('ROUND', 'LOC_LABEL', 'SEX_LABEL')]
        ppr_numbers[['Goals909090bysex']] <- goal909090bysex

        # aggregate over sex by community type
        tmp <- tmp[, sum(value), by=c('ROUND', 'LOC_LABEL', 'variable')]
        goal909090 <- tmp[, {
            z1 <- (variable == 'HIV_N')
            V1[!z1]/V1[z1]
        }, by=c('ROUND', 'LOC_LABEL')]
        ppr_numbers[['Goals909090']] <- goal909090

        # what should lambda (ie out of study level of suppression compared to in-study)
        # be in order to achieve goals 
        tmp <- out[ variable %like% 'HIV_N|HOMOGENEOUS|ZERO' & GROUP =='overall']
        tmp <- tmp[, .(N=sum(value)), by=c('ROUND', 'LOC_LABEL', 'SEX_LABEL', 'variable')]
        goal909090_lambda <- tmp[,{
            z.homog <- (variable %like% 'HOMOGENEOUS')
            z.zero <- (variable %like% 'ZERO')
            z.hivn <- (variable %like% 'HIV_N')
            hivn.out.of.study <- N[z.homog] - N[z.zero]
            # number needed to suppress out of study/ number
            hivn.left.to.suppress <- (N[z.hivn] * .9^3 - N[z.zero])
            required_suppression = hivn.left.to.suppress/hivn.out.of.study
        }, by=c('ROUND', 'LOC_LABEL', 'SEX_LABEL')]
        ppr_numbers[['OOSSuppressionNeededFor909090']] <- goal909090_lambda

        tmp <- out[ variable %like% 'HIV_N|HOMOGENEOUS' & GROUP !='overall']
        tmp1 <- tmp[, {
            z1 <- (variable == 'HIV_N')
            value[!z1]/value[z1] < .65
        }, by=c('ROUND', 'SEX_LABEL','LOC_LABEL', 'GROUP')]
        ppr_numbers[['Mgroupsbelow65supp']] <- tmp1[SEX_LABEL == 'M', ]

        tmp1 <- tmp[, {
            z1 <- (variable == 'HIV_N')
            value[z1]-value[!z1]
        }, by=c('ROUND', 'SEX_LABEL','LOC_LABEL', 'GROUP')]
        more_unsupp_f_than_m <- tmp1[,{
            V1[SEX_LABEL == 'F'] > V1[SEX_LABEL == 'M']
            }, by=c('ROUND','LOC_LABEL', 'GROUP')]
        ppr_numbers[['moreNunsuppInFvsM']] <- more_unsupp_f_than_m 

        ppr_numbers
    }

    plot.unaids.goals(out)
    out
}
