.p <- function(x)
    paste0(round(100*x, 2), '%')

dfacets <- list(
    sex = setNames(c('Male', 'Female'), c('M', 'F')),
    comm = setNames(c('Fishing', 'Inland'), c('fishing', 'inland') )
)

plot.all.gps <- function( loc='fishing')
{
    # Plotting helpers
    # ________________

    .plot.hiv.prevalence <- function(DFIT, DRAW, include.point=TRUE)
    {
        ppDT <- copy(DRAW)
        cols <- c('M', 'CU', 'CL')
        ppDT[, (cols) := binconf(HIV_N, N, return.df=T), ]

        tmp <- DFIT[ROUND == 16, .(M, AGE_LABEL, SEX_LABEL)]
        # y_upper_limit <- ppDT[, max(M, na.rm=T) + 0.01]
        y_upper_limit <- fifelse(loc == 'fishing', yes=.75, no=.4)

        ggplot(DFIT, aes(x=AGE_LABEL, colour=SEX_LABEL)) + 		
            geom_ribbon(aes( ymin=CL, ymax=CU, fill=SEX_LABEL), alpha=0.2, col=NA) +
            geom_line(aes( y=M)) +
            geom_line(data=tmp,aes(y=M), linetype='dotted') +
            scale_x_continuous( expand=c(0,0) ) + 
            scale_y_continuous(labels=scales:::percent, expand=c(0,0), limits=c(0,y_upper_limit)) +
            scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
            scale_fill_manual(values=c('M'='royalblue3','F'='deeppink2')) +
            facet_wrap(~ ROUND, ncol=4, labeller=labeller(ROUND=round_labs)) + {
            if(include.point) geom_point(data=ppDT, aes(y=M, pch=SEX_LABEL, size=N))
            } + 
            scale_size(range = c(0, 3)) +
            theme_default() + 
            labs(x='\nage at visit (years)', 
                y='HIV prevalence\n(95% credible interval)', 
                pch='gender', fill='gender', color='gender',
                size='population size'
            )
    }

    .plot.supp.hiv <- function(DFIT, DRAW, include.point=TRUE)
    {
        # DFIT = dfit[LOC_LABEL==loc]; DRAW =  draw[LOC_LABEL==loc])
        ppDT <- copy(DRAW)
        cols <- c('M', 'CU', 'CL')
        ppDT[, (cols) := binconf(HIV_N - VLNS_N, HIV_N, return.df=T)]
        tmp <- DFIT[ROUND == 16, .(M, AGE_LABEL, SEX_LABEL)]
        y_upper_limit <- ppDT[, max(M, na.rm=T) + 0.01]

        ggplot(DFIT, aes(x=AGE_LABEL, colour=SEX_LABEL)) + 
            geom_ribbon(aes( ymin=CL, ymax=CU, fill=SEX_LABEL), alpha=0.2, col=NA) +
            geom_line(aes( y=M)) +
            geom_line(data=tmp,aes(y=M), linetype='dotted') +
            geom_hline(yintercept=c(0.9^2, 0.95^2), linetype='dashed') +
            # geom_richtext(size=2.5, fill=NA, label.color=NA, color='black',
            #               aes(15, 0.9^2, 
            #                   label = 'UNAIDS 0.90<sup>2</sup>', 
            #                   vjust = 'bottom', hjust='left')) + 
            # geom_richtext(size=2.5, fill=NA, label.color=NA, color='black',
            #               aes(15, 0.95^2, 
            #                   label = 'UNAIDS 0.95<sup>2</sup>', 
            #                   vjust = 'bottom', hjust='left')) + 
            scale_x_continuous( expand=c(0,0) ) + 
            scale_y_continuous(labels=scales:::percent, expand=c(0,0), limits=c(0,1)) +
            scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
            scale_fill_manual(values=c('M'='royalblue3','F'='deeppink2')) + {
                if(include.point) geom_point(data=ppDT, aes(y=M, pch=SEX_LABEL, size=N)) 
            } +  
            facet_wrap(~ ROUND, ncol=4, labeller=labeller(ROUND=round_labs)) +
            scale_size(range = c(0, 3)) +
            theme_default() +
            labs(x='\nage at visit (years)', 
                 y='HIV+ individuals with suppressed viral load\n(95% credible interval)\n', 
                 colour='gender', fill='gender',
                 pch='gender', 
                 size='population size'
            )
    }

    .plot.unsupp.pop <- function(DFIT, DRAW, include.point=TRUE)
    {
        ppDT <- copy(DRAW)
        cols <- c('M', 'CU', 'CL')
        ppDT[, (cols) := binconf( VLNS_N, N, return.df=T)]

        tmp <- DFIT[ROUND == 16, .(M, AGE_LABEL, SEX_LABEL)]
        y_upper_limit <- fifelse(loc=='fishing', yes=.45, no=.15)

        p <- ggplot(DFIT, aes(x=AGE_LABEL,colour=SEX_LABEL)) + 		
            geom_ribbon(aes(ymin=CL, ymax=CU, fill=SEX_LABEL), colour=NA, alpha=0.2) +			
            geom_line(aes(x=AGE_LABEL, y=M)) +
            geom_line(data=tmp,aes(y=M), linetype='dotted') + {
                if(include.point) geom_point(data=ppDT, aes(y=M, pch=SEX_LABEL, size=N))
            } +
            scale_x_continuous( expand=c(0,0) ) + 
            scale_y_continuous(labels=scales:::percent, limits=c(0, y_upper_limit ), expand=c(0,0)) +
            scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
            scale_fill_manual(values=c('M'='royalblue3','F'='deeppink2')) +
            scale_size(range = c(0, 3)) +
            facet_wrap(~ ROUND, ncol=4, labeller=labeller(ROUND=round_labs)) +
            theme_default() +
            labs(x='\nage at visit (years)', 
                 y='population with unsuppressed viral load\n(95% credible interval)\n', 
                 pch='gender', colour='gender', fill='gender',
                 size="population size")
            p
    }

    # Load HIV+ prev data
    # ___________________

    tmp <- load.summarised.draws.from.rdafiles("prevalence", files=rda_files, include.raw=TRUE)
    data_raw <- tmp[[1]] |> subset(LOC_LABEL == loc)
    dfit <- tmp[[2]] |> subset(LOC_LABEL == loc)
    p1 <- .plot.hiv.prevalence(dfit, data_raw)
    p1.noraw <- .plot.hiv.prevalence(dfit, data_raw, include.point=FALSE)

    tmp <- load.summarised.draws.from.rdafiles("suppAmongInfected", files=rda_files, include.raw=TRUE)
    data_raw <- tmp[[1]] |> subset(LOC_LABEL == loc)
    dfit <- tmp[[2]] |> subset(LOC_LABEL == loc)
    p2 <- .plot.supp.hiv(dfit, data_raw)
    p2.noraw <- .plot.supp.hiv(dfit, data_raw, include.point=FALSE)

    tmp <- load.summarised.draws.from.rdafiles("suppAmongPop", files=rda_files, include.raw=TRUE)
    data_raw <- tmp[[1]] |> subset(LOC_LABEL == loc)
    dfit <- tmp[[2]] |> subset(LOC_LABEL == loc)
    p3 <- .plot.unsupp.pop(dfit, data_raw)
    p3.noraw <- .plot.unsupp.pop(dfit, data_raw, include.point=FALSE)

    # Merge together
    mod <-  theme(strip.text.x = element_blank()) + rremove("xlab")
    p <- ggarrange(p1 + rremove("xlab"),
              p2 + mod,
              p3 + mod, 
              align="v",
              ncol=1,
              common.legend = TRUE,
              legend='bottom'
    ) |>  annotate_figure( top=text_grob(community_dictionary[['long']][loc], size=12) )
    p.noraw <- ggarrange(p1.noraw + rremove("xlab"),
              p2.noraw + mod,
              p3.noraw + mod, 
              align="v",
              ncol=1,
              common.legend = TRUE,
              legend='bottom'
    ) |>  annotate_figure( top=text_grob(community_dictionary[['long']][loc], size=12))
    list(withraw=p, noraw=p.noraw)
}

make.table.unaids.goals <- function(age_group_width=5)
{
    # Helpers
    .p <- scales::label_percent(accuracy=.01)

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

    dprev <- load.summarised.draws.from.rdafiles('prevalence', files=dfiles , include.raw=FALSE)
    dsuppinf <- load.summarised.draws.from.rdafiles('suppAmongInfected', files=dfiles , include.raw=FALSE)

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
        # DT <- copy(out)
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
    max.round <- DT[, max(ROUND)]
    dunaids <- DT[ROUND==max.round]
    dunaids_allround <- copy(DT)
    
    dprev <- load.summarised.draws.from.rdafiles("prevalence", files=rda_files, include.raw=FALSE)
    dsupp_pop <- load.summarised.draws.from.rdafiles("suppAmongPop", files=rda_files, include.raw=FALSE)

    # Get the posterior estimate for the number of
    # HIV cases & of viraemic (unsuppressed) individuals among eligible. 
    .get.postN <- function(DT,DU=dunaids_allround, nmN=NA, nmM=NA)
    {
        stopifnot(! is.na(nmN) & ! is.na(nmM))
        tmp1 <- DT[, .(ROUND, LOC_LABEL, SEX_LABEL, AGE_LABEL=AGEYRS, M)]
        merge(
              DU, tmp1,
              by=c('ROUND', 'LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL')
        ) -> DU
        DU[, DUMMY:=M*ELIGIBLE]
        setnames(DU, 'DUMMY', nmN)
        setnames(DU, 'M', nmM)
        return(DU)
    }
    dunaids_allround <- .get.postN( dprev, nmN='HIV_N_M', nmM='PREVALENCE_M')
    dunaids_allround <- .get.postN( dsupp_pop, nmN='HIV_UNSUPP_M', nmM='UNSUPP_PROP_M')

    # Consistency:  can't be more unsupp than positive
    dunaids_allround[HIV_N_M < HIV_UNSUPP_M, 
        HIV_UNSUPP_M := {
            if(.N > 0){
                cat('Forcing posterior unsuppressed to be lower than HIV+ cases in', .N, 'cases ...\n')
                HIV_N_M
        }
    }]


    
    # Histogram of unsuppressed participants.
    #________________________________________

    dplot1 <- copy(dunaids_allround) 
    dplot1[, `:=` (
        PREVALENCE_M = NULL, 
        UNSUPP_PROP_M= NULL,
        ROW_ID=NULL
    )]
    age_breaks <- seq(15, 50, by=5)
    age_labels <- paste0(age_breaks, '-',shift(age_breaks-1, -1))[- length(age_breaks)]
    dplot1[, `:=` (
        SEX_LABEL = fifelse(SEX_LABEL == 'F', 'Female', 'Male'),
        AGE_GROUP  =  cut(
            AGE_LABEL, 
            breaks=age_breaks,
            labels=age_labels, 
            include.lowest = TRUE
        ))]
    dplot1 <- dplot1[, lapply(.SD, sum), by=c('ROUND', 'LOC_LABEL', 'SEX_LABEL', 'AGE_GROUP')]
    dplot1[, SEX_AGE_GROUP := paste(SEX_LABEL,AGE_GROUP) ]

    p1 <- dplot1 |>
        prettify_labels() |> 
        ggplot(aes(x = ROUND_LAB, group = SEX_LABEL, fill = SEX_AGE_GROUP )) + 
            geom_col( aes(y = HIV_UNSUPP_M),  position='dodge' ) + 
        facet_wrap(~LOC_LABEL) + 
        scale_fill_manual(values = palettes$sexagegroup) + 
        scale_y_continuous(expand = expansion(mult=c(0, .1))) +
        guides(fill=guide_legend(ncol = 7, nrow= 2, byrow=TRUE, key_width=.5, key_height=.5)) +
        labs( 
            x = 'Interview Round', 
            y = 'Estimated Number of unsuppressed census eligible individuals',
            fill = 'Gender and age group'
        )  + nm_reqs + 
        theme(
            legend.text = element_text(size = 5),
            legend.key.size = unit(.4, 'cm')
        )
    filename <- paste0('main_estimated_nsupp_by_loc_sex_round_age.pdf')
    ggsave_nature(p=p1, filename=filename,LALA=file.path(vl.out.dir, 'figures') , 
        w=18, h=13, add_reqs = FALSE)


    p2 <- dplot1 |>
        prettify_labels() |> 
        ggplot(aes(x = ROUND_LAB, group = SEX_LABEL, fill = SEX_AGE_GROUP )) + 
            geom_col( aes(y = VLNS_N),  position='dodge' ) + 
        facet_wrap(~LOC_LABEL) + 
        scale_fill_manual(values = palettes$sexagegroup) + 
        scale_y_continuous(expand = expansion(mult=c(0, .1))) +
        guides(fill=guide_legend(ncol = 7, nrow= 2, byrow=TRUE, key_width=.5, key_height=.5)) +
        labs( 
            x = 'Interview Round', 
            y = 'Raw Number of unsuppressed participants',
            fill = 'Gender and age group'
        )  + nm_reqs + 
        theme(
            legend.text = element_text(size = 5),
            legend.key.size = unit(.4, 'cm')
        )
    filename <- paste0('main_raw_nsupp_by_loc_sex_round_age.pdf')
    ggsave_nature(p2, filename=filename,LALA=file.path(vl.out.dir, 'figures') ,
        w=18, h=13, add_reqs = FALSE)


    # Focus on last round
    #____________________

    dunaids <- dunaids_allround[ROUND == max.round]

    # Summarise HIV status among participants and non-participants
    cols <- c('ROUND', 'LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL', 'ELIGIBLE', 'N', 'HIV_N_M', 'HIV_UNSUPP_M')
    tmp <- subset(dunaids, select=cols)
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
    # NOTE: (?) If we want to change HIV homogeneity assumptions, here:
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
    
    p <- dplot |> 
        subset(ROUND == max.round) |> 
        prettify_labels() |> 
        ggplot(aes(x=AGE_LABEL)) +
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
            labs(x='Age', 
                 y='Census Eligible Population',
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


plot.empirical.prob.of.suppression.with.age <- function(DT=dcens)
{
    # DT <- copy(dcens)
    by_vars <- c('ROUND', 'SEX_LABEL', 'LOC_LABEL')
    dcumulative <- DT[, {
        stopifnot(!is.unsorted(AGE_LABEL))
        out <- list(
            AGE_LABEL = AGE_LABEL,
            CUMULATIVE_N = cumsum(N),
            CUMULATIVE_HIV_N = cumsum(HIV_N),
            CUMULATIVE_VLNS_N = cumsum(VLNS_N),
            SUM_N = sum(N),
            SUM_HIV_N = sum(HIV_N),
            SUM_VLNS_N = sum(VLNS_N)
        )
        with(out, {
            binints0 <- Hmisc::binconf(n= CUMULATIVE_HIV_N, x=CUMULATIVE_VLNS_N, return.df = TRUE)
            binints1 <- Hmisc::binconf(n= CUMULATIVE_N, x=CUMULATIVE_VLNS_N, return.df = TRUE)
            binints2 <- Hmisc::binconf(n= SUM_VLNS_N, x=CUMULATIVE_VLNS_N, return.df = TRUE)
            binints3 <- Hmisc::binconf(n= SUM_HIV_N, x=CUMULATIVE_HIV_N, return.df = TRUE)
            names(binints0) <- paste0(c('M', 'CL', 'CU'), '_AMONG_POS')
            names(binints1) <- paste0(c('M', 'CL', 'CU'), '_AMONG_PARTICIPANTS')
            names(binints2) <- paste0(c('M', 'CL', 'CU'), '_AMONG_TOTNS')
            names(binints3) <- paste0(c('M', 'CL', 'CU'), '_AMONG_TOTHIV')
            out <- cbind(out, binints0, binints1, binints2, binints3)
        })
    }, by=by_vars] |> 
        prettify_labels()

    .plot.probability.suppression.below.age.in.group <- function(DT, .y, .ymin, .ymax, group='')
    {
        DT |>
            ggplot(aes(x=AGE_LABEL, y={{.y}}, ymin={{.ymin}}, ymax={{.ymax}}, color=SEX_LAB)) + 
                geom_point() +
                geom_linerange() +
                facet_grid( SEX_LAB + LOC_LABEL ~ ROUND, labeller = labeller(ROUND=round_labs, LOC_LABEL=community_dictionary[['short']])) + 
                theme_default() +
                scale_color_manual(values=palettes$sex, labels=sex_dictionary2) + 
                scale_y_continuous(expand=c(0,0), labels=scales::percent ) +
                labs(
                    x='Age', 
                    y=paste0('Empirical probability of being unsuppressed at age <= x (among ', group,')'),
                    color='Gender'
                )
    }

    .plot.cdf.of.age.given.group <- function(DT, .y, .ymin, .ymax, group, include.uncertainty=FALSE)
    {

        DT |> ggplot(aes(x=AGE_LABEL, y={{ .y }}, ymin= {{ .ymin }}, ymax= {{ .ymax }}, color=ROUND_LAB)) + 
            geom_line(aes( linetype = ROUND_LAB )) + {
                if(include.uncertainty) geom_linerange()
            } +
            facet_grid( LOC_LABEL ~ SEX_LAB, labeller =  labeller(LOC_LABEL=community_dictionary[['short']]) ) +
            theme_default() +
            scale_y_continuous(expand=c(0,0), labels=scales::percent ) +
            # scale_color_manual() + 
            viridis::scale_color_viridis(discrete=TRUE, begin=.2, end=1) + 
            labs(
                x='Age', 
                y=paste0('Empirical CDF of', group), 
                linetype='Interview round',
                color='Interview round'
            )
    }

    .perform.kstest.across.rounds <- function(DT, VAR)
    {
        DT[, {
            lapply(16:19, function(r) {
                rep.int(x=AGE_LABEL[ROUND == r], times=get(VAR)[ROUND == r] )
            }) -> samples 
            names(samples) <- 16:19
            out <- CJ(unique(ROUND), unique(ROUND))[V1 < V2]
            out[, `:=` (IDX1=V1-15, IDX2=V2-15)]
            out[, KS.PVALUE := ks.test(samples[[IDX1]], samples[[IDX2]])$p.value, by=c('IDX1', 'IDX2') ]
            out
        }, by=c('SEX_LABEL', 'LOC_LABEL')]
    }

    p1 <- .plot.probability.suppression.below.age.in.group(dcumulative, .y=M_AMONG_POS, .ymin=CL_AMONG_POS, .ymax=CU_AMONG_POS, group='HIV positive participants')
    filename <- 'prob_supp_belowage_among_hivp.pdf'
    ggsave2(p1, file=filename, LALA=vl.out.dir.figures , w=25, h=20, u='cm')

    p2 <- .plot.probability.suppression.below.age.in.group(dcumulative, .y=M_AMONG_PARTICIPANTS, .ymin=CL_AMONG_PARTICIPANTS, .ymax=CU_AMONG_PARTICIPANTS, group='all participants')
    filename <- 'prob_supp_belowage_among_parts.pdf'
    ggsave2(p2, file=filename, LALA=vl.out.dir.figures , w=25, h=20, u='cm')

    p3 <- .plot.cdf.of.age.given.group(dcumulative, .y=M_AMONG_TOTNS, group='non-suppression')
    filename <- 'ecdf_age_nonsupp.pdf'
    ggsave2(p3, file=filename, LALA=vl.out.dir.figures , w=25, h=20, u='cm')

    p4 <- plot.cdf.of.age.given.group(dcumulative, .y=M_AMONG_TOTHIV,  group='HIV positives')
    filename <- 'ecdf_age_positives.pdf'
    ggsave2(p4, file=filename, LALA=vl.out.dir.figures , w=25, h=20, u='cm')


    # however this does not account for prevalence and participation.
    # What is the probability of a participant to be suppressed 

    # Kolmo-Smirno to test temporal differences in these distributions.
    ks.vlns <- .perform.kstest.across.rounds(DT=DT, VAR='VLNS_N')[, GROUP := 'VLNS_N']
    ks.hivn <- .perform.kstest.across.rounds(DT=DT, VAR='HIV_N')[, GROUP := 'HIV_N']
    ks.part <- .perform.kstest.across.rounds(DT=DT, VAR='N')[, GROUP := 'N']
    ks.results <- rbind(ks.vlns, ks.hivn, ks.part)

    if(0)
    {
        ks.vlns[, table(KS.PVALUE < .05)]
        ks.vlns[KS.PVALUE < .05]
        ks.hivn[, table(KS.PVALUE < .05)]
        ks.hivn[KS.PVALUE < .05]
        ks.part[, table(KS.PVALUE < .05)]
        ks.part[KS.PVALUE < .05]
    }

    list(
        p_supp_hiv=p1, 
        p_supp_parts=p2, 
        p_cdf_supp=p1, 
        p_cdf_hiv=p2, 
        ks=ks.results
    )
}

plot.rakai.map <- function(.size=3, labs=FALSE){

    require(ggmap)
    require(maps)
    require(ggsn)

    # load comm ids
    comms <- fread(path.community.idx)
    names(comms) <- toupper(names(comms))
    comms[, COMM_NUM := as.character(COMM_NUM)]

    # load gps and merge to comms
    env_gps <- new.env()
    load(path.community.gps, envir=env_gps)
    comms <- with(env_gps$comgps, {
        data.table(
            COMM_NUM=as.character(COMM_NUM),
            longitude=as.numeric(longitude),
            latitude=as.numeric(latitude)
    )}) |> merge(x=comms,y=_, by="COMM_NUM") 
    comms[, TYPE := fifelse(COMM_ID %like% 'i', 'inland', 'fishing') ]

    # if 2 comms have the same COMM_ID, take average of lat/long
    comms <- comms[, .(
        TYPE = unique(TYPE),
        latitude = mean(latitude),
        longitude = mean(longitude)
    ), by=COMM_ID]

    # get population sizes.
    tmp <- readRDS(path.comm.censsize)
    comms <- merge(comms, tmp, by.x=c('COMM_ID', 'TYPE'), by.y=c('COMM_IDX', 'TYPE'))

    with(comms, c(
            left=min(longitude) - .05,
            right=max(longitude) + .05,
            bottom=min(latitude) - .05,
            top=max(latitude) + .05
    )) -> box

    ggmap::get_stamenmap(
        bbox = box,
        maptype = 'watercolor',
        zoom = 11
    ) -> map

    .breaks <- list(NULL, waiver)[[labs + 1]]
    .xlab <- list(NULL, "Longitude (°E)")[[labs + 1]]
    .ylab <- list(NULL, "Latitude (°N)")[[labs + 1]]

    prettify_labels(comms)
    p <- ggmap(map) + 
        geom_point(data=comms, 
            aes(x=longitude, y=latitude, 
                shape=LOC_LAB, # size=POP_SIZE,
                color=LOC_LAB, fill=LOC_LAB
            ), size=.size) +
        scale_shape_manual(values=c('inland'=21, 'fishing'=25, "Inland"=21, "Fishing"=25)) +
        scale_color_manual(values=c('black', 'black')) +
        scale_fill_manual(values=palettes$comm) +
        scale_y_continuous(breaks = .breaks, expand=c(0,0)) +
        scale_x_continuous(breaks = .breaks, expand=c(0,0)) +
        my_labs(color=NULL, fill=NULL, shape=NULL, x=.xlab, y=.ylab) +
        theme(
            legend.position = c(1,1),
            legend.justification = c(1,1),
            legend.background = element_rect(fill='white', color='black'),
            legend.direction = 'horizontal'
        )  +
        ggsn::scalebar(
            x.min=box['left'] + .02,
            x.max=box['right'] - 0.01, 
            y.min=box['bottom'] + .03,
            y.max=box['top'] - .01 ,
            transform = TRUE, 
            dist=10, dist_unit = 'km',
            location = 'bottomleft', 
            border.size = .1,
            height=0.01,
            st.bottom = FALSE,
            st.dist=0.04,
            st.size=2)
    p
}

plot.relative.suppression.vs.round16.ratio <- function(DT= dincreasessupp)
{
    
    dplot <- copy(DT)
    prettify_labels(dplot)

    ggplot(dplot, aes(x=AGEYRS, color=SEX_LAB, fill=SEX_LAB)) + 
        geom_hline(yintercept=1, linetype=2, color='grey80') +
        geom_hline(yintercept=4, linetype=2, color='grey80') +
        geom_ribbon(aes(ymin=CL, ymax=CU), alpha=0.2, color=NA) +
        geom_line(aes(y=M)) +
        facet_grid(LOC_LAB  ~ ROUND_LAB, scales="free_y", labeller= labeller(ROUND_LAB = round_labs) ) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(trans = scales::log2_trans(),
            breaks = scales::trans_breaks(paste0("log2"), function(x) 2^x),
            labels = scales::trans_format(paste0("log2"), scales::label_math(expr = 2^.x))
        ) +
        scale_color_manual(values=palettes$sex, labels=sex_dictionary2) + 
        scale_fill_manual(values=palettes$sex, labels=sex_dictionary2) + 
        theme_default() +
        my_labs(x='', y='Posterior ratio of viral suppression among PLHIV relative to Round 16') +
        NULL
}

plot.relative.suppression.vs.round16.diff <- function(DT = dincreasessupp_diff)
{

    dplot <- copy(dincreasessupp)
    prettify_labels(dplot)

    ggplot(dplot, aes(x=AGEYRS, color=SEX_LAB, fill=SEX_LAB, linetype=LOC_LAB)) + 
        geom_hline(yintercept=1, linetype=2, color='grey80') +
        # geom_ribbon(aes(ymin=CL, ymax=CU), alpha=0.2, color=NA) +
        geom_line(aes(y=CL)) +
        facet_grid(~ ROUND_LAB, scales="free_y", labeller= labeller(ROUND_LAB = round_labs) ) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_log10(
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        scale_color_manual(values=palettes$sex, labels=sex_dictionary2) + 
        scale_fill_manual(values=palettes$sex, labels=sex_dictionary2) + 
        theme_default() +
        my_labs(x='', y='Posterior ratio of viral suppression among PLHIV relative to Round 16') +
        NULL
}


plot_2yaxis_hist_lines <- function(DThist, DTline, sec_name="Contribution to HIV prevalence"){

    ALPHA = .5; DODGE = 1
    
    # rescale 2nd data frame for the secondary axis
    .sec_axis_scale <-  max(DThist$CU) / max(DTline$CU) 
    by_cols <- DThist[, sapply(.SD, function(x) typeof(x) == "character")] |> which() |> names()
    DThist_sa <- DThist[, lapply(.SD, function(x) x/.sec_axis_scale), .SDcols = c("M", "CL", "CU"), by=c(by_cols, "AGEYRS")]
    ggplot(DTline, aes(x=AGEYRS, y=M, ymin=CL, ymax=CU, fill=SEX_LAB, color=SEX_LAB)) +
        geom_col(data=DThist_sa, alpha=ALPHA, color='grey80', linewidth=.2 , position=position_dodge(width = DODGE))+
        geom_linerange(data=DThist_sa, position=position_dodge(width = DODGE)) +
        geom_line() +
        geom_ribbon(alpha=ALPHA, color=NA) +
        scale_y_continuous(
            labels = scales::label_percent(), 
            expand = expansion(mult = c(0, .1)),
            sec.axis = sec_axis(
                trans= ~ .*.sec_axis_scale,
                labels = scales::label_percent(),
                name=sec_name) 
        ) +
        scale_x_continuous(expand = c(0,0), breaks= c(seq(15, 45, 5), 50)) + 
        scale_fill_manual(values=palettes$sex, labels=sex_dictionary2) +  
        scale_color_manual(values=palettes$sex, labels=sex_dictionary2) +  
        facet_grid( . ~ LOC_LAB )  +
        my_labs(y = "HIV prevalence by age") +
        theme_default() +
        NULL
}

plot_prevalenceandcontrid <- function(
    DTprev,
    DTcontrib,
    sec_name="Contribution to HIV prevalence in age group",
    legend.key.size=unit(0.5, 'cm'),
    sec_axis_scale = NA_real_,
    extra_fig = NULL,
    CrI=TRUE,
    slides=FALSE){

    ALPHA = .5; DODGE = .5

    n_loc <- DTprev[, uniqueN(LOC)]
    n_sex <- DTprev[, uniqueN(SEX)]

    REQS <- if(slides==TRUE){
        slides_reqs
    }else{
        nm_reqs
    }

    .plot.facet <- function(DT){

        naturemed_reqs()

        # if empty data table do not plot:
        if(nrow(DT) == 0) return(NULL)

        # plot settings
        .loc_lab <- unique(DT$LOC_LAB)
        .sex_lab <- unique(DT$SEX_LAB)
        .ymax <- c(.38, .65)[.loc_lab %like% "Fishing" + 1]
        .sec_axis_scale <- fcoalesce( 
            sec_axis_scale,
            0.035 /.ymax
        )
        .M.intcode <- as.integer(.sex_lab == "Male") + 1

        # rescale 2nd y-axis
        DT[ LABEL == "Contribution to HIV prevalence", c("M", "CL", "CU") := {
            stopifnot("probably wrong label"=.N > 0)
             .(M / .sec_axis_scale, CL / .sec_axis_scale, CU / .sec_axis_scale )
        }]

        # relabel
        tmp_labs <- c(
            `HIV prevalence`="HIV prevalence in age group",
            `Contribution to HIV prevalence`= sec_name
        )
        DT[, LABEL2 := tmp_labs[ LABEL ] ]
        DT[, LABEL2 := factor(LABEL2, levels=tmp_labs, ordered=TRUE)]
        DTline <- subset(DT, LABEL == "Contribution to HIV prevalence")
        DThist <- subset(DT, LABEL == "HIV prevalence")
        .breaks <- c(unique(DThist$LABEL2), unique(DTline$LABEL2))
        # return(DThist[,range(M)])


        .facet_formula <- if(n_loc == 1){
            as.formula(LOC_LAB ~ SEX_LAB)
        }else{
            as.formula(LOC_LAB ~ SEX_LAB)
        }
        
        # rescale 2nd data frame for the secondary axis
        ggplot(DT, aes(x=AGEYRS, y=M, ymin=CL, ymax=CU, fill=LABEL2, color=LABEL2)) +
            geom_line( data=DTline) + 
            geom_col(
                data=DThist, 
                alpha=ALPHA, color=NA,
                position=position_dodge(width = DODGE)
            ) + { 
                if(CrI){ geom_linerange(position=position_dodge(width=DODGE), linewidth = .2) } else {NULL}
            } +
            scale_y_continuous(
                labels = scales::label_percent(), 
                limits=c(0, c(.38, .68)[.loc_lab %like% 'Fishing' + 1] ),
                expand = expansion(mult = c(0, 0)),
                sec.axis = sec_axis(
                    trans= ~ .*.sec_axis_scale,
                    labels = scales::label_percent(),
                    name=list(NULL, sec_name)[[.M.intcode]] )
            ) +
            scale_x_continuous(expand = c(0,0), breaks= c(seq(15, 45, 5), 50)) + 
            scale_fill_manual(values=rev(palettes$minimal2), breaks = .breaks)  +
            scale_color_manual(values=rev(palettes$minimal2), breaks= .breaks)  +
            facet_wrap( .facet_formula, 
                labeller=labeller(LOC_LAB=community_dictionary$longest, SEX_LAB=sex_dictionary2, .multi_line=FALSE)
            )  +
            my_labs(
                y = list("HIV prevalence by age", NULL)[[.M.intcode]], 
                x=NULL, 
                color="",
                fill=""
            ) +
            theme_default(
                strip.placement = "outside",
                legend.key.size=legend.key.size, 
                legend.spacing.x = legend.key.size,
                plot.margin=margin(t=0,b=0,l=0,r=0, unit="cm")
            ) +
            REQS
    }

    # bind
    DTprev[, LABEL := "HIV prevalence"]
    DTcontrib[, LABEL := "Contribution to HIV prevalence"]
    DT <- rbind(DTprev, DTcontrib)
    prettify_labels(DT)

    facets_dts <- list(
        subset(DT, LOC_LAB == "Inland" & SEX_LAB == "Female"),
        subset(DT, LOC_LAB == "Inland" & SEX_LAB == "Male"),
        subset(DT, LOC_LAB == "Fishing" & SEX_LAB == "Female"),
        subset(DT, LOC_LAB == "Fishing" & SEX_LAB == "Male")
    )
    # remove facets in case they are not needed
    ps <- lapply(facets_dts, .plot.facet)

    ps <- ps[! sapply(ps, is.null)]

    if(length(ps) == 4){
        p <- ggarrange(
            ps[[1]], ps[[2]],
            ps[[3]], ps[[4]],
            ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
        p
    }else{
        p <- ggarrange(
            ps[[2]] + theme(strip.text.y = element_blank(), axis.title.y.right = element_blank()),
            ps[[1]] + theme(strip.text.y = element_blank()) + labs(y=""), 
            nrow=1, common.legend = TRUE, legend="bottom")
        N <- fifelse(is.null(extra_fig), 1, 2)
        p <- ggarrange(extra_fig, p, nrow=1, ncol=N, widths=c(1, 2)[1:N])
    }
    return(p)
}

plot_suppandcontrib <- function(
    DTprev1,
    DTcontrib1, 
    sec_name="Contribution to viraemic population",
    prevalence.label="FILL IT!",
    viraemia.label="Contribution to viraemic population",
    legend.key.size=unit(0.4, 'cm'),
    remove.legend=FALSE,
    sec_axis_scale=NA_real_,
    slides=FALSE,
    CrI=TRUE,
    UNAIDS=TRUE
){
    if(UNAIDS){
        require(geomtextpath)
    }

    ALPHA = .5; DODGE = 1

    n_loc <- DTprev1[, uniqueN(LOC)]
    n_sex <- DTprev1[, uniqueN(SEX)]
    REQS <- c(nm_reqs, slides_reqs)[slides==TRUE + 1]

    .plot.facet <- function(DT){

        # if empty data table do not plot:
        if(nrow(DT) == 0) return(NULL)

        naturemed_reqs()


        .sec_axis_scale <- fcoalesce( 
            sec_axis_scale,
            DT[, { max(CU[LABEL == viraemia.label]) / max(CU[LABEL == prevalence.label]) },]
        )
        cat(.sec_axis_scale, "\n")


        DT[ LABEL == viraemia.label, c("M", "CL", "CU") := {
            stopifnot("probably wrong label"=.N > 0)
             .(M / .sec_axis_scale, CL / .sec_axis_scale, CU / .sec_axis_scale )
        }]

        # prettify labels
        tmp_labs <- c(
            paste(prevalence.label),
            sec_name
            # "Distribution of individuals with viraemic viral load"
        )
        names(tmp_labs) <- c(prevalence.label, viraemia.label)
        DT[, LABEL2 := tmp_labs[LABEL] ]
        DT[, LABEL2 := factor(LABEL2, levels=tmp_labs, ordered=TRUE)]

        # conditional faceting order
        facet_formula <- if(n_loc == 1){
            as.formula(LOC_LAB ~ SEX_LAB)
        }else{
            as.formula(SEX_LAB + LOC_LAB ~ .)
        }

        # rescale 2nd data frame for the secondary axis
        ggplot(DT, aes(x=AGEYRS, y=M, ymin=CL, ymax=CU, fill=LABEL2)) +
            geom_col(alpha=ALPHA, color='grey80', position=position_dodge(width = DODGE))+
            { if(CrI){ geom_linerange(position=position_dodge(width = DODGE)) } else {NULL} } +
            scale_y_continuous(
                labels = scales::label_percent(), 
                limits = c(0, c(.55, .77)[as.integer(CrI == TRUE) + 1] ),
                expand = expansion(mult=c(0,0)),
                sec.axis = sec_axis(
                    trans = ~ .*.sec_axis_scale,
                    labels = scales::label_percent(),
                    name = sec_name) 
            ) + {
                if(UNAIDS){
                    geom_texthline(
                        yintercept=1-.95^3, color='red', linetype='dashed', 
                        label="UNAIDS 95-95-95", vjust=0.5, hjust=1) 
                }else{ NULL }
            } +
            scale_x_continuous(expand = c(0,0), breaks= c(seq(15, 45, 5), 50)) + 
            scale_fill_manual(values=palettes$minimal3)  +
            scale_color_manual(values=palettes$minimal3)  +
            facet_grid( facet_formula, labeller=labeller(LOC_LAB=community_dictionary$longest, SEX_LAB=sex_dictionary2) )  +
            my_labs(y = prevalence.label, x="", color="", fill="") +
            theme_default(
                strip.placement = "outside",
                legend.key.size=legend.key.size, 
                legend.spacing.x = legend.key.size
            ) +
            REQS 

    }

    DTprev <- copy(DTprev1)
    DTcontrib <- copy(DTcontrib1)

    # bind
    DTprev[, LABEL := prevalence.label]
    # DTprev[, (c('M', 'CL', 'CU')) := lapply(.SD, function(x) 1-x), .SDcols=c('M', 'CL', 'CU')]
    DTcontrib[, LABEL := viraemia.label]
    DT <- rbind(DTprev, DTcontrib)
    prettify_labels(DT)

    facets_dts <- list(
        subset(DT, LOC_LAB == "Inland" & SEX_LAB == "Female"),
        subset(DT, LOC_LAB == "Inland" & SEX_LAB == "Male"),
        subset(DT, LOC_LAB == "Fishing" & SEX_LAB == "Female"),
        subset(DT, LOC_LAB == "Fishing" & SEX_LAB == "Male")
    )
    # remove facets in case they are not needed
    ps <- lapply(facets_dts, .plot.facet)
    ps <- ps[! sapply(ps, is.null)]

    if(length(ps) == 4){
        p <- ggarrange(
            ps[[1]], ps[[2]],
            ps[[3]], ps[[4]],
            ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
        p
    }else{
        .theme <- function(x){
            theme(
                strip.text.y=element_blank(),
                strip.text.x=element_text()
            )
        }
        p <- ggarrange(
            ps[[2]] + theme(strip.text.y = element_blank(), axis.title.y.right = element_blank()),
            ps[[1]] + theme(strip.text.y = element_blank()) + labs(y=""), 
            ncol=2, nrow=1, common.legend = TRUE, legend=c("bottom", "none")[remove.legend + 1])
    }

    return(p)
}

plot.uganda.map <- function(zoom="medium", maptype='toner-lite', labs=TRUE){

    require(ggmap)

    stopifnot(zoom %in% c("far", "medium","close"))

    .xlab <- list(NULL, "Longitude (°E)")[[labs + 1]]
    .ylab <- list(NULL, "Latitude (°N)")[[labs + 1]]

    box_uganda <- if( zoom == "far" ){
        c(top=6.158199, left=24.681059, bottom=- 6.314563, right=42.125359)
    }else if( zoom == "medium" ){
        c(top=1.5, left=29, bottom= -8, right=40)
    }else if( zoom == "close" ){
        c(top=NA_real_, left=NA_real_, bottom=NA_real_, right=NA_real_)
    }

    env_gps <- new.env()
    load(path.community.gps, envir=env_gps)
    box <- with(env_gps$comgps, {
        data.table(top=max(latitude), left=min(longitude), bottom=min(latitude), right=max(longitude))
    }) 
    
    uganda <- get_stamenmap( 
        bbox=box_uganda, 
        maptype = maptype,
        source = "stamen",
        zoom = 6,
    )
    ggmap(uganda)  +
        geom_polygon(data = map_data("lakes"), aes(x = long, y = lat, group = group), fill = "blue", alpha 
= 0.2) +
        geom_rect(xmin=box$left, xmax=box$right, ymin=box$bottom, ymax=box$top, fill=NA, color='red') + 
        theme_classic() +
        labs( x=.xlab ,y=.ylab)
}

plot.comparison.ftp.nonftp.and.all <- function(env=env_list, DT=djoint, model="run-gp-supp-hiv", round=19, y_upper_limit=NA){
    

    ALPHA = .3
    cols <- c("LOC_LAB", "SEX_LAB", "AGEYRS", "M", "CL", "CU", "FTP_LAB")
    model_lab <- fcase(
        model == "run-gp-supp-hiv", "nsinf.by.age",
        model == "run-gp-supp-pop", "nspop.by.age",
        model == "run-gp-prevl", "prev.hiv.by.age",
        default = NA_character_)
    ylab <- fcase(
        model == "run-gp-supp-hiv", "Prevalence of viral suppression among HIV positive individuals",
        model == "run-gp-supp-pop", "Prevalence of viraemia among population",
        model == "run-gp-prevl", "Prevalence of HIV",
        default = NA_character_
    )
    stopifnot( "Unknown model label" = !is.na(model_lab) )

    # subset to data of interest
    comp <- rbind.ftpstatus.datatables.by.round(model_lab, round, envir_list = env) |> 
        prettify_labels() |> 
        setnames("AGE_LABEL", "AGEYRS") |>
        subset(select=cols, AGEYRS %between% c(15, 49))
    merged <- djoint |> 
        subset(MODEL == model & ROUND == round) |> 
        prettify_labels() |> 
        subset(select=setdiff(cols, "FTP_LAB")) 
    merged[, FTP_LAB := "Weighted average"]


    # make plot
    p <- rbind(comp, merged) |> ggplot(aes(x=AGEYRS, y=M, ymin=CL, ymax=CU, linetype=FTP_LAB)) + 
        geom_ribbon(data=merged, color=NA, alpha=ALPHA) +
        geom_line() +
        scale_y_continuous(labels=scales:::percent, expand=c(0,0), limits=c(0, y_upper_limit)) +
        scale_x_continuous(expand = c(0,0), breaks= c(seq(15, 45, 5), 50)) + 
        scale_linetype_manual(values=c(2,3,1)) + 
        facet_grid(SEX_LAB ~ LOC_LAB) +
        my_labs(y=ylab) +
        theme_default() + 
        NULL
    return(p)
}

plot.comparison.prevalence.fishinginland.oneround <- function(DT, model, round, ylim=NA){

    require(geomtextpath)
    stopifnot(model %in% names(model_dictionary) )

    dplot <- subset(DT, MODEL == model & ROUND == round) |> 
        prettify_labels()
    dplot[ , LOC_LAB := community_dictionary$longest[LOC_LAB] ]

    ggplot(dplot, aes(x=AGEYRS, y=M, ymin=CL, ymax=CU, fill=LOC_LAB, color=LOC_LAB)) + 
        geom_line() +
        geom_ribbon(alpha=ALPHA, color=NA) + {
            if(MODEL == 'run-gp-supp-hiv'){
                geom_texthline(
                    yintercept=.95^3, color='red', linetype='dashed', 
                    label="UNAIDS 95-95-95", vjust=0.5, hjust=1) 
            }else{ NULL }
        } +
        facet_grid(~ SEX_LAB) +
        scale_color_manual(values=palettes$comm)  +
        scale_fill_manual(values=palettes$comm)  +
        scale_x_continuous(expand = c(0,0), breaks= c(seq(15, 45, 5))) + 
        scale_y_continuous(labels=scales::label_percent(), limits=c(0,ylim), expand=c(0,0)) + 
        theme_default() +
        my_labs(y=model_dictionary[model])
}

plot.main.suppression.among.plhiv <- function(DT=djoint, hist=TRUE, unaids=FALSE, rev){

    .ALPHA=.3; .HIST=TRUE; .LINEWIDTH=.2
    .ylab <- fifelse(rev, 
        yes= "Prevalence of viraemia in PLHIV by age group",
        no= "Prevalence of viral suppression in PLHIV by age group",
    )
    dplot <- subset(DT, ROUND %in% c(16,19) & MODEL == 'run-gp-supp-hiv')
    dplot <- prettify_labels(dplot)
    if(rev){
        dplot[, `:=` (M = 1-M, CL=1-CU, CU=1-CL)]
    }

    .main <- function(.hist=hist, .unaids=unaids){
        .pd <- position_dodge(width=1)
        .yint <- fifelse(rev, yes=1-.95^3, no=.95^3)
        out <- list(
            if(!hist){ geom_ribbon(alpha=.ALPHA, color=NA) },
            if(!hist){ geom_line(linewidth = .LINEWIDTH) },
            if(hist){ geom_col(position=.pd) },
            if(hist){ geom_linerange(position=.pd, color='grey40') },
            if(unaids) {
                geomtextpath::geom_texthline(
                    yintercept = .yint,
                    color = 'purple', 
                    linetype='dashed', 
                    label = 'UNAIDS 95-95-95', 
                    vjust=.5, hjust=fifelse(rev, yes=1, no=0),
                    linewidth=.LINEWIDTH,   
                    size=2.5
                )
            },
            NULL
        )
        out[! sapply(out, is.null)]
    }
    out <- .main(TRUE)

    ggplot(dplot, 
        aes( x=AGEYRS, y=M, ymin=CL, ymax=CU, fill=SEX_LAB, color=SEX_LAB)
    ) + 
        .main() +
        facet_grid( 
            LOC_LAB ~ ROUND_LAB, 
            labeller=labeller(ROUND_LAB=round_labs, LOC_LAB=community_dictionary$longest),
        ) +
        scale_y_percentagef(.02) + 
        scale_x_continuous(breaks=seq(15, 60, 5), expand=c(0,0)) + 
        scale_color_manual(values=palettes$sex, breaks = sex_dictionary2) +
        scale_fill_manual(values=palettes$sex, breaks= sex_dictionary2) +
        theme_default() +
        my_labs(y=.ylab) 
}

plot_propofpop_of_viraemic_byagesex_stratbycommround <- function(DT, colorby="ROUND", cri=TRUE){

    # helpers for plot flexibility
    colorby <- match.arg(colorby, c("ROUND_LAB", "SEX_LAB")) |> invisible()
    expr_col <- ifelse(colorby=="ROUND_LAB", expr(ROUND_LAB), expr(SEX_LAB))
    .scales <- function(x={colorby=="ROUND_LAB"}){
        out <- list( 
            scale_y_continuous(expand=expansion(mult=c(0,.1)), labels=scales::label_percent()), 
            scale_x_continuous(expand=c(0,0)), 
            if(!x){ scale_color_manual(values=palettes$sex, labels=sex_dictionary2)}, 
            if(!x){ scale_fill_manual(values=palettes$sex, labels=sex_dictionary2)}, 
            if(x){ scale_fill_viridis_d(begin=0, end = .7, option = "magma" )}, 
            if(x){ scale_color_viridis_d(begin=0, end = .7, option = "magma" )}, 
            NULL
        )
        out[ ! sapply(out, is.null) ]
    }

    # data manipulation: 
    dplot <- subset(DT, MODEL == 'run-gp-supp-pop') |>
        set(j = c('IL', 'IU'), value = NULL) |> 
        merge(dcens[, - c("ELIGIBLE", "AGEGROUP") ], by=c("ROUND", "LOC", 'SEX', "AGEYRS")) |> 
        prettify_labels() 
    cols <- c("M", "CL", "CU")
    dplot[, 
        (cols) := lapply(.SD, function(x, y) x*y, y=proportions(ELIGIBLE_SMOOTH)),
        .SDcols = cols,
        by=c("ROUND", "LOC_LAB")
    ]

    ggplot(dplot,   
        aes(x=AGEYRS, y=M, ymin=CL, ymax=CU,
            fill=eval(expr_col), color=eval(expr_col),
            linetype=ROUND_LAB, group=ROUND_LAB)
        ) + { if(cri){
            geom_ribbon(linetype=0, alpha=.1)
        }else{NULL} } +
        geom_line() +
        theme_default() + 
        facet_grid(
            LOC_LAB~SEX_LAB,
            scales="free_y",
            labeller=labeller(
                SEX_LAB=sex_dictionary2,
                LOC_LAB=community_dictionary$longest
        ) ) + 
        .scales() +
        my_labs(
            y="Viraemic individuals in age group as a proportion of the census-eligible population",
            color=my_labs_dictionary[colorby], fill=my_labs_dictionary[colorby], linetype=my_labs_dictionary["ROUND_LAB"]
        ) +
        nm_reqs
}
