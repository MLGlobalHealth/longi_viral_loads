process.hiv.testing.data <- function(DT, rounds=NULL, include.comm=FALSE)
{
    DTT <- copy(DT); 
    names(DTT) <- toupper(names(DTT)) 

    DT[, list( 
        STUDY_ID=STUDY_ID,
        ROUND=round2numeric(ROUND),
        COMM_NUM=COMM_NUM,
        TEST_EVER=fcase(
            RHIVEVER == 1, TRUE,
            RHIVEVER == 2, FALSE,
            default = NA),
        # assign to midpoint (and 'more than 4' == 5 )
        TEST_YEAR_AGO=fcase(
            HIVPERIOD == 1, .5,
            HIVPERIOD == 2, 1.5,
            HIVPERIOD == 3, 3,
            HIVPERIOD == 4, 5,
            default = NA),
        TEST_LAST_RESULT=fcase(
            HIVRSLT == 1, 'negative',
            HIVRSLT == 2, 'positive',
            default =NA
        )
    ), ] -> out

    if(!is.null(rounds))
    {
        stopifnot(rounds %in% unique(out$ROUND))
        out <- subset(out, ROUND %in% rounds)
    }
    if(! include.comm)
        out[, COMM_NUM := NULL]

    return(out)
}

make.testing.frequency.metric.among.negatives <- function(DT=dall2, include.new.diagnoses=FALSE)
{
    # for every community, average the time since last test across HIV negative
    # if no time since last test, set 10 years
    dnegatives <- subset(DT, HIV_STATUS == 0 & !is.na(TEST_EVER))

    # set NA test year ago as 10, don't know how meaningful this is
    dnegatives[is.na(TEST_YEAR_AGO), TEST_YEAR_AGO := 10]

    if(include.new.diagnoses)
    {
        stop('TODO')
    }

    out <- dnegatives[, list( 
        TEST_YA_MEAN=mean(TEST_YEAR_AGO)
        ) , by=c('COMM_NUM', 'SEX', 'ROUND')]

    return(out)
}

make.hivsuppression.metric.among.positives <- function(DT=dall2)
{
    # for every community, average the time since last test across HIV negative
    # if no time since last test, set 10 years
    dpositives <- subset(DT, HIV_AND_VL == 1)

    out <- dpositives[, list( 
        VL_SUPP_MEAN=mean(VL_COPIES <= VIREMIC_VIRAL_LOAD)
        ) , by=c('COMM_NUM', 'SEX', 'ROUND')]

    return(out)
}

plot.testing.answers.proportions <- function(
    DT=dall2,
    var,
    group='negatives',
    round.aggregate=TRUE,
    percentage=FALSE,
    chi2=FALSE)
{
    # DT=dall2; var='TEST_LAST_RESULT'; group='negatives'; round.aggregate=TRUE; percentage=FALSE
    # rm(DT, var, group, round.aggregate, percentage)

    stopifnot(var %in% names(DT))
    stopifnot(group %in% c('negatives', 'all', 'negatives_and_newdiagnoses', 'newdiagnoses'))
    by_cols <- c('ROUND', 'COMM_NUM','SEX', var)

    fill_lab <- fcase(
        var=="TEST_YEAR_AGO", "How long ago did you test?",
        var=="TEST_EVER", "Ever tested",
        var=="TEST_LAST_RESULT", "Result of last test",
        var=="TEST_LAST_YEAR", "Reported a test in the last year",
        default = "TOFILL!"
    )
    y_lab <- "N"

    dplot <- copy(DT)

    # select population denominator
    if(group == 'negatives')
    {
        dplot <- subset(dplot, HIV_STATUS==0)
        y_lab <- paste0(y_lab, " (among HIV negatives)") 
    }else if(group == 'negatives_and_newdiagnoses'){
        dplot <- subset(dplot, HIV_STATUS == 0 | NEW_DIAGN == TRUE)
        y_lab <- paste0(y_lab, '(among HIV negatives and new diagnoses)')
    }else if(group == 'newdiagnoses'){
        dplot <- subset(dplot, NEW_DIAGN == TRUE)
        y_lab <- paste0(y_lab, '(among new diagnoses only)')
    }else if(TRUE){
        warning('TODO')
    }

    dplot[, COMM_NUM := as.factor(COMM_NUM)]
    dplot[, (var) := lapply(.SD, as.factor) , .SDcols=var]
    dplot <- dplot[, .N ,by=by_cols]

    # do we need to aggregate over rounds?
    if(round.aggregate == TRUE)
    {
        aggregated <- dplot[, 
        list(ROUND='AGGR', N=sum(N)),
        by=setdiff(by_cols, 'ROUND')]
        dplot <- rbind(dplot, aggregated)
    }

    # compute proportions
    if(percentage)
        dplot[, P := N/sum(N), by=setdiff(by_cols, var)]

    if(chi2)
    {
        # are distributions different across groups?
        get.chi2.pvalue <- function(tab)
            summary(tab)$p.value
        pvalues <- dall2[ , .(P=table(COMM_NUM, get(var)) |> get.chi2.pvalue()), by=c('SEX', 'ROUND')]  
        pvalues[, P_LAB:=as.character(round(P,3))]
        pvalues[as.numeric(P_LAB) == '0', P_LAB := "<.001"]
        pvalues[, P_LAB := paste("chi2 p:", P_LAB)]
        pvalues[, `:=` (X=-Inf, Y=Inf, hjust='left', vjust='up')]
        pvalues |> prettify_labels()
    }

    prettify_var <- function(DT)
    {
        if(var %in% names(dall_dictionaries))
        {
            set(DT, NULL, var, as.character(DT[, get(var)]))
        }
        return(DT)
    }

    p <- dplot |> 
        prettify_labels() |>
        prettify_var() |>
        ggplot(aes(fill=NA)) +
        geom_col(aes_string(x="COMM_NUM", fill=var, y=fifelse(percentage, yes="P", no="N"))) + {
        if(chi2){
            geom_text(data=pvalues, aes(x=X, y=Y, label=P_LAB, vjust=1, hjust=0)) 
        }} +
        scale_color_discrete() +
        scale_y_continuous(limits=c(0,NA), expand=c(0,0, .1, 0))+
        facet_grid(SEX_LAB~ROUND_LAB, scales='free_y') +
        theme(legend.position='bottom',
              axis.text.x = element_text(angle = 45, hjust = 1, vjust=1),
        ) +
        guides(fill = guide_legend(override.aes = list(size = 1))) +
        labs(x="Community number", fill=fill_lab, y=y_lab) +
        NULL 

    if(var %in% c('TEST_EVER', 'TEST_LAST_RESULT', 'TEST_LAST_YEAR' ) )
    {
        p <- p + scale_fill_manual(
            values=palettes[['bool']],
        )
        p

    }else{
        p <- p + scale_fill_viridis_d(
            na.value = "grey80",
            labels=labs_from_dictionaries(dall_dictionaries[[var]])
        )
    }
}

get.avgsupp.avgtest <- function(DT, by_cols, group='negatives')
{
    stopifnot(group %in% c('negatives', 'negatives_and_newdiagnoses', 'newdiagnoses'))

    out <- DT[, {
        idx.posvl <- which(HIV_AND_VL == 1)

        idx.hivneg <- which(HIV_STATUS==0)
        idx.newdiagn <- which(NEW_DIAGN==TRUE)
        idx2 <- c()
        if(group %like% 'negatives')
            idx2 <- c(idx2, idx.hivneg)
        if(group %like% 'newdiagnoses')
            idx2 <- c(idx2, idx.newdiagn)
        idx2 <- unique(idx2)

        list( 
            GROUP=group,
            PROP_VLSUPP=mean(VL_COPIES[idx.posvl] <= VIREMIC_VIRAL_LOAD),
            PROP_TEST=mean(TEST_EVER[idx2], na.rm=TRUE),
            PROP_TEST_LESS1=mean(TEST_YEAR_AGO[idx2]<=1, na.rm=TRUE),
            PROP_TEST_LESS2=mean(TEST_YEAR_AGO[idx2]<=2, na.rm=TRUE),
            PROP_TEST_MORE5=mean(TEST_YEAR_AGO[idx2]>=4.5 | is.na(TEST_YEAR_AGO), na.rm=TRUE)
        )}, by=by_cols ] 

    if(! 'ROUND' %in% by_cols)
        out[, ROUND := "AGGR"]

    return(out)
}

plot.suppression.vs.testing.answers.proportions <- function(DT,.group, type)
{
    stopifnot(type %in% c('evertest','test1ya', 'test5ya'))

    .xlab <-  gsub(pattern="_and_", replacement=' and ', .group) |>
        gsub( pattern="negatives", replacement="HIV negatives") |> 
        gsub( pattern="newdiagnoses", replacement="new diagnoses") 

    dplot <- get.avgsupp.avgtest(DT, by_cols = c('COMM_NUM','FC','ROUND', 'SEX'), group=.group)
    dplot_agg <- get.avgsupp.avgtest(DT, by_cols = c('COMM_NUM','FC','SEX'), group=.group)
    dplot <- rbind(dplot, dplot_agg) |> prettify_labels()

    if(type=='evertest')
    {
        p <- ggplot(dplot[ROUND == 'AGGR'], aes(y=PROP_VLSUPP, color=SEX_LAB)) +
            geom_point(aes(x=PROP_TEST)) +
            scale_x_continuous(limits=c(0,1)) +
            scale_y_continuous(limits=c(0,1)) +
            scale_color_manual(values=palettes[['sex']]) +
            theme(legend.position='bottom') +
            labs(
                x=paste0("Proportion of ",.xlab," reporting prior tests"),
                y="Proportion of HIV positives with suppressed viral load",
                color="Gender")
        return(p)
    }else if(type %like% 'test[0-9]ya'){

        var <- fcase(
            type == 'test1ya', "PROP_TEST_LESS1",
            type == 'test5ya', "PROP_TEST_MORE5"
        )

        .xlab2 <- fcase(
            type == 'test1ya', 'reporting test in the last year',
            type == 'test5ya', 'reporting no test or test more than 5 years ago'
        )

        p <- dplot |> 
            prettify_labels() |> 
            ggplot(aes_string(y='PROP_VLSUPP',x=var, color='FC')) +
            geom_text(aes(label=COMM_NUM), size=2) +
            ggpubr::stat_cor() +
            facet_grid(SEX_LAB~ROUND_LAB) +
            scale_x_continuous(limits=c(0,1)) +
            scale_y_continuous(limits=c(0,1)) + 
            scale_color_manual(values=palettes[['comm']]) +
            theme(legend.position='bottom')+
            labs(
                x=paste('Proportion of',.xlab, .xlab2),
                y='Proportion of Virally suppressed among positives',
                color="Community type")
        return(p)
    }
}
