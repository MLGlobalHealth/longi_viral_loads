inverse.logit <- function(x) 
    exp(x) / (exp(x) + 1)

simple.binomial.tests.by.community <- function(DT) {
    # DT <- copy(vlc)
    # COMM_NUM = 0 for Total across comunities
    by_cols <- c('ROUND', 'SEX')
    tmp0 <- DT[, .(NUM=N*PVLNS_MEAN, DEN=N*PHIV_MEAN), by=c(by_cols, 'COMM_NUM')] |> 
        cube(j=.(P=sum(NUM)/sum(DEN), NUM=sum(NUM),DEN=sum(DEN) ), by=c(by_cols, 'COMM_NUM')) |> 
        subset( !( is.na(ROUND) | is.na(SEX) ))
    tmp0 <- merge( 
        tmp0[! is.na(COMM_NUM)], 
        tmp0[is.na(COMM_NUM), .(P_ROUND_SEX= P) ,by=by_cols],
        by=by_cols
    )

    # now can we test whether the community empirical proportion is significantly lower.
    alpha_bonferroni <- function(alpha, n)
        1 - (1 -alpha)^(1/n)

    # let's try computing p-values with binom.test
    interval_test <- tmp0[COMM_NUM != 0, .( 
        P_ROUND_SEX = P_ROUND_SEX, 
        NUM=NUM,
        DEN=DEN,
        P_VALUE = binom.test(x=NUM, n=DEN, p=P_ROUND_SEX, alternative="two.sided")$p.value
    ), by=c(by_cols, 'COMM_NUM')] 
    interval_test[, REJECT:=P_VALUE <= alpha_bonferroni(.05, .N), by=by_cols ]
    interval_test[REJECT == TRUE]

    # now look at it 'longitudinally': combine p-values through Fisher's method
    interval_test[, P_COMB_FISHER:=poolr::fisher(P_VALUE)$p , by=c('COMM_NUM', 'SEX')] 
    if(0)
    {
        p <- interval_test |> 
            merge(dcomm[, .(COMM_NUM, TYPE)]) |> 
            subset(P<.05) |>
            ggplot(aes(x=ROUND, y=P_VALUE, linetype=TYPE,  color=as.factor(COMM_NUM)), size=P<.05 ) + 
                # geom_point() +
                facet_grid( SEX ~ .) +
                geom_line(aes(group=interaction(COMM_NUM, SEX))) +
                scale_size_manual(values=c(`TRUE`=5, `FALSE`=.05)) +
                scale_y_percentage +
                theme_default()
        

        vlc |> 
            subset( COMM_NUM %in% c(2, 8, 38, 74, 106, 370, 401, 602, 771, 773) ) |>
            ggplot( aes(x=PHIV_MEAN, y=PVLNSofHIV_MEAN)) +
            scale_x_continuous(labels=scales:::percent) +
            scale_y_continuous(labels=scales:::percent) +
            geom_hline(data=interval_test, aes(yintercept=P_ROUND_SEX)) +
            geom_errorbar(aes( ymin=PVLNSofHIV_CL, ymax=PVLNSofHIV_CU), alpha=0.2) +
            geom_errorbarh(aes( xmin=PHIV_CL, xmax=PHIV_CU), alpha=0.2) +
            geom_point(aes( colour=FC2)) +
            geom_text(aes( label=COMM_NUM), size=2) +
            facet_grid(ROUND~SEX) +
            # scale_colour_manual(values=palettes$comm) + 
            scale_colour_manual(values= palettes$comm2) +
            theme_default() + 
            my_labs(color='Community type')

        ## in both
        # 2 = generally good level of supp 
        # 8 = generally high
        # 38 = big fishing village, high sample size
        # 602 = nothing seems too extreme to me, but generally good?
        # 771 = big fishing village, high sample size
        # 773 = quite high in both men and women 
    }

    interval_test
}

plot.prev.viraemia.amonghiv.by.comm.round <- function(DT, MEANS=NULL)
{
    DT |> ggplot(aes(x = ROUND, y=PVLNSofHIV_MEAN, color=FC)) +
        geom_line(aes(group = COMM_NUM)) + {
        if(! is.null(MEANS))
            geom_line(
                data=MEANS,
                aes(x=ROUND, y=P_ROUND_SEX),
                color='black',
                linetype='dashed'
        ) 
        } +
        geom_text(aes(label = COMM_NUM), color='black', size = 4) +
        facet_grid(~SEX) +
        scale_y_continuous(labels = scales:::percent) +
        scale_colour_manual(values = palettes$comm) +
        theme_default() +
        my_labs( x = "Survey Round" )
    
}

get.glm.data <- function(DT)
{

    # get community data
    dcomm <- .get.dcomm()

    DT <- copy(dall)
    DT <- .preprocess.ds.oli(DT)
    DT[ COMM_NUM==22, COMM_NUM:=1 ]

    .f <- function(x,y)
        as.vector( unname ( binconf(sum(x), length(y) )  ) )

    vlc <- DT[, {
        z  <- .f( HIV_STATUS == 1, HIV_STATUS)
        z2 <- .f( VLNS == 1, VLNS )
        z3 <- .f( VLNS == 1, which(HIV_STATUS == 1) )
        list(FC=FC[1],
             N= length(HIV_STATUS),
             HIV_N = sum(HIV_STATUS == 1),
             VLNS_N = sum(VLNS == 1))		
    }, by=c('ROUND', 'COMM_NUM','SEX', 'AGEYRS')]

    # setkey(vlc, SEX, PHIV_MEAN)
    vlc[, SEX := factor(SEX, levels=c('M', 'F'), labels=c('men', 'women'))]
    vlc <- merge(vlc, dcomm[, .(COMM_NUM, FC2=TYPE) ] , by='COMM_NUM', all.x=TRUE)

    # additional columns
    comm_lvls <- vlc[, sort(unique(COMM_NUM)), by='FC2'][, V1, ]
    vlc[, SEX_LABEL := fifelse(SEX == 'women', yes='F', no='M')]
    vlc[, AGEGROUP := cut(AGEYRS, breaks=c(15, 24.01, 25, 34.01, 35, 50.01)-.01 ) ]
    .gs <- function(x) { y <- gsub(']|\\(', ' ', x);gsub(',', '_', y) }
    vlc[, AGEGROUP := .gs(AGEGROUP)]
    vlc
}

fit.glm.model <- function(formula, suffix, .outdir=glm.out.dir )
{
    # Loads or runs a stan model.

    stopifnot(is.character(suffix))
    filename <- file.path(.outdir, suffix )

    FileExists <- file.exists(filename)
    HasRandomEffects <- length(lme4::findbars(formula))>0
    NamesFormula <- unique(all.names(formula))
    NamesFormula <- NamesFormula[NamesFormula %like% '[A-z]']
    NamesFormulaExist <- all(NamesFormula %in% names(dglm))

    if( FileExists  )
    {
        cat('Loading fit...\n')
        fit <- readRDS(filename)

    }else if( NamesFormulaExist & HasRandomEffects ){

        fit <- stan_glmer(
            data=dglm,
            formula=formula,
            weights=HIV_N,
            prior_intercept = normal( prior.pars$intercept.mean, 5),
            family=binomial(link='logit'))
        saveRDS(fit, filename)


    }else if( NamesFormulaExist & ! HasRandomEffects ){

        fit <- stan_glm(
            data=dglm,
            formula=formula,
            weights=HIV_N,
            prior_intercept = normal( prior.pars$intercept.mean, 5),
            family=binomial(link='logit'))
        saveRDS(fit, filename)

    }else{
        stop('Model not found, and could not be run due to missing columns.\n')
    }

    return(fit)
}

analyse_glm_model <- function(glm, prefix)
{
    # TODO: prefix is not doing anythng atm
    # what need to be analysed:
    # - convergence
    # - pair plots 
    # - posterior predictive checks. 

    # glm <- copy(glm_5dec)
    summ <- summary(glm) |> as.data.table(keep.rownames = TRUE)

    # Get diagnostics and print
    summ[ ! rn %in% 'log-posterior', {
            z1 <- which.min(n_eff);
            z2 <- which.max(Rhat);
            list( NEFF=n_eff[z1], V_NEFF=rn[z1], RHAT=Rhat[z2], V_RHAT=rn[z2])
    } ] -> diagns

    cat('Minimum effective sample size:', diagns$NEFF, 'for parameter', diagns$V_NEFF, '\n')
    cat('Maximum Rhat:',diagns$RHAT, 'for parameter', diagns$V_RHAT, '\n')

    
    # get a table of the fit
    if ( 0 )
        modelsummary(glm, output = 'DT') # not optimal 4 models with multiple pars.

    # bayesplot::available_ppc()

    ## TOCLEAN 
    # POSTERIOR PREDICTIVE CHECKS
    # https://mc-stan.org/rstanarm/reference/plot.stanreg.html
    p_check <- pp_check(glm_model_choice)
    filename <- file.path('comm_stanglm_diagnostics_ppcheck_density.pdf')
    ggsave2(p_check, file=filename, w=9, h=8)

    pp_check(glm_model_choice, plotfun = "boxplot", nreps = 10, notch = FALSE)
    pp_check(glm_model_choice, plotfun = "scatter_avg") # y vs. average yrep
    pp_check(glm_model_choice, plotfun = "hist", nreps=3) # y vs. average yrep

    # Posterior vs prior
    p <- posterior_vs_prior( glm_model_choice) +
        guides(color=FALSE) + theme(axis.text.x = element_text(angle = 90))
    filename <- file.path( 'comm_stanglm_diagnostics_postvsprior.pdf')
    ggsave2(p, file=filename, w=9, h=6)

    # Statistical 'significance': probability of direction
    p_direction(glm_model_choice)
    library('sjPlot')
    sjPlot::tab_model(glm_model_choice)

    tmp <- modelsummary(glm_model_choice, statistic='mad', output='markdown')
    tmp
}

agresti.coull.exploration <- function(DT)
{
    # Get agresti-coull intervals (adding a ".25 prior")
    cols <- c("HIV_N", "VLNS_N", "N")
    by_cols <- c("FC2", "COMM_NUM", "AGEGROUP", "SEX_LABEL", "ROUND")

    DT[, lapply(.SD, sum), .SDcols = cols, by = by_cols ][,
        Hmisc::binconf(
            n = HIV_N + 4, x = VLNS_N + 1,
            method = "wilson", return.df = TRUE
    ), by = by_cols ] -> dvl_agresti

    dvl_agresti |> ggplot(aes(x = ROUND, y = logit(PointEst), color = FC2, linetype = FC2)) +
        geom_line(aes(group = COMM_NUM)) +
        facet_grid(AGEGROUP ~ SEX_LABEL) +
        labs(
            x = "Round",
            y = "Logit suppression",
            title = "Community-level suppression trajectories",
            subtitle = "(Agresti-coull intervals adding 1 VLNS of 4 HIV_pos)"
        ) -> p1
    filename <- "edaglm_logitprev_vs_round_by_agesexcomm.pdf"
    ggsave2(p1, file = filename, LALA = glm.out.dir, w = 9, h = 8)

    # aggregate by = something else
    dvl_agresti[,
        {
            z19 <- which(ROUND == "19")
            z16 <- which(ROUND == "16")
            list(Delta = PointEst[z19] - PointEst[z16])
        },
        by = c("FC2", "COMM_NUM", "AGEGROUP", "SEX_LABEL")
    ][,
        .(DIFF = median(Delta)),
        by = c("FC2", "SEX_LABEL")
    ] |> dcast(FC2 ~ SEX_LABEL, value.var = "DIFF")

    dvl_agresti |> ggplot(aes(x = as.factor(COMM_NUM), y = PointEst, color = as.factor(ROUND))) +
        geom_point() +
        facet_grid(SEX_LABEL ~ .) +
        theme(axis.text.x = element_text(angle = 90)) +
        labs(x = "Community Label", y = "Proportion of suppressed among study participants", color = "Round") -> p
    filename <- "edaglm_proprsupp_by_community_round_sex.pdf"
    ggsave2(p, file = filename, LALA = glm.out.dir, w = 9, h = 8)

    list(DT=dvl_agresti, p1=p1, p=p)
}

get.effective.sample.size.glm <- function(fit)
{
    tmp <- summary(fit) |> 
        as.data.frame() |>
        as.data.table(keep.rownames = TRUE)

    # exclude log-posterior
    tmp <- tmp[!rn %in% "log-posterior", {
        z1 <- which.min(n_eff)
        z2 <- which.max(Rhat)
        list(n = n_eff[z1], pn = rn[z1], R = Rhat[z2], pR = rn[z2])
    }]
    cat("Minimum effective sample size:", tmp$n, "for parameter", tmp$pn, "\n")
    cat("Maximum Rhat:", tmp$R, "for parameter", tmp$pR, "\n")
    tmp
}



run.bayesian.anova <- function(DT, formula)
{

    anova_formula <- VLNS_N / HIV_N ~ 
        (1 | FC) + 
        (1 | FC:COMM_NUM) + 
        (1 | SEX ) + 
        (1 | AGEGROUP ) + 
        (1 | SEX : AGEGROUP)

    fit <- glmer(
        data= dglm, 
        formula = ,
        weights = HIV_N,  
        family = binomial(link = "logit")
    )

    fit.sim <- arm::sim(fit)

    fit
    str(fit.sim) 
    fit.sim@fixef |> str()
    fit.sim@ranef |> str()

    fit.sim@ranef$FC |> str()
    fit.sim@ranef[['FC']][, , 1] 
    fit.sim@ranef$FC[[2]]
    fit.sim@ranef$FC[[3]]

}


make.gelman.anova.plto <- function(fit)
{
    ps <- c(CL=0.025, IL=0.25,M=0.5,IU=0.75,CU=0.975)

    ## NO, also need estimates of beta...
    #draws_sigma <- as.data.frame(glm_fullanova, regex_pars='Sigma')
    #sigma_quantiles <- draws_sigma |> 
    #    lapply(quantile, probs=ps ) |> 
    #    as.data.table() |> 
    #    t()
    #names(sigma_quantiles) <- names(ps)

    draws <- as.data.frame(glm_fullanova, regex_pars="\\(Intercept\\) ")

    # need to find each source of variation. 
    .find.source.variation.label <- function(lst)
    {
        lapply(lst, function(x){
            tmp <- stringr::str_split(x,":")[[1]]
            paste(tmp[1:(length(tmp)/2)], collapse = ':')
        }) |> unlist()
    }
    names(draws) <- gsub( "^b\\[\\(Intercept\\) |\\]$", "",names(draws) )  
    group_label <-  .find.source.variation.label(names(draws))

    # keeping same notation as gelman
    J_m <- table(group_label)
    # assuming the most basic constraitnts
    df_m <- J_m - 1 

    # Now I would need to implement the formula with the C_ms and use the mean/variance draws.
    
}

fit.all.glm.models <- function()
{
    # GLM on slopes
    glm_slopes <- fit.glm.model(
        suffix = "supppofpart_slopes_fit.rds",
        formula = VLNS_N / HIV_N ~ ROUNDi * (AGEGROUP + SEX + FC)
    )

    # GLM anova like, without random effects
    glm_anova <- fit.glm.model(
        formula = VLNS_N / HIV_N ~ ROUND * AGEGROUP * SEX * FC,
        suffix = "supppofpart_anova_fit.rds"
    )

    # without random effects (base)
    glm_base <- fit.glm.model(
        formula = VLNS_N / HIV_N ~ FC + ROUND:AGESEX + (1 | COMM_NUM),
        suffix = "suppofpart_base_fit.rds"
    )

    # with random effects by community
    glm_reffs <- fit.glm.model(
        formula = VLNS_N / HIV_N ~ FC + ROUND:AGESEX + (1 | COMM_NUM),
        suffix = "suppofpart_glmer_fit.rds"
    )

    # model from the 5th december
    glm_5dec <- fit.glm.model(
        formula = VLNS_N / HIV_N ~ SEX_LABEL + AGEGROUP + (1 | ROUND:FC) + (1 | COMM_NUM),
        suffix = "suppofpart_5thdec.rds"
    )

    # Take "most significant effects" from the anova model
    Glm3MostSignificant <- fit.glm.model(
        formula = VLNS_N / HIV_N ~ (1 | AGEGROUP:SEX_LABEL) + (1 | ROUND:FC) + (1 | COMM_NUM),
        suffix = "suppofpart_6mostsignificantanova.rds"
    )

    # Full Anova
    FullAnovaFormula <- VLNS_N / HIV_N ~
        (1 | SEX_LABEL) + (1 | ROUND) + (1 | AGEGROUP) + (1 | FC) +
        (1 | SEX_LABEL:ROUND) + (1 | SEX_LABEL:AGEGROUP) + (1 | SEX_LABEL:FC) + (1 | ROUND:AGEGROUP) + (1 | ROUND:FC) + (1 | AGEGROUP:FC) +
        (1 | SEX_LABEL:ROUND:AGEGROUP) + (1 | SEX_LABEL:ROUND:FC) + (1 | ROUND:AGEGROUP:FC) +
        (1 | SEX_LABEL:ROUND:AGEGROUP:FC)

    glm_fullanova <- fit.glm.model(
        formula = FullAnovaFormula,
        suffix = "suppofpart_fullanova.rds"
    )

    glm_8may <- fit.glm.model(
        formula = VLNS_N / HIV_N ~ (SEX_LABEL*ROUND) + AGEGROUP + FC + (1|COMM_NUM),
        suffix="suppofpart_sexround_age_fc_rancomm.rds"
    )
}

compare.random.effects <- function(fit)
{

    stopifnot(
        class(fit) %like% "glmerMod|stanreg|lmerMod" |> any()
    )
        
    ran_eff <- ranef(fit, whichel="COMM_NUM", condVar=FALSE) |> as.data.table()
    std_dev <- VarCorr(fit)$COMM_NUM |> as.numeric() |> sqrt()

    ran_eff |> as.data.table()
    p <- ggplot(ran_eff, aes(x=as.factor(grp), y=condval)) + 
        geom_point() +
        theme_default() + 
        geom_hline(aes(yintercept = +1.96*std_dev ), linetype='dotted') +
        geom_hline(aes(yintercept = -1.96*std_dev ), linetype='dotted') +
        labs(x='Community Number', y='Median Community-Level random effect') 
        NULL

    list(plot = p, ran_eff=ran_eff, std_dev=std_dev )
}
