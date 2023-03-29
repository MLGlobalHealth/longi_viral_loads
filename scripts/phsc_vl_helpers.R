get.dall <- function(path, make_flowchart=TRUE)
{
    # Load data: exclude round 20 as incomplete
    dall <- fread(path)
    dall <- unique(dall[ROUND >= 16 & ROUND <= 19])

    # remove useless cols:
    dall[, CURR_ID := NULL]
    dall <- unique(dall)

    # study dataset keys
    key_cols <- c('STUDY_ID', 'ROUND')
    setkeyv(dall, key_cols)
    stopifnot(dall[, .N, by=key_cols][, all(N == 1)])

    # rename variables according to Oli's old script + remove 1 unknown sex
    setnames(dall,
             c('HIV_VL', 'COMM'),
             c('VL_COPIES', 'FC') )
    dall[, HIV_AND_VL := ifelse( HIV_STATUS == 1 & !is.na(VL_COPIES), 1, 0)]

    # remove individuals without sex reported
    dall <- dall[! SEX=='']

    # make study flowchart if packages are installed 
    if('DiagrammeR' %in% installed.packages() & make_flowchart)
        make.study.flowchart(dall)

    return(dall)
}

.get.dcomm <- function()
{
    path <- file.path(indir.deepanalyses.xiaoyue, 'PANGEA2_RCCS1519_UVRI', 'community_names.csv')
    dcomm <- fread(path)
    .f <- function(x)
        fcase(x %like% '^f', 'fishing', x %like% '^t', 'trading', x %like% '^a', 'agrarian')
    dcomm[, TYPE:=.f(COMM_NUM_A)]
    dcomm[, COMM_NUM:=as.integer(COMM_NUM_RAW)]
    dcomm
}

make.study.flowchart <- function(DT)
{
    
    library(DiagrammeR)
    library(DiagrammeRsvg)
    library(rsvg)
    library(htmltools)
    ls_fchart <- list()

    # Get numbers for flowchart
    ls_fchart[['all']] <- DT[, .(uniqueN(STUDY_ID), .N)]
    ls_fchart[['pos']] <- DT[HIV_STATUS == 1, .(uniqueN(STUDY_ID), .N)]
    ls_fchart[['vls']] <- DT[HIV_AND_VL == 1, .(uniqueN(STUDY_ID), .N)]
    ls_fchart <- lapply(ls_fchart, function(DT) DT[, lapply(.SD, formatC, big.mark=',')])

    idx <- DT[, unique(ROUND)]
    statement_for_round <- function(r)
    {
      .statement <- function(tot,x, y, z)
      {
        .r <- function(a) 
          paste0(formatC(a, big.mark = ','), ' (',format(round(100*a/tot, 2), nsmall=2), '%)')
        paste0(
          "Round ",r, "\n", 
          formatC(tot, big.mark = ','), " total VL measurements:\\\\l &#8226; ",
          .r(x), " were vireamic.\\\\l &#8226; ",
          .r(y), " had low level viraemia.\\\\l &#8226; ",
          .r(z), " were non-viraemic.\\\\l")
        
      }
      DT[HIV_AND_VL == 1 & ROUND == r, 
         .statement(tot=.N,
                    x=sum(VL_COPIES > 1000),
                    y=sum(VL_COPIES %between% c(200, 1000)),
                    z=sum(VL_COPIES < 1000)) ]
    }
    ls_fchart_round <- lapply(idx, statement_for_round)
    names(ls_fchart_round) <- paste0('r',idx)

    filename <- file.path(out.dir, 'flowchart_numbers.pdf')
    DiagrammeR::grViz(
                      "
                      digraph graph2 {

                          # Graph attributes
                          graph [splines=line, nodesep=.2]
                          node [fontsize=15]
                          edge [arrowsize=.5]

                          # node definitions with substituted label text
                          node [shape = rectangle, width = 7, height=.4]
                          a [label = '@@1']
                          b [label = '@@2']
                          c [label = '@@3']

                          node [shape = rectangle, width = 1]
                          r16 [label='@@4']
                          r17 [label='@@5']
                          r18 [label='@@6']
                          r19 [label='@@7']

                          node [shape = point, width=.001, height=.001, label='']
                          h; h1; h2; h3; h4

                          # Edge statement
                          a -> b -> c
                          edge [arrowhead=none]
                          c -> h 
                          {rank=same; arrowhead=none; h1 -> h2 ; h2 -> h; h -> h3; h3 -> h4}

                          edge[ tailclip = true, headclip=true, arrowhead=normal]
                          h1 -> r16
                          h2 -> r17
                          h3 -> r18
                          h4 -> r19
                      }

                      [1]: paste0(ls_fchart[['all']][[1]], ' study participants, accounting for ', ls_fchart[['all']][[2]], ' visits.')
                      [2]: paste0(ls_fchart[['pos']][[1]], ' HIV positive participants, accounting for ', ls_fchart[['pos']][[2]], ' visits.')
                      [3]: paste0(ls_fchart[['vls']][[1]], ' HIV positive participants with VL measurements,\\naccounting for ', ls_fchart[['vls']][[2]], ' visits.')
                      [4]: ls_fchart_round$r16
                      [5]: ls_fchart_round$r17
                      [6]: ls_fchart_round$r18
                      [7]: ls_fchart_round$r19
                      ")  |> export_svg() |> charToRaw() |> rsvg_pdf(filename)
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
    # what need to be analysed:
    # - convergence
    # - pair plots 
    # - posterior predictive checks. 

    glm <- copy(glm_5dec)
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
    filename <- file.path('comm_stanglm_diagnostics_ppcheck_density.png')
    ggsave2(p_check, file=filename, w=9, h=8)

    pp_check(glm_model_choice, plotfun = "boxplot", nreps = 10, notch = FALSE)
    pp_check(glm_model_choice, plotfun = "scatter_avg") # y vs. average yrep
    pp_check(glm_model_choice, plotfun = "hist", nreps=3) # y vs. average yrep

    # Posterior vs prior
    p <- posterior_vs_prior( glm_model_choice) +
        guides(color=FALSE) + theme(axis.text.x = element_text(angle = 90))
    filename <- file.path( 'comm_stanglm_diagnostics_postvsprior.png')
    ggsave2(p, file=filename, w=9, h=6)

    # Statistical 'significance': probability of direction
    p_direction(glm_model_choice)
    library('sjPlot')
    sjPlot::tab_model(glm_model_choice)

    tmp <- modelsummary(glm_model_choice, statistic='mad', output='markdown')
    tmp


}

.get.prior.ranges <- function(stan.data, DT, shape, scale)
{
    # Extracts alpha and rho from stan.data
    # Then computes 95% intervals assuming inv.gamma and normal prior distributions
    cols <- grep('alpha_hyper', names(stan.data), value=TRUE)
    tmp <- unlist(as.vector(stan.data[cols]))
    tmp <- data.table(par = names(tmp), val = as.vector(tmp))

    .f <- function(reg,x)
        as.integer(gsub('^.*?_([0-9])([0-9])$', reg, x))

    tmp[, `:=`(
               SEX = .f('\\2', par),
               LOC = .f('\\1', par),
               GP_hyper_par = gsub('^.*?(alpha|rho).*?$','\\1', par)
               )]

    tmp <- rbind(
                 tmp,
                 data.table(par=c('rho_hyper_scale', 'rho_hyper_shape'),
                            val=c(scale, shape), SEX=NA,
                            LOC=NA, GP_hyper_par='rho')
    )
    tmp

    ps <- c(0.025,0.25,0.5,0.75,0.975)
    cols <- c('CL','IL','M','IU','CU')

    # tmp[GP_hyper_par=='rho', (cols) := as.list(invgamma::qinvgamma(p=ps, shape=val, scale=val)), by=val]

    qhalfnorm <- function(p, sd)
    {
        x <- abs(rnorm(1000000, mean=0, sd=sd))
        as.list(quantile(x,ps))
    }

    qinvgamma <- function(p, sh, sc)
    {
        x <- MCMCpack::rinvgamma(1000000, shape=sh, scale=sc)
        #CHECK ps2 <- c(0.01, 0.99)
        as.list(quantile(x,ps))
    }

    tmp[GP_hyper_par=='alpha', (cols) := qhalfnorm(p=ps, sd=val), by=val]
    tmp[GP_hyper_par=='rho', (cols) := qinvgamma(p=ps, sh=shape, sc=scale) ]

    tmp1 <- tmp[GP_hyper_par=='rho', list(
                                          GP_hyper_par='rho',
                                          SEX=c(0,0,1,1),
                                          LOC=c(0,1,0,1),
                                          CL=CL[1],
                                          CU=CU[1])]

    cols <- names(tmp1)
    tmp <- rbind(tmp1, tmp[GP_hyper_par=='alpha', ..cols])

    tmp1 <- unique(DT[, .(SEX,SEX_LABEL,LOC,LOC_LABEL)])
    tmp <- merge(tmp1, tmp, by=c('SEX','LOC'))
    tmp[, col:='prior']
    tmp
}

.make.stan.data.gp <- function(DTsd, 
                            num.var=NA,
                            den.var=NA,
                            rho_hyper_lower_bound=1, rho_hyper_upper_bound=35/2)
{
        stopifnot(length(num.var==1) & length(den.var==1))

        tmp <- seq(DTsd[, min(AGE_LABEL)], DTsd[, max(AGE_LABEL)+1], .5)
        tmp1 <- which(tmp%%1==0.5)

        # num.var='HIV_N'
        # den.var='N'
        stan.data <- list(
                x_predict = tmp,
                N_predict = length(tmp),
                observed_idx = tmp1,
                N_observed = length(tmp1),
                y_observed_00 = DTsd[SEX==0 & LOC==0, ..num.var][[1]],
                y_observed_10 = DTsd[SEX==1 & LOC==0, ..num.var][[1]],
                y_observed_01 = DTsd[SEX==0 & LOC==1, ..num.var][[1]],
                y_observed_11 = DTsd[SEX==1 & LOC==1, ..num.var][[1]],
                total_observed_00 = DTsd[SEX==0 & LOC==0, ..den.var][[1]],
                total_observed_10 = DTsd[SEX==1 & LOC==0, ..den.var][[1]],
                total_observed_01 = DTsd[SEX==0 & LOC==1, ..den.var][[1]],
                total_observed_11 = DTsd[SEX==1 & LOC==1, ..den.var][[1]],
                alpha_hyper_par_00 = 2,
                alpha_hyper_par_10 = 2,
                alpha_hyper_par_01 = 2,
                alpha_hyper_par_11 = 2,
                rho_hyper_lower_bound = rho_hyper_lower_bound,
                rho_hyper_upper_bound = rho_hyper_upper_bound 
        )

        stan.data
}

.plot.gp.hyperparameters <- function(GP, PR)
{
        p <- ggplot(GP, aes(colour=col,
                            x= paste0(GP_hyper_par, ' ', LOC_LABEL, ' ', SEX_LABEL)
                            )) + 
                geom_errorbar(data=PR, aes(ymin=CL, ymax=CU))  +
                geom_linerange(aes(ymin=CL, ymax=CU)) +
                geom_point(aes(y=M)) +
                coord_flip() +
                theme_bw() +
                theme(legend.position='bottom') +
                labs(x='GP hyperparameter\n', y='', colour='',
                     title='95% Credible Intervals')
        p
}

.write.CIs <- function(m, l, u, d=1, percent=F)
{
        .r <- function(x) round(x*multiplier, digits=d)
        if(percent)
        {
                multiplier <- 100
                out <- paste0( .r(m),'% (', .r(l),'% - ', .r(u) ,'%)' )
        }else{
                multiplier <- 1
                paste0( .r(m),' [', .r(l),'-', .r(u) ,']' )
        }
}

date2numeric <- function( x )
{
	if(!class(x)%in%c('Date','character'))	return( x )
	x	<- as.POSIXlt(x)
	tmp	<- x$year + 1900
	x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
	x	
}

.preprocess.ds.oli <- function(DT, rm.na.vl=TRUE)
{
        # DT <- copy(dall)
        DT <- subset(DT, AGEYRS <= 50)

        # remove HIV+ individuals with missing VLs and 
        if(rm.na.vl)
                DT <- subset(DT, HIV_STATUS==0 | HIV_AND_VL==1)

        # consider only ARVMED for infected
        set(DT, DT[, which(ARVMED==1 & HIV_STATUS==0)], 'ARVMED', 0) 
        
        # define VL_COPIES for uninfected
        set(DT, NULL, 'VLC', DT$VL_COPIES)
        set(DT, DT[,which(HIV_STATUS==0)], 'VLC', 0)
        
        # define undetectable VL (machine-undetectable)
        # define suppressed VL (according to WHO criteria)	
        set(DT, NULL, 'VLU', DT[, as.integer(VLC<VL_DETECTABLE)])
        set(DT, NULL, 'VLS', DT[, as.integer(VLC<VIREMIC_VIRAL_LOAD)])
        set(DT, NULL, 'VLD', DT[, as.integer(VLC>=VL_DETECTABLE)])
        set(DT, NULL, 'VLNS', DT[, as.integer(VLC>=VIREMIC_VIRAL_LOAD)])
        set(DT, NULL, 'HIV_AND_VLD', DT[, as.integer(VLD==1 & HIV_AND_VL==1)])
        
        # reset VLC below machine detectable to 0
        set(DT, DT[, which(HIV_AND_VL==1 & VLU==1)], 'VLC', 0)
        setkey(DT, ROUND, FC, SEX, AGEYRS)

        DT
}

.preprocess.make.vla <- function(DT, select=c('N', 'HIV_N', 'VLNS_N', 'ARV_N', 'VL_MEAN', 'VL_SD', 'VL_MEAN_SD', 'HIV_N'))
{
	tmp <- seq.int(min(DT$AGEYRS), max(DT$AGEYRS))
        tmp1 <- DT[, sort(unique(ROUND))]

	vla <- as.data.table(expand.grid(ROUND=tmp1,
                                         FC=c('fishing','inland'),
                                         SEX=c('M','F'),
                                         AGEYRS=tmp))
	vla <- vla[, {		
                z <- which(DT$ROUND==ROUND & DT$FC==FC & DT$SEX==SEX & DT$AGEYRS==AGEYRS)	
                list(N          = length(z),
                     HIV_N      = sum(DT$HIV_STATUS[z]==1),
                     VLNS_N     = sum(DT$VLNS[z]==1),
                     ARV_N      = sum(DT$ARVMED[z]==0 & DT$HIV_STATUS[z]==1 & !is.na(DT$ARVMED[z])),
                     VL_MEAN    = mean(DT$VLC[z]),
                     VL_SD      = sd(DT$VLC[z]),
                     VL_MEAN_SD = sd(DT$VLC[z]) / sqrt(length(z)),
                     HIV_N      = sum(DT$HIV_STATUS[z]==1)
                )
        }, by=names(vla)]

	setnames(vla, c('FC','SEX','AGEYRS'), c('LOC_LABEL','SEX_LABEL','AGE_LABEL'))
	vla[, LOC:= as.integer(LOC_LABEL=='fishing')]
	vla[, SEX:= as.integer(SEX_LABEL=='M')]
	vla[, AGE:= AGE_LABEL-14L]
	vla[, ROW_ID:= seq_len(nrow(vla))]

        # Extract selected fields
        cols <- c('ROUND', 'LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL', 'LOC', 'SEX', 'AGE', 'ROW_ID')
        cols <- unique(c(cols, select))
        vla[, ..cols]
}

.preprocess.make.vla.2 <- function(DT, select=c('N', 'PHIV_MEAN', 'PHIV_CL', 'PHIV_CU', 'PVLNS_MEAN', 'PVLNS_CL', 'PVLNS_CU', 'PVLNSofHIV_MEAN', 'PVLNSofHIV_CL', 'PVLNSofHIV_CU'))
{
	tmp <- seq.int(min(DT$AGEYRS), max(DT$AGEYRS))
        tmp1 <- DT[, sort(unique(ROUND))]
	vla <- as.data.table(expand.grid(ROUND=tmp1, 
                                         FC=c('fishing','inland'),
                                         SEX=c('M','F'),
                                         AGEYRS=tmp))

        .f <- function(x,y) as.vector( unname ( binconf( sum(x), length(y) )))
	ans <- vla[, {		
				z <- which(DT$ROUND == ROUND & DT$FC==FC & DT$SEX==SEX & DT$AGEYRS<=(AGEYRS+2) & DT$AGEYRS>=(AGEYRS-1))				
                                z2 <- .f(DT$HIV_STATUS[z]==1, z)
                                z3 <- .f(DT$VLNS[z]==1, z)
                                z4 <- .f(DT$VLNS[z]==1, which(DT$HIV_STATUS[z]==1))
                                # z5 <- .f(DT$ARVMED[z]==0 & DT$HIV_STATUS[z] & !is.na(DT$ARVMED), which(DT$hiv_status[z]==1 & !is.na(DT$arvmed[z]))) 

				list(
                                     N= length(z),
                                     PHIV_MEAN= z2[1],
                                     PHIV_CL= z2[2],
                                     PHIV_CU= z2[3],				 
                                     PVLNS_MEAN= z3[1],
                                     PVLNS_CL= z3[2],
                                     PVLNS_CU= z3[3],
                                     PVLNSofHIV_MEAN= z4[1],
                                     PVLNSofHIV_CL= z4[2],
                                     PVLNSofHIV_CU= z4[3]#,
                                     # PARVofHIV_MEAN= z5[1],
                                     # PARVofHIV_CL= z5[2],
                                     # PARVofHIV_CU= z5[3]
                                )				


        }, by=names(vla)]

	set(ans, NULL, 'SEX', ans[, factor(SEX, levels=c('M','F'), labels=c('men','women'))])

        cols <- c(names(vla), select)
        return(ans[, ..cols])
}

vl.get.eligible.round17<- function()
{
    require(data.table)
    
    infile <- "~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/rakai_elibility.rda"
    outfile.base <- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/"
    load(infile)
    
    # subset to data of interest
    de <- as.data.table(eldat)
    de <- subset(de, status%in%c('_Participated','Away','Blood refusal','Missing data','Other','Refused','urine sample'))
    de <- subset(de, visit==17)

}

vl.vlprops.by.comm.gender.loc<- function(DT, write.csv=FALSE)
{
    # DT <- copy(dall)
    DT <- .preprocess.ds.oli(DT)
    # merge two communities that fully overlap, so we have 40 communities in the end 
    DT[ COMM_NUM==22, COMM_NUM:=1 ]

    # get.community types
    dcomm <- .get.dcomm()

    # calculate HIV prevalence and proportion not suppressed of HIV+ by community and gender
    .f <- function(x,y)
        as.vector( unname ( binconf(sum(x), length(y) )  ) )

    vlc <- DT[, {
        z  <- .f( HIV_STATUS == 1, HIV_STATUS)
        z2 <- .f( VLNS == 1, VLNS )
        z3 <- .f( VLNS == 1, which(HIV_STATUS == 1) )
        list(FC=FC[1],
             N= length(HIV_STATUS),
             PHIV_MEAN= z[1],
             PHIV_CL= z[2],
             PHIV_CU= z[3],
             PVLNS_MEAN= z2[1],
             PVLNS_CL= z2[2],
             PVLNS_CU= z2[3],
             PVLNSofHIV_MEAN= z3[1],
             PVLNSofHIV_CL= z3[2],
             PVLNSofHIV_CU= z3[3],
             VLC_MEAN= mean(VLC))
    }, by=c('ROUND', 'COMM_NUM','SEX')]

    .f <- function(m, l, u)
    {
        .r <- function(x) round(x*100, digits=1)
        paste0( .r(m),' [', .r(l),'-', .r(u) ,']' )
    }

    vlc[, PHIV_L:= .f(PHIV_MEAN, PHIV_CL, PHIV_CU)]
    vlc[, PVLNS_L:= .f(PVLNS_MEAN, PVLNS_CL, PVLNS_CU)]
    vlc[, PVLNSofHIV_L:= .f(PVLNSofHIV_MEAN, PVLNSofHIV_CL, PVLNSofHIV_CU)]

    setkey(vlc, SEX, PHIV_MEAN)
    vlc[, SEX := factor(SEX, levels=c('M', 'F'), labels=c('men', 'women'))]
    names(dcomm)
    vlc <- merge(vlc, dcomm[, .(COMM_NUM, FC2=TYPE) ] , by='COMM_NUM', all.x=TRUE)


    p_inf <- ggplot(vlc) +
        scale_x_continuous(labels=scales:::percent) +
        scale_y_continuous(labels=scales:::percent) +
        geom_errorbar(aes(x=PHIV_MEAN, ymin=PVLNSofHIV_CL, ymax=PVLNSofHIV_CU), alpha=0.2) +
        geom_errorbarh(aes(y=PVLNSofHIV_MEAN, xmin=PHIV_CL, xmax=PHIV_CU), alpha=0.2) +
        geom_point(aes(x=PHIV_MEAN, y=PVLNSofHIV_MEAN, colour=FC2)) +
        geom_text(aes(x=PHIV_MEAN, y=PVLNSofHIV_MEAN, label=COMM_NUM), size=2) +
        facet_grid(ROUND~SEX) +
        # scale_colour_manual(values=palettes$comm) + 
        scale_colour_manual(values= palettes$comm2) +
        theme_bw() +
        theme(legend.position='bottom') + 
        labs(x='\nHIV prevalence', 
             y='proportion unsuppressed HIV among infected\n', 
             colour='community type')
        p_inf


    filename <- file.path('220729_hivnotsuppofhiv_vs_hivprev_by_round_gender_fishinland.pdf')
    ggsave2(p_inf, file=filename, LALA=glm.out.dir, w=9, h=12)


    p_pop <- ggplot(vlc) +
        scale_x_continuous(labels=scales:::percent) +
        scale_y_continuous(labels=scales:::percent) +
        geom_errorbar(aes(x=PHIV_MEAN, ymin=PVLNS_CL, ymax=PVLNS_CU), alpha=0.2) +
        geom_errorbarh(aes(y=PVLNS_MEAN, xmin=PHIV_CL, xmax=PHIV_CU), alpha=0.2) +
        geom_point(aes(x=PHIV_MEAN, y=PVLNS_MEAN, colour=FC)) +
        geom_text(aes(x=PHIV_MEAN, y=PVLNS_MEAN, label=COMM_NUM), size=2) +
        facet_grid(ROUND~SEX) +
        scale_colour_manual(values=palettes$comm2) + 
        theme_bw() +
        theme(legend.position='bottom') + 
        labs(x='\nHIV prevalence', 
             y='proportion unsuppressed HIV among population\n', 
             colour='community type')

    filename <- file.path('220729_hivnotsuppofpop_vs_hivprev_by_round_gender_fishinland.pdf')
    ggsave2(p_pop, file=filename,LALA=glm.out.dir, w=9, h=12)

    if(write.csv)
    {
        # write results to file
        filename <- file.path(glm.out.dir,
                              '220729_hivnotsuppofhiv_vs_hivprev_by_round_gender_fishinland.csv')
        fwrite(vlc, file=filename)
    }

    list(DT=vlc, p_among_pop=p_pop, p_among_inf=p_inf)
}

vl.prevalence.by.gender.loc.age.gp <- function(DT, refit=FALSE, vl.out.dir.=vl.out.dir) 
{
    cat('\n\n--- Analysing HIV+ Prevalence ---\n\n')

    # DT <- copy(dall); refit=FALSE
    DT <- .preprocess.ds.oli(DT)

    tmp <- c('N', 'HIV_N', 'VLNS_N', 'ARV_N')
    vla <- .preprocess.make.vla(DT, select=tmp)

    # Stan file locations
    file.stan <- file.path(path.stan, 'vl_binomial_gp.stan')
    
    .fit.stan.and.plot.by.round <- function(DT, itr=10e3, wrmp=5e2, chns=1, cntrl=list(max_treedepth= 15, adapt_delta= 0.999))
    {
        #  DT <- copy(vla[ROUND == 16]); vl.out.dir.=vl.out.dir
        round <- DT[, unique(ROUND)]
        stopifnot(length(round) == 1)

        cat('Fitting stan model for round ', round, '\n')
        stan.data <- .make.stan.data.gp(DTsd=DT, num.var='HIV_N', den.var='N')
        filename <- file.path(vl.out.dir., 
                              paste0("hivprevalence_gp_stanfit_round",round,"_220729.rds"))

        if(file.exists(filename) & refit == FALSE)
        {
            cat('Loading previously run HMC... \n')
            fit <- readRDS(filename)
        }else{
            stan.model <- stan_model(file=file.stan, model_name= 'gp_all')
            fit <- sampling(stan.model, data=stan.data, 
                            iter=itr, warmup=wrmp,
                            chains=chns, control=cntrl)
            saveRDS(fit, file=filename)
        }

        # Analyse posterior
        # _________________

        # Extract summary excluding hyper parameters
        dsum <- summary(fit)$summary
        idx <- rownames(summary(fit)$summary) %like% 'rho_hyper_par'
        dsum <- dsum[!idx,]
        
        cat('The minimum effective sample size is:\n')
        cat(min( dsum[, 'n_eff'], na.rm=T ), '\n')
        cat('Out ot ', nrow(dsum), 'parameters\n')
        cat(sum(dsum[, 'n_eff' ] < 4000, na.rm=TRUE), ' had n_eff < 4000,\n')
        cat(sum(dsum[, 'n_eff' ] < 1000, na.rm=TRUE), ' had n_eff < 1000.\n')

        re <- rstan::extract(fit)
        ps <- c(0.025,0.25,0.5,0.75,0.975)

        # extract hyperparams rho		
        tmp <- grep('rho_[0-9]|alpha_[0-9]', names(re), value=TRUE)
        .f <- function(x) transpose(as.data.table(quantile(x, probs=ps)))
        
        tmp <- lapply(re[tmp], .f)
        tmp <- rbindlist(tmp, idcol='par')
        names(tmp) <- c('GP_hyper_par','CL','IL','M','IU','CU')

        .f <- function(reg, x)
                gsub('^([a-z]+)_([0-9])([0-9])',reg,x)
        tmp <- tmp[ ! GP_hyper_par %like% 'bound']
        tmp[, SEX := as.integer(.f('\\2',GP_hyper_par))]
        tmp[, LOC := as.integer(.f('\\3',GP_hyper_par))]
        tmp[, GP_hyper_par := .f('\\1',GP_hyper_par)]

        tmp1 <- unique(DT[, .(SEX,SEX_LABEL,LOC,LOC_LABEL)])
        prev.hiv.gp.pars <- merge(tmp1, tmp, by=c('SEX','LOC'))
        prev.hiv.gp.pars[, col:='posterior']
        prev.hiv.gp.pars[, GP_hyper_par:=gsub('_[0-1][0-1]', '',GP_hyper_par)]

        tmp <- .get.prior.ranges(stan.data, DT=tmp1, 
                shape=re$rho_hyper_par_shape2[1],
                scale=re$rho_hyper_par_scale2[1])

        p <- .plot.gp.hyperparameters(prev.hiv.gp.pars, tmp)
        filename=paste0('220729_hivprevalence_gppars_round',round,'.pdf')
        ggsave2(p, file=filename, w=6, h=3)

        # make use stan.data for PPC.
        # ___________________________

        ppDT <- copy(DT)
        cols <- c('M', 'CU', 'CL')
        ppDT[, (cols) := binconf(HIV_N, N, return.df=T)]

        # make prevalence plot by age
        # ___________________________

        .f <- function(x) apply(x, 2, quantile, probs=ps)
        tmp <- cbind( .f(re$p_predict_00), .f(re$p_predict_10),
                      .f(re$p_predict_01), .f(re$p_predict_11))
        rownames(tmp) <- c('CL','IL','M','IU','CU')
        tmp <- as.data.table(reshape2::melt(tmp))
        prev.hiv.by.age <- dcast.data.table(tmp, Var2~Var1, value.var='value')
        tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
        prev.hiv.by.age <- cbind(tmp, prev.hiv.by.age) 
        prev.hiv.by.age <- merge(unique(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), prev.hiv.by.age, by=c('SEX','LOC'))

        p <- ggplot(prev.hiv.by.age, aes(x=AGE_LABEL, colour=SEX_LABEL)) + 		
                geom_ribbon(aes( ymin=CL, ymax=CU, fill=SEX_LABEL), alpha=0.2, col=NA) +
                geom_line(aes( y=M)) +
                scale_x_continuous( expand=c(0,0) ) + 
                scale_y_continuous(labels=scales:::percent, expand=c(0,0), limits=c(0,.75)) +
                scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
                scale_fill_manual(values=c('M'='royalblue3','F'='deeppink2')) +
                facet_wrap(~LOC_LABEL, ncol=2) +
                theme_bw() +
                theme(legend.position='bottom') +
                labs(x='\nage at visit (years)', 
                     y='HIV prevalence (95% credibility interval)\n', 
                     colour='gender', fill='gender')

        filename=paste0('220729_hivprevalence_vs_age_by_gender_fishinland_gp_round',round,'.pdf')
        ggsave2(p, file=filename, w=6, h=5)

        p <- p  + 
                geom_point(data=ppDT, aes(y=M, pch=SEX_LABEL, size=N)) +
                # geom_linerange(data=ppDT, aes(ymin=CL, ymax=CU), linetype='dotted' ) +
                scale_size(range = c(0, 3)) +
                labs(pch='gender', linetype='gender', size='population size')

        filename=paste0('220729_hivprevalence_vs_age_by_gender_fishinland_data_gp_round',round,'.pdf')
        ggsave2(p, file=filename, w=6, h=5)


        # extract basic prevalence estimates
        # __________________________________

        ps <- c(0.025, 0.5, 0.975)
        tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
        rp <- as.data.table(reshape2::melt(tmp))
        tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict-0.5, SEX=c(0,1), LOC=c(0,1)))
        tmp[, Var2:= seq_len(nrow(tmp))]
        rp <- merge(tmp, rp, by='Var2')
        setnames(rp, c('Var1','value'), c('iterations','P'))

        rp <- merge(rp, DT, by=c('LOC','SEX','AGE_LABEL'))	
        rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
        rp <- rp[, list(Q= quantile(P, probs=ps), P=c('CL','M','CU')), by=c('LOC','SEX')]
        rp <- dcast.data.table(rp, LOC+SEX~P, value.var='Q')

        .r <- function(x) round(x*100, digits=2)
        .f <- function(x,y,z) paste0( .r(x),'% (', .r(y), '% - ', .r(z), '%')
        rp[, LABEL:= .f(M, CL, CU)]
        rp <- merge(unique(subset(DT, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
        prev.hiv.by.sex.loc <- copy(rp)
        # prev.hiv.by.sex.loc

        # extract prevalence ratio female:male and male:female
        # ____________________________________________________

        tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
        rp <- as.data.table(reshape2::melt(tmp))
        tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict-0.5, SEX=c(0,1), LOC=c(0,1)))
        tmp[, Var2:= seq_len(nrow(tmp))]
        rp <- merge(tmp, rp, by='Var2')
        setnames(rp, c('Var1','value'), c('iterations','P'))
        rp <- merge(rp, DT, by=c('LOC','SEX','AGE_LABEL'))	
        rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
        rp <- merge(unique(subset(DT, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
        rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations~SEX_LABEL, value.var='P')
        rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC','iterations')]
        rp <- melt(rp, id.vars=c('LOC','iterations'))[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC','variable')]
        rp <- dcast.data.table(rp, LOC+variable~P, value.var='Q')
        rp[, LABEL:= .f(M, CL, CU)]
        rp <- merge(unique(subset(DT, select=c(LOC,LOC_LABEL))), rp, by=c('LOC'))
        prevratio.hiv.by.loc <- copy(rp)

        # plot prevalence ratio F:M and M:F by age
        # ________________________________________________________

        tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
        rp <- as.data.table(reshape2::melt(tmp))
        tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
        tmp[, Var2:= seq_len(nrow(tmp))]
        rp <- merge(tmp, rp, by='Var2')
        setnames(rp, c('Var1','value'), c('iterations','P'))
        rp <- merge(rp, unique(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), by=c('LOC','SEX'))	
        rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations+AGE_LABEL~SEX_LABEL, value.var='P')
        rp[, PR_FM:= F/M]
        rp[, PR_MF:=M/F]
        rp <- melt(rp, id.vars=c('LOC_LABEL','AGE_LABEL','iterations'), measure.vars=c('PR_FM','PR_MF'))[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','AGE_LABEL','variable')]
        rp <- dcast.data.table(rp, LOC_LABEL+AGE_LABEL+variable~P, value.var='Q')
        prevratio.hiv.by.loc.age <- copy(rp)
        p <- ggplot(subset(prevratio.hiv.by.loc.age, variable=='PR_FM')) + 	
                geom_hline(yintercept=1, linetype=2) +
                geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=LOC_LABEL), alpha=0.2) +
                geom_line(aes(x=AGE_LABEL, y=M)) +
                scale_x_continuous( expand=c(0,0) ) +
                scale_y_log10(expand=c(0,0)) +
                coord_cartesian(ylim=c(.5,50)) +
                facet_wrap(~LOC_LABEL, ncol=2) +
                theme_bw() +
                labs(x='\nage at visit (years)', 
                     y='female to male HIV prevalence ratio\n(95% credibility interval)\n')

        filename=paste0('220729_hivprevalenceratio_vs_age_by_fishinland_stan_round',round,'.pdf')
        ggsave2(p, file=filename, w=8, h=5)

        # extract if diff in F:M prevalence risk ratio in fish vs inland
        # ______________________________________________________________

        tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
        rp <- as.data.table(reshape2::melt(tmp))
        tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict-0.5, SEX=c(0,1), LOC=c(0,1)))
        tmp[, Var2:= seq_len(nrow(tmp))]
        rp <- merge(tmp, rp, by='Var2')
        setnames(rp, c('Var1','value'), c('iterations','P'))
        rp <- merge(DT, rp, by=c('LOC','SEX','AGE_LABEL'))
        rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC_LABEL','SEX_LABEL','iterations')]	
        rp <- dcast.data.table(rp, LOC_LABEL+iterations~SEX_LABEL, value.var='P')
        rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC_LABEL','iterations')]
        rp <- dcast.data.table(rp, iterations~LOC_LABEL, value.var='PR_FM')
        rp[, PR_FM_D:= inland-fishing]	
        rp[, list(Q= quantile(PR_FM_D, probs=ps), P=c('CL','M','CU'))]
        # Save diff. qttities
        filename=paste0("220729f_hivprevalence_round",round,".rda")
        save(DT, re, prev.hiv.by.age, prevratio.hiv.by.loc,
             prev.hiv.by.sex.loc, prevratio.hiv.by.loc.age,
             file=file.path(vl.out.dir., filename))


        #	make table version suppressed
        .f1 <- function(x) sprintf('%2.1f', x*100)
        .f2 <- function(x) sprintf('%2.2f', x*100)
        .p1 <- function(x,y,z)  paste0(.f1(x), '(', .f1(y), '-', .f1(z))
        .p2 <- function(x,y,z)  paste0(.f2(x), '(', .f2(y), '-', .f2(z))

        prev.hiv.by.age[, LABEL:= .p1(M, CL, CU) ]
        set(prev.hiv.by.age, NULL, 'SEX_LABEL', prev.hiv.by.age[, factor(as.character(SEX_LABEL), levels=c('F','M'))])
        prevratio.hiv.by.loc.age[, LABEL2:= .p2(M, CL, CU)]
        dt <- subset(prev.hiv.by.age, AGE_LABEL %in% c(20.5,25.5,30.5,35.5,40.5,45.5))	
        dt <- dcast.data.table(dt, LOC_LABEL+AGE_LABEL~SEX_LABEL, value.var='LABEL')
        tmp <- subset(prevratio.hiv.by.loc.age,
                      variable=="PR_FM" & AGE_LABEL%in%c(20.5,25.5,30.5,35.5,40.5,45.5),
                      c(LOC_LABEL, AGE_LABEL, LABEL2))
        dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
        tmp <- subset(prevratio.hiv.by.loc.age,
                      variable=="PR_MF" & AGE_LABEL %in% c(20.5,25.5,30.5,35.5,40.5,45.5),
                      c(LOC_LABEL, AGE_LABEL, LABEL2))
        setnames(tmp, 'LABEL2', 'PR_MF')
        dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))

        filename = file.path(vl.out.dir., paste0("220729f_hivprevalence_round",round,".csv"))
        fwrite(dt, row.names=FALSE, file=filename)

        cat('Round', round, ': done.\n')
        return(TRUE)
    }
    
    if(parallelise){
        # Use parallel backend if run locally
        foreach(
                r = vla[, unique(ROUND)],
                .combine='c'
                ) %dopar% {
            cat('Running Round', r, '\n')
            .fit.stan.and.plot.by.round(vla[ ROUND ==r, ])
        } -> tmp

    }else{

        # else do not
        for( r in args$round)
            .fit.stan.and.plot.by.round(vla[ ROUND == r, ]) 
    }

    return(tmp)
}

vl.prevalence.by.gender.loc.age.icar<- function(DT)
{
        # DT <- copy(dall)
        DT <- .preprocess.ds.oli(DT)
	
        tmp <- c('N', 'HIV_N', 'VLNS_N', 'ARV_N')
        vla <- .preprocess.make.vla(DT, select=tmp)

        # Stan file locations
        file.stan.1 <- file.path(path.stan, 'vl_prevalence_by_gender_loc_age_icar_1.stan')
        file.stan.2 <- file.path(path.stan, 'vl_prevalence_by_gender_loc_age_icar_2.stan')

        .fit.stan.and.plot.by.round <- function(DT, iter=20e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999))
        {
            round <- DT[, unique(ROUND)]
            stopifnot(length(round) == 1)
            cat('Fitting stan model for round ', round, '\n')

            # Specify stan data
            #__________________
            stan.data <- list(
                              N     = nrow(vla),
                              TOTAL = vla[,N],
                              K     = vla[,HIV_N],
                              AGE_N = vla[, max(AGE)],
                              AGE   = vla[, AGE],
                              SEX   = vla[, SEX],
                              LOC   = vla[, LOC]
            )

            # second order RW prior,
            node1 =  c(vla[, seq.int(1, max(AGE)-1L)], vla[, seq.int(1, max(AGE)-2L)])
            node2 =  c(vla[, seq.int(2, max(AGE))], vla[, seq.int(3, max(AGE))])
            tmp <- sort(node1, index.return=TRUE)$ix
            stan.data$node1 <- node1[tmp]
            stan.data$node2 <- node2[tmp]
            stan.data$N_edges <-  length(stan.data$node1)


            # Fit stan model and save
            # _______________________

            stan.model <- stan_model(file=file.stan.2, model_name= 'icar2')

            fit <- sampling(stan.model, data=stan.data,
                            iter=iter, warmup=warmup, chains=chains,
                            control = control)

            #save(fit, file=file.path(vl.out.dir., "hivprevalence_icar_stanfit_200428.rda"))		# trends by age quite rough, using Cauchy prior on sigma
            #save(fit, file=file.path(vl.out.dir., "hivprevalence_icar_stanfit_200428c.rda"))	# trends by age still quite rough, using N(0,0.1) prior on sigma
            filename <- paste0("hivprevalence_icar_stanfit_round",round,"_220729.rds")
            saveRDS(fit, file=file.path(vl.out.dir., filename))

            min( summary(fit)$summary[, 'n_eff'] )
            re <- rstan::extract(fit)
            ps <- c(0.025,0.5,0.975)

            quantile(re$sigma_loc0, probs=ps)

            # make prevalence plot by age
            # ___________________________
            tmp <- apply(re$p, 2, quantile, probs=ps)
            rownames(tmp) <- c('CL','M','CU')
            tmp <- as.data.table(reshape2::melt(tmp))	
            prev.hiv.by.age <- cbind(vla, dcast.data.table(tmp, Var2~Var1, value.var='value'))	

            p <- ggplot(prev.hiv.by.age) + 		
                geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2) +
                geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL)) +
                scale_x_continuous( expand=c(0,0) ) + 
                scale_y_continuous(label=scales:::percent) +
                scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
                facet_wrap(~LOC_LABEL, ncol=2) +
                theme_bw() +
                labs(x='\nage at visit (years)', 
                     y='HIV prevalence (95% credibility interval)\n', 
                     colour='gender', 
                     linetype='location')

            filename=paste0('220729_hivprevalence_vs_age_by_gender_fishinland_icar_round',round,'.pdf')
            filename=file.path(vl.out.dir., filename)
            cat('Saving', filename, '...\n')
            filename2 <- gsub('pdf$','png',filename)
            ggsave(p, filename=filename2, w=6, h=5)
            ggsave(p, filename=filename, w=6, h=5)

            # extract basic prevalence estimates
            # __________________________________
            rp <- as.data.table(reshape2::melt( re$p ))
            setnames(rp, 2:3, c('ROW_ID','P')) 
            rp <- merge(rp, vla, by='ROW_ID')
            rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
            rp <- rp[, list(Q= quantile(P, probs=ps), P=c('CL','M','CU')), by=c('LOC','SEX')]
            rp <- dcast.data.table(rp, LOC+SEX~P, value.var='Q')
            rp[, LABEL:= paste0(round(M, d=2),'% (',round(CL, d=2),'% - ',round(CU,d=2),'%)') ]
            rp <- merge(unique(subset(vla, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
            prev.hiv.by.sex.loc <- copy(rp)

            # extract prevalence ratio female:male and male:female
            # ____________________________________________________
            rp <- as.data.table(reshape2::melt( re$p ))
            setnames(rp, 2:3, c('ROW_ID','P')) 
            rp <- merge(rp, vla, by='ROW_ID')
            rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
            rp <- merge(unique(subset(vla, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
            rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations~SEX_LABEL, value.var='P')
            rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC','iterations')]
            rp <- melt(rp, id.vars=c('LOC','iterations'))[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC','variable')]
            rp <- dcast.data.table(rp, LOC+variable~P, value.var='Q')
            rp[, LABEL:= paste0(round(M, d=2),' (',round(CL, d=2),' - ',round(CU,d=2),')') ]
            rp <- merge(unique(subset(vla, select=c(LOC,LOC_LABEL))), rp, by=c('LOC'))
            prevratio.hiv.by.loc <- copy(rp)

            # plot prevalence ratio female:male and male:female by age
            # ________________________________________________________

            rp <- as.data.table(reshape2::melt( re$p ))
            setnames(rp, 2:3, c('ROW_ID','P')) 
            rp <- merge(rp, vla, by='ROW_ID')
            rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations+AGE+AGE_LABEL~SEX_LABEL, value.var='P')
            rp[, PR_FM:= F/M]
            rp[, PR_MF:=M/F]
            rp <- melt(rp, id.vars=c('LOC_LABEL','AGE_LABEL','iterations'), measure.vars=c('PR_FM','PR_MF'))[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','AGE_LABEL','variable')]
            rp <- dcast.data.table(rp, LOC_LABEL+AGE_LABEL+variable~P, value.var='Q')
            prevratio.hiv.by.loc.age <- copy(rp)
            p <- ggplot(subset(prevratio.hiv.by.loc.age, variable=='PR_FM')) + 	
                geom_hline(yintercept=1, linetype=2) +
                geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=LOC_LABEL), alpha=0.2) +
                geom_line(aes(x=AGE_LABEL, y=M)) +
                scale_x_continuous( expand=c(0,0) ) + 			
                facet_wrap(~LOC_LABEL, ncol=2) +
                theme_bw() +
                labs(x='\nage at visit (years)', 
                     y='female to male HIV prevalence ratio\n(95% credibility interval)\n')

            filename=paste0('220729_hivprevalenceration_vs_age_by_fishinland_icar_round',round,'.pdf') 
            filename=file.path(vl.out.dir., filename)
            cat('Saving', filename, '...\n')
            filename2 <- gsub('pdf$','png',filename)
            ggsave(p, filename=filename2, w=10, h=5)
            ggsave(p, filename=filename, w=10, h=5)

            # extract if difference in F:M prevalence risk ratio in fishing vs inland
            # _______________________________________________________________________

            rp <- as.data.table(reshape2::melt( re$p ))
            setnames(rp, 2:3, c('ROW_ID','P')) 
            rp <- merge(rp, vla, by='ROW_ID')
            rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
            rp <- merge(unique(subset(vla, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
            rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations~SEX_LABEL, value.var='P')
            rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC_LABEL','iterations')]
            rp <- dcast.data.table(rp, iterations~LOC_LABEL, value.var='PR_FM')
            rp[, PR_FM_D:= inland-fishing]	
            rp[, list(Q= quantile(PR_FM_D, probs=ps), P=c('CL','M','CU'))]

            filename=file.path(vl.out.dir., paste0('hivprevalence_round',round,'220729.rda'))
            cat('Saving', filename, '...\n')
            save(vla, re, prev.hiv.by.age, prevratio.hiv.by.loc, 
                 prev.hiv.by.sex.loc, prevratio.hiv.by.loc.age,
                 file=filename)

            TRUE
        }

        foreach(
                r = vla[, unique(ROUND)],
                .combine='c'
        ) %dopar% {
                cat('Running Round', r, '\n')
                .fit.stan.and.plot.by.round(vla[ ROUND == r, ], iter=10e3)

        } -> tmp
        tmp

        return(tmp)
}

vl.meanviralload.by.gender.loc.age.icar<- function(DT) 
{

    # DT <- copy(dall)
    DT <- .preprocess.ds.oli(DT)

    tmp <- seq.int(min(DT$AGEYRS), max(DT$AGEYRS))
    tmp1 <- DT[, sort(unique(ROUND))]
    vla <- as.data.table(expand.grid(ROUND=tmp1,
                                     FC=c('fishing','inland'),
                                     SEX=c('M','F'),
                                     AGEYRS=tmp))
    vla <- vla[, {	
        z <- which(DT$ROUND==ROUND, DT$FC==FC & DT$SEX==SEX & DT$AGEYRS==AGEYRS)	
        list(N          = length(z),
             VL_MEAN    = mean(DT$VLC[z]),
             VL_SD      = sd(DT$VLC[z]),
             VL_MEAN_SD = sd(DT$VLC[z]) / sqrt(length(z)),
             HIV_N      = sum(DT$HIV_STATUS[z]==1),
             VLNS_N     = sum(DT$VLNS[z]==1),
             ARV_N      = sum(DT$ARVMED[z]==0 & DT$HIV_STATUS[z]==1 & !is.na(DT$ARVMED[z]) )
             )}, by=names(vla)]

    setnames(vla, c('FC','SEX','AGEYRS'), c('LOC_LABEL','SEX_LABEL','AGE_LABEL'))
    vla[, LOC:= as.integer(LOC_LABEL=='fishing')]
    vla[, SEX:= as.integer(SEX_LABEL=='M')]
    vla[, AGE:= AGE_LABEL-14L]
    vla[, ROW_ID:= seq_len(nrow(vla))]

    p <- ggplot(vla, aes(x=AGE_LABEL, fill=SEX_LABEL, color=SEX_LABEL) ) + 		
        geom_ribbon(aes(x=AGE_LABEL, ymin=VL_MEAN-2*VL_MEAN_SD, ymax=VL_MEAN+2*VL_MEAN_SD, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2, colour=NA) +
        geom_line(aes(x=AGE_LABEL, y=VL_MEAN, colour=SEX_LABEL)) +
        scale_x_continuous( expand=c(0,0) ) + 
        scale_y_log10() +
        scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
        scale_fill_manual(values=c('M'='royalblue3','F'='deeppink2')) +
        facet_grid(ROUND~LOC_LABEL, scales='free') +
        theme_bw() +
        labs(x='\nage at visit (years)', 
             y='mean viral load\n(95% credibility interval)\n', 
             colour='gender', fill='gender', linetype='location')

    filename = '220729d_mvl_vs_age_by_gender_fishinland_raw.pdf'
    filename = file.path(vl.out.dir., filename)
    cat('Saving', filename, '...\n')
    filename2 <- gsub('pdf$','png',filename)
    ggsave(p, filename=filename2, w=6, h=10)
    ggsave(p, filename=filename, w=6, h=10)

    # Stan file locations
    file.stan.1 <- file.path(path.stan, 'vl_meanviralload_by_gender_loc_age_icar.stan')
    file.stan.2 <- file.path(path.stan, 'vl_meanviralload_by_gender_loc_age_icar2.stan')

    stan.model <- stan_model(file=file.stan.1, model_name= 'icar_age_interactions')

    .fit.stan.and.plot.by.round <- function(DT , iter=20e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999))
    {
        # DT <- copy(vla[ROUND == 16] )
        round <- DT[, unique(ROUND)]
        stopifnot(length(round) == 1)
        cat('Fitting stan model for round ', round, '\n')
        DT[is.na(VL_MEAN_SD), VL_MEAN_SD := 0]

        # second order RW prior
        node1 <-  c(DT[, seq.int(1, max(AGE)-1L)], DT[, seq.int(1, max(AGE)-2L)])
        node2 <-  c(DT[, seq.int(2, max(AGE))], DT[, seq.int(3, max(AGE))])
        tmp <- sort(node1, index.return=TRUE)$ix

        stan.data <- list(
                          N       = nrow(DT),
                          MEAN    = DT[,VL_MEAN],
                          SD      = pmax(1, DT[,VL_MEAN_SD]),
                          AGE_N   = DT[, max(AGE)],
                          AGE     = DT[, AGE],
                          SEX     = DT[, SEX],
                          LOC     = DT[, LOC],
                          node1 = node1[tmp],
                          node2 = node2[tmp],
                          N_edges = length(tmp)
        )

        fit <- sampling(stan.model, data=stan.data, iter=iter, warmup=warmup, chains=chains, control=control)

        filename=paste0('mvlinpop_icar_stanfit_round',round,'220729b.rds')
        cat('Saving', filename, '...\n')
        saveRDS(fit, file=file.path(vl.out.dir., filename) )

        min( summary(fit)$summary[, 'n_eff'] )
        re <- rstan::extract(fit)
        ps <- c(0.025,0.5,0.975)

        #	make prevalence plot by age
        tmp <- apply(re$mu, 2, quantile, probs=ps)
        rownames(tmp) <- c('CL','M','CU')
        tmp <- as.data.table(reshape2::melt(tmp))	
        mvl.by.age <- cbind(vla, dcast.data.table(tmp, Var2~Var1, value.var='value'))	

        p <- ggplot(mvl.by.age) + 		
            geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2) +
            geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL)) +
            scale_x_continuous( expand=c(0,0) ) + 
            #scale_y_log10() +
            scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
            facet_wrap(~LOC_LABEL, ncol=2) +
            theme_bw() +
            labs(x='\nage at visit (years)', 
                 y='mean viral load\n(95% credibility interval)\n', 
                 colour='gender', 
                 linetype='location')

        filename=paste0('2207289_mvl_vs_age_by_gender_fishinland_stan_round',round,'.pdf')
        filename=file.path(vl.out.dir., filename)
        filename2 <- gsub('pdf$','png',filename)
        ggsave(p, filename=filename2, w=6, h=5)
        ggsave(p, filename=filename, w=6, h=5)

        TRUE
    }

    if(parallelise)
    {
        # Use parallel backend if run locally
        foreach(
                r = vla[, unique(ROUND)],
                .combine='c'
                ) %dopar% {
            cat('Running Round', r, '\n')
            .fit.stan.and.plot.by.round(vla[ ROUND ==r, ])
        } -> tmp

    }else{

        # else do not
        for( r in args$round)
               .fit.stan.and.plot.by.round(vla[ ROUND == r, ]) 
    }

    return(tmp)
}

vl.suppofinfected.by.gender.loc.age.gp<- function(DT, refit=FALSE, vl.out.dir.=vl.out.dir)
{
    cat('\n\n--- Analyse suppressed among infected ---\n\n')

    # DT <- copy(dall); refit=FALSE; vl.out.dir.=vl.out.dir
    DT  <- .preprocess.ds.oli(DT)
    vla <- .preprocess.make.vla(DT, select=c('N', 'HIV_N', 'VLNS_N', 'ARV_N'))
    # NOTE: ARV_N == 0

    file.stan.1 <- file.path(path.stan, 'vl_binomial_gp.stan')

    .fit.stan.and.plot.by.round <- function(DT, itr=10e3, wrmp=5e2, chns=1, cntrl = list(max_treedepth= 15, adapt_delta= 0.999))
    {
        # DT <- copy(vla[ROUND == 16] )
        round <- DT[, unique(ROUND)]
        stopifnot(length(round) == 1)
        cat('Fitting stan model for round ', round, '\n')

        DT[, VLSUP_N := HIV_N - VLNS_N,]
        stan.data <- .make.stan.data.gp(DT, num.var='VLSUP_N', den.var='HIV_N')
        filename <- paste0( '220729f_notsuppAmongInfected_gp_stan_round',round,'.rds')
        filename <- file.path(vl.out.dir., filename)

        if(file.exists(filename) & refit == FALSE)
        {
            cat('Loading previously run HMC... \n')
            fit <- readRDS(filename)
        }else{
            stan.model <- stan_model(file.stan.1, model_name='gp_all')	 
            fit <- sampling(stan.model, data=stan.data, 
                            iter=itr, warmup=wrmp,
                            chains=chns, control=cntrl)
            saveRDS(fit, file=filename)
        }

        # compare to self-report
        # ______________________

        arv_bool <- DT[, any(ARV_N) > 0]
        if(arv_bool)
        {
            DT[, HIV_NARV_N := HIV_N-ARV_N]
            stan.data.2 <- .make.stan.data.gp(DT, num.var='HIV_NARV_N', den.var='HIV_N')

            filename <- paste0('220729f_notARVAmongInfected_gp_stan_round',round,'.rds')
            filename <- file.path(vl.out.dir., filename)

            if(file.exists(filename) & refit == FALSE)
            {
                cat('Loading previously run HMC... \n')
                fit2 <- readRDS(filename)
            }else{
                fit2 <- sampling(stan.model, data=stan.data, 
                                 iter=iter, warmup=warmup,
                                 chains=chains, control=control)
                saveRDS(fit2, file=filename)
            }
        }

        # Analyse posterior
        # _________________

        # Extract summary excluding hyper parameters
        dsum <- summary(fit)$summary
        idx <- rownames(summary(fit)$summary) %like% 'rho_hyper_par'
        dsum <- dsum[!idx,]
        
        cat('The minimum effective sample size is:\n')
        cat(min( dsum[, 'n_eff'], na.rm=T ), '\n')
        cat('Out ot ', nrow(dsum), 'parameters\n')
        cat(sum(dsum[, 'n_eff' ] < 4000, na.rm=TRUE), ' had n_eff < 4000,\n')
        cat(sum(dsum[, 'n_eff' ] < 1000, na.rm=TRUE), ' had n_eff < 1000.\n')

        re <- rstan::extract(fit)
        # re2 <- rstan::extract(fit2)
        ps <- c(0.025,0.5,0.975)
        tmp <- summary(fit)$summary
        tmp[grepl('^p_predict_',rownames(tmp)),]

        # extract hyperparams rho
        # _______________________

        ps <- c(0.025,0.25,0.5,0.75,0.975)
        .f <- function(x) transpose(as.data.table(quantile(x, probs=ps)))

        tmp <- c('rho_00','rho_10','rho_01','rho_11','alpha_00','alpha_10','alpha_01','alpha_11')
        tmp <- lapply(re[tmp], .f)
        tmp <- rbindlist(tmp, idcol='GP_hyper_par')
        setnames(tmp, paste0('V', 1:5), c('CL','IL','M','IU','CU') )

        .f <- function(reg, x) gsub('^([a-z]+)_([0-9])([0-9])',reg,x)
        tmp[, SEX:= as.integer(.f( '\\2', GP_hyper_par))]
        tmp[, LOC:= as.integer(.f( '\\3', GP_hyper_par))]
        tmp[, GP_hyper_par:= .f('\\1',GP_hyper_par)]

        tmp1 <- unique(DT[, .(SEX,SEX_LABEL,LOC,LOC_LABEL)])
        nsinf.gp.pars <- merge(tmp1, tmp, by=c('SEX','LOC'))
        nsinf.gp.pars[, col:='posterior']
        tmp <- .get.prior.ranges(stan.data, DT=tmp1,
                                 shape=re$rho_hyper_par_shape2[1],
                                 scale=re$rho_hyper_par_scale2[1])
        p <- .plot.gp.hyperparameters(nsinf.gp.pars, tmp)
        filename <- paste0('220729f_notsuppAmongInfected_gppars_round',round,'.pdf')
        ggsave2(p, file=filename, w=6, h=3)


        # make use stan.data for PPC.
        # ___________________________

        ppDT <- copy(DT)
        cols <- c('M', 'CU', 'CL')
        ppDT[, (cols) := binconf(HIV_N - VLNS_N, HIV_N, return.df=T)]

        # make prevalence plot by age
        # ___________________________

        .f <- function(x) apply(x, 2, quantile, probs=ps)
        tmp <- cbind(.f(re$p_predict_00), .f(re$p_predict_10),
                     .f(re$p_predict_01), .f(re$p_predict_11))

        rownames(tmp) <- c('CL','IL','M','IU','CU')
        tmp <- as.data.table(reshape2::melt(tmp))

        nsinf.by.age <- dcast.data.table(tmp, Var2~Var1, value.var='value')
        tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
        nsinf.by.age <- cbind(tmp, nsinf.by.age) 
        nsinf.by.age <- merge(unique(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), nsinf.by.age, by=c('SEX','LOC'))
        nsinf.by.age[, STAT:='VLNS']

        nsinf.by.age2 <- copy(nsinf.by.age)
        if(arv_bool)
        {
            .f <- function(x) apply(x, 2, quantile, probs=ps)
            tmp <- cbind( .f(re2$p_predict_00), .f(re2$p_predict_10),
                         .f(re2$p_predict_01), .f(re2$p_predict_11))
            rownames(tmp) <- c('CL','IL','M','IU','CU')
            tmp <- as.data.table(reshape2::melt(tmp))
            nainf.by.age <- dcast.data.table(tmp, Var2~Var1, value.var='value')
            tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
            nainf.by.age <- cbind(tmp, nainf.by.age) 

            tmp1 <- unique(DT[, .(SEX,SEX_LABEL,LOC,LOC_LABEL)])
            nainf.by.age <- merge(tmp1, nainf.by.age, by=c('SEX','LOC'))
            nainf.by.age[, STAT:='VLNA']

            tmp <- subset(nainf.by.age, select=c(SEX, LOC, AGE_LABEL, M, CL, CU))
            setnames(tmp, c('M','CL','CU'), c('M2','CL2','CU2'))
            nsinf.by.age2 <- merge(nsinf.by.age, tmp, by=c('SEX','LOC','AGE_LABEL'))			
        }

        p <- ggplot(nsinf.by.age2, aes(x=AGE_LABEL, colour=SEX_LABEL)) + 		
            geom_ribbon(aes(ymin=CL, ymax=CU, fill=SEX_LABEL), alpha=0.2, col=NA) +
            geom_line(aes(y=M)) +
            geom_hline(yintercept=c(0.9^3, 0.95^3), linetype='dashed') +
            scale_x_continuous( expand=c(0,0) ) + 
            scale_y_continuous(labels=scales:::percent, expand=c(0,0), limits=c(0,1)) +
            scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
            scale_fill_manual(values=c('M'='royalblue3','F'='deeppink2')) +
            facet_wrap(~LOC_LABEL, ncol=2) +
            theme_bw() +
            theme(legend.position='bottom') + 
            labs(x='\nage at visit (years)', 
                 y='HIV+ individuals with suppressed viral load\n(95% credibility interval)\n', 
                 colour='gender', fill='gender')

        filename <- paste0('220729f_suppAmongInfected_vs_age_by_gender_fishinland_gp_round',round,'.pdf')
        ggsave2(p, file=filename, w=6, h=5)

        p <- p + 
            geom_point(data=ppDT, aes(y=M, pch=SEX_LABEL, size=N)) +
            # geom_linerange(data=ppDT, aes(ymin=CL, ymax=CU), linetype='dotted' ) +
            scale_size(range = c(0, 3)) +
            labs(pch='gender', linetype='gender', size='population size')

        filename <- paste0('220729f_suppAmongInfected_vs_age_by_gender_fishinland_data_gp_round',round,'.pdf')
        ggsave2(p, file=filename, w=6, h=5)


        # p <- ggplot(tmp) + 		
        #         geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, fill=SEX_LABEL), alpha=0.2) +
        #         geom_hline(yintercept=c(0.9^3, 0.95^3)) +
        #         geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL)) +
        #         geom_line(aes(x=AGE_LABEL, y=M2, colour=SEX_LABEL), linetype=2) +
        #         scale_x_continuous( expand=c(0,0) ) + 
        #         scale_y_continuous(labels=scales:::percent, expand=c(0,0), limits=c(0,1)) +
        #         scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
        #         scale_fill_manual(values=c('M'='royalblue3','F'='deeppink2')) +
        #         facet_wrap(~LOC_LABEL, ncol=2) +
        #         theme_bw() +
        #         theme(legend.position='bottom') +
        #         labs(x='\nage at visit (years)', 
        #              y='HIV+ individuals with suppressed viral load\n(95% credibility interval)\n', 
        #              colour='gender', fill='gender')

        # filename <- paste0('220729f_suppAmongInfected_vs_age_by_gender_fishinland_stan_v2_round',round,'.pdf')
        # ggsave2(p, file=filename, w=6, h=5)	

        tmp <- copy( nsinf.by.age)
        if(arv_bool)
        {
            tmp <- rbind(nsinf.by.age, nainf.by.age, fill=TRUE)
        }

        p <- ggplot(tmp, aes(x=AGE_LABEL, colour=SEX_LABEL)) + 		
            geom_ribbon(aes(ymin=CL, ymax=CU, fill=SEX_LABEL), colour=NA, alpha=0.2) +
            geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL, linetype=STAT)) +
            scale_x_continuous( expand=c(0,0) ) + 
            scale_y_continuous(labels=scales:::percent) +
            scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
            scale_fill_manual(values=c('M'='royalblue3','F'='deeppink2')) +
            scale_linetype_manual(values=c('VLNS'='solid','VLNA'='dotdash')) +
            facet_grid(SEX_LABEL~LOC_LABEL) +
            theme_bw() +
            theme(legend.position='bottom')
        labs(x='\nage at visit (years)', 
             y='HIV+ individuals with unsuppressed viral load\n(95% credibility interval)\n', 
             colour='gender', fill='gender')

        # p <- p + 
        #         geom_point(data=ppDT, aes(y=M, pch=SEX_LABEL, size=N)) +
        #         # geom_linerange(data=ppDT, aes(ymin=CL, ymax=CU), linetype='dotted' ) +
        #         scale_size(range = c(0, 3)) +
        #         labs(pch='gender', linetype='gender', size='population size')

        filename <- paste0('220729f_suppAmongInfected_vs_age_by_gender_fishinland_gp_v3_round',round,'.pdf')
        ggsave2(p, file=filename, w=9, h=8)

        # extract basic not supp estimates
        # ________________________________

        DT[, list(N=sum(VLNS_N), 
                  P=sum(VLNS_N) / sum(HIV_N),
                  N2=sum(HIV_N)-sum(VLNS_N), 
                  P2= 1-sum(VLNS_N) / sum(HIV_N)), by=c('LOC_LABEL','SEX_LABEL')]

        ps <- c(0.025,0.5,0.975)
        tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
        rp <- as.data.table(reshape2::melt( tmp ))
        setnames(rp, 1:3, c('iterations','ROW_ID','P')) 
        tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
        tmp[, ROW_ID:= seq_len(nrow(tmp))]
        tmp <- merge(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL,N)), tmp, by=c('SEX','LOC','AGE_LABEL'), all=TRUE)
        set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N]) 
        set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N])
        set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
        set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
        set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])
        set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])	
        rp <- merge(tmp, rp, by=c('ROW_ID'))	
        rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC_LABEL','SEX_LABEL','iterations')]
        rp <- rp[, list(Q= quantile(P, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','SEX_LABEL')]
        rp <- dcast.data.table(rp, LOC_LABEL+SEX_LABEL~P, value.var='Q')
        rp[, LABEL:= paste0(round(M*100, d=2),'% (',round(CL*100, d=2),'% - ',round(CU*100,d=2),'%)') ]
        nsinf.by.sex.loc <- copy(rp)


        # extract risk ratio of suppressed VL female:male and male:female
        # _______________________________________________________________

        tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
        rp <- as.data.table(reshape2::melt( tmp ))
        setnames(rp, 1:3, c('iterations','ROW_ID','P')) 
        tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
        tmp[, ROW_ID:= seq_len(nrow(tmp))]
        tmp <- merge(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL,N)), tmp, by=c('SEX','LOC','AGE_LABEL'), all=TRUE)
        set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N]) 
        set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N])
        set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
        set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
        set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])
        set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])	
        rp <- merge(tmp, rp, by=c('ROW_ID'))	
        rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC_LABEL','SEX_LABEL','iterations')]
        rp <- dcast.data.table(rp, LOC_LABEL+iterations~SEX_LABEL, value.var='P')
        rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC_LABEL','iterations')]
        rp <- melt(rp, id.vars=c('LOC_LABEL','iterations'))
        rp <- rp[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','variable')]
        rp <- dcast.data.table(rp, LOC_LABEL+variable~P, value.var='Q')
        rp[, LABEL:= .write.CIs(M, CL, CU, d=2)]
        nsinf.by.loc <- copy(rp)	

        # extract risk ratio of unsuppressed VL female:male and male:female by age
        # ________________________________________________________________________

        tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
        rp <- as.data.table(reshape2::melt( tmp ))
        setnames(rp, 1:3, c('iterations','ROW_ID','P')) 
        tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
        tmp[, ROW_ID:= seq_len(nrow(tmp))]
        tmp <- merge(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL,N)), tmp, by=c('SEX','LOC','AGE_LABEL'), all=TRUE)
        set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N]) 
        set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N])
        set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
        set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
        set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])
        set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])	
        rp <- merge(tmp, rp, by=c('ROW_ID'))	
        rp <- dcast.data.table(rp, LOC_LABEL+iterations+AGE_LABEL~SEX_LABEL, value.var='P')
        rp[, PR_FM:= F/M]
        rp[, PR_MF:= M/F]
        rp <- melt(rp, id.vars=c('LOC_LABEL','AGE_LABEL','iterations'), measure.vars=c('PR_FM','PR_MF'))
        rp <- rp[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','AGE_LABEL','variable')]	
        rp <- dcast.data.table(rp, LOC_LABEL+AGE_LABEL+variable~P, value.var='Q')
        rp[, LABEL:= .write.CIs(M, CL, CU, d=2)]
        nsinf.ratio.by.loc.age <- copy(rp)


        # extract if F:M risk ratio of unsupp VL is diff in fishing vs inland
        # ___________________________________________________________________

        tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
        rp <- as.data.table(reshape2::melt( tmp ))
        setnames(rp, 1:3, c('iterations','ROW_ID','P')) 
        tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
        tmp[, ROW_ID:= seq_len(nrow(tmp))]
        tmp <- merge(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL,N)), tmp, by=c('SEX','LOC','AGE_LABEL'), all=TRUE)
        set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N]) 
        set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N])
        set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
        set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
        set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])
        set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])	
        rp <- merge(tmp, rp, by=c('ROW_ID'))	
        rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC_LABEL','SEX_LABEL','iterations')]
        rp <- dcast.data.table(rp, LOC_LABEL+iterations~SEX_LABEL, value.var='P')
        rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC_LABEL','iterations')]
        rp <- dcast.data.table(rp, iterations~LOC_LABEL, value.var='PR_FM')
        rp[, PR_FM_D:= fishing-inland]	
        rp[, list(Q= quantile(PR_FM_D, probs=ps), P=c('CL','M','CU'))]


        filename <- paste0('220729f_suppAmongInfected_round',round,'.rda')

        if(!arv_bool)
        {
            nainf.by.age <- re2 <- data.table()

        }
        save(DT, re, re2, nainf.by.age,
             nsinf.by.age, nsinf.by.sex.loc,
             nsinf.by.loc, nsinf.ratio.by.loc.age,
             file=file.path(vl.out.dir.,filename))


        # make table version suppressed
        # _____________________________

        nsinf.by.age[, LABEL:= paste0(sprintf('%2.1f',M*100),' (',sprintf('%2.1f',CL*100),' - ',sprintf('%2.1f',CU*100),')') ]
        set(nsinf.by.age, NULL, 'SEX_LABEL', nsinf.by.age[, factor(as.character(SEX_LABEL), levels=c('F','M'))])	
        nsinf.ratio.by.loc.age[, LABEL2:= paste0(sprintf('%2.2f',M),' (',sprintf('%2.2f',CL),' - ',sprintf('%2.2f',CU),')') ]
        dt <- subset(nsinf.by.age, AGE_LABEL%in%c(20.5,25.5,30.5,35.5,40.5,45.5))	
        dt <- dcast.data.table(dt, LOC_LABEL+AGE_LABEL~SEX_LABEL, value.var='LABEL')
        tmp <- subset(nsinf.ratio.by.loc.age, variable=="PR_FM" & AGE_LABEL%in%c(20.5,25.5,30.5,35.5,40.5,45.5), c(LOC_LABEL, AGE_LABEL, LABEL2))
        dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
        tmp <- subset(nsinf.ratio.by.loc.age, variable=="PR_MF" & AGE_LABEL%in%c(20.5,25.5,30.5,35.5,40.5,45.5), c(LOC_LABEL, AGE_LABEL, LABEL2))
        setnames(tmp, 'LABEL2', 'PR_MF')
        dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))

        filename <- paste0('220729f_suppamonginfected_round',round, '.csv')
        fwrite(dt, row.names=FALSE, file=file.path(vl.out.dir.,filename))	

        TRUE
    }

    if(parallelise)
    {
        # Use parallel backend if run locally
        foreach(
                r = vla[, unique(ROUND)],
                .combine='c'
                ) %dopar% {
            cat('Running Round', r, '\n')
            .fit.stan.and.plot.by.round(vla[ ROUND ==r, ])
        } -> tmp

    }else{
        # else do not
        for( r in args$round)
            .fit.stan.and.plot.by.round(vla[ ROUND == r, ]) 
    }

    return(tmp)
}

vl.suppofinfected.by.gender.loc.age.icar<- function(DT, refit=FALSE)
{
    # DT <- copy(dall)
	vl.out.dir <- file.path(vl.out.dir)
    DT <- .preprocess.ds.oli(DT)

    vla <- .preprocess.make.vla(DT, select=c('N', 'HIV_N', 'VLNS_N', 'ARV_N'))

    # Stan file locations
    file.stan.1 <- file.path(path.stan, 'vl_suppofinfected_by_gender_loc_age_icar_1.stan')

    # list.files(path.stan, pattern='suppofinfected')
    stan.model1 <- stan_model(file.stan.1, model_name= 'icar_age_interactions')

    .fit.stan.and.plot.by.round <- function(DT , iter=20e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999))
    {
        # DT <- copy(vla[ROUND == 16] )
        round <- DT[, unique(ROUND)]
        stopifnot(length(round) == 1)
        cat('Fitting stan model for round ', round, '\n')

        # second order RW prior
        node1 <-  c(DT[, seq.int(1, max(AGE)-1L)], DT[, seq.int(1, max(AGE)-2L)])
        node2 <-  c(DT[, seq.int(2, max(AGE))], DT[, seq.int(3, max(AGE))])
        tmp <- sort(stan.data$node1, index.return=TRUE)$ix

        stan.data <- list(
                          N = nrow(DT),
                          TOTAL = DT[,HIV_N],
                          K = DT[,VLNS_N],
                          AGE_N = DT[, max(AGE)],
                          AGE = DT[, AGE],
                          SEX = DT[, SEX],
                          LOC = DT[, LOC],
                          node1 = node1[tmp],
                          node2 = node2[tmp],
                          N_edges =  length(node1)
        )
        fit <- sampling(stan.model, data=stan.data, iter=iter, warmup=warmup, chains=chains, control=control) 
        # trends by age quite rough, using Cauchy prior on sigma
        # trends by age still quite rough, using N(0,0.1) prior on sigma

        filename=paste0('notsuppAmongInfected_icar_stan_round',round,'_220729.rds')
        filename=file.path(vl.out.dir, filename)

        if(file.exists(filename) & refit == FALSE)
        {
            cat('Loading previously run HMC... \n')
            fit <- readRDS(filename)
        }else{
            fit <- sampling(stan.model, data=stan.data, 
                            iter=iter, warmup=warmup,
                            chains=chains, control=control)
            saveRDS(fit, file=filename)
        }

        # second order RW prior
        node1 <-  c(DT[, seq.int(1, max(AGE)-1L)], DT[, seq.int(1, max(AGE)-2L)])
        node2 <-  c(DT[, seq.int(2, max(AGE))], DT[, seq.int(3, max(AGE))])
        tmp <- sort(node1, index.return=TRUE)$ix

        stan.data <- list(
                          N = nrow(DT),
                          TOTAL = DT[,HIV_N],
                          K = DT[,ARV_N],
                          AGE_N = DT[, max(AGE)],
                          AGE = DT[, AGE],
                          SEX = DT[, SEX],
                          LOC = DT[, LOC],
                          node1 = node1[tmp],
                          node2 = node2[tmp],
                          N_edges =  length(node1)
        )

        filename=paste0('notARVAmongInfected_icar_stan_round',round,'_220729.rds')
        filename=file.path(vl.out.dir, filename)

        if(file.exists(filename) & refit == FALSE)
        {
            cat('Loading previously run HMC... \n')
            fit2 <- readRDS(filename)
        }else{
            fit2 <- sampling(stan.model, data=stan.data, 
                             iter=iter, warmup=warmup,
                             chains=chains, control=control)
            saveRDS(fit2, file=filename)
        }

        re <- rstan::extract(fit)
        re2 <- rstan::extract(fit2)
        ps <- c(0.025,0.5,0.975)


        # make prevalence plot by age
        # ___________________________
        tmp <- apply(re$p, 2, quantile, probs=ps)
        rownames(tmp) <- c('CL','M','CU')
        tmp <- as.data.table(reshape2::melt(tmp))	
        nsprev.by.age <- cbind(DT, dcast.data.table(tmp, Var2~Var1, value.var='value'))
        nsprev.by.age[, STAT:='VLNS']
        tmp <- apply(re2$p, 2, quantile, probs=ps)
        rownames(tmp) <- c('CL','M','CU')
        tmp <- as.data.table(reshape2::melt(tmp))	
        naprev.by.age <- cbind(DT, dcast.data.table(tmp, Var2~Var1, value.var='value'))
        naprev.by.age[, STAT:='VLNA']
        tmp <- subset(naprev.by.age, select=c(ROW_ID, M, CL, CU))
        setnames(tmp, c('M','CL','CU'), c('M2','CL2','CU2'))
        tmp <- merge(nsprev.by.age, tmp, by='ROW_ID')	

        p <- ggplot(tmp) + 		
            geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2) +
            # Do we need the interaction above?
            geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL)) +
            scale_x_continuous( expand=c(0,0) ) + 
            scale_y_continuous(label=scales:::percent) +
            scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
            # scale_fill_manual(values=c('M'='royalblue3','F'='deeppink2')) +
            facet_wrap(~LOC_LABEL, ncol=2) +
            theme_bw() +
            labs(x='\nage at visit (years)', 
                 y='HIV+ individuals with unsuppressed viral load\n(95% credibility interval)\n', 
                 colour='gender')

        filename <- paste0('220729d_notsuppAmongInfected_vs_age_by_gender_fishinland_icar_v1_round',round,'.pdf')
        ggsave(p, file=filename, w=6, h=5)	

        p <- ggplot(tmp) + 		
            geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2) +
            geom_hline(yintercept=c(1-0.9^3, 1-0.95^3)) +
            geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL)) +
            geom_line(aes(x=AGE_LABEL, y=M2, colour=SEX_LABEL), linetype=2) +
            scale_x_continuous( expand=c(0,0) ) + 
            scale_y_continuous(label=scales:::percent) +
            scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
            facet_wrap(~LOC_LABEL, ncol=2) +
            theme_bw() +
            labs(x='\nage at visit (years)', 
                 y='HIV+ individuals with unsuppressed viral load\n(95% credibility interval)\n', 
                 colour='gender')

        filename <- paste0('220729d_notsuppAmongInfected_vs_age_by_gender_fishinland_icar_v2_round',round,'.pdf')
        ggsave2(p, file=filename, w=6, h=5)	

        p <- ggplot(tmp) + 		
            geom_ribbon(aes(x=AGE_LABEL, ymin=1-CU, ymax=1-CL, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2) +
            geom_hline(yintercept=c(0.9^3, 0.95^3)) +
            geom_line(aes(x=AGE_LABEL, y=1-M, colour=SEX_LABEL)) +
            geom_line(aes(x=AGE_LABEL, y=1-M2, colour=SEX_LABEL), linetype=2) +
            scale_x_continuous( expand=c(0,0) ) + 
            scale_y_continuous(label=scales:::percent) +
            scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
            facet_wrap(~LOC_LABEL, ncol=2) +
            theme_bw() +
            labs(x='\nage at visit (years)', 
                 y='HIV+ individuals with suppressed viral load\n(95% credibility interval)\n', 
                 colour='gender')
        filename <- paste0('220729d_suppAmongInfected_vs_age_by_gender_fishinland_icar_v2_round',round,'.pdf')
        ggsave2(p, file=file.path(vl.out.dir,filename), w=6, h=5)

        tmp <- rbind(nsprev.by.age, naprev.by.age, fill=TRUE)

        p <- ggplot(tmp) + 		
            geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=interaction(SEX_LABEL,LOC_LABEL,STAT)), alpha=0.2) +
            geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL, linetype=STAT)) +
            scale_x_continuous( expand=c(0,0) ) + 
            scale_y_continuous(label=scales:::percent) +
            scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
            scale_linetype_manual(values=c('VLNS'='solid','VLNA'='dotdash')) +
            facet_grid(SEX_LABEL~LOC_LABEL) +
            theme_bw() +
            labs(x='\nage at visit (years)', 
                 y='HIV+ individuals with unsuppressed viral load\n(95% credibility interval)\n', 
                 colour='gender')

        filename <- paste0('220729d_notsuppAmongInfected_vs_age_by_gender_fishinland_icar_v3_round',round,'.pdf')
        ggsave(p, file=filename, w=9, h=8)
        
        p <- ggplot(tmp) + 		
            geom_ribbon(aes(x=AGE_LABEL, ymin=1-CU, ymax=1-CL, group=interaction(SEX_LABEL,LOC_LABEL,STAT)), alpha=0.2) +
            geom_line(aes(x=AGE_LABEL, y=1-M, colour=SEX_LABEL, linetype=STAT)) +
            scale_x_continuous( expand=c(0,0) ) + 
            scale_y_continuous(label=scales:::percent) +
            scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
            scale_linetype_manual(values=c('VLNS'='solid','VLNA'='dotdash')) +
            facet_grid(SEX_LABEL~LOC_LABEL) +
            theme_bw() +
            labs(x='\nage at visit (years)', 
                 y='HIV+ individuals with suppressed viral load\n(95% credibility interval)\n', 
                 colour='gender')
                
        filename <- paste0('220729d_notsuppAmongInfected_vs_age_by_gender_fishinland_v4_round',round,'.pdf')
        ggsave(p, file=vl.out.dir, w=9, h=8)
                

        # extract basic not supp estimates
        # ________________________________

        DT[, list(N=sum(VLNS_N), P=sum(VLNS_N) / sum(HIV_N), N2=sum(HIV_N)-sum(VLNS_N), P2= 1-sum(VLNS_N) / sum(HIV_N)), by=c('LOC_LABEL','SEX_LABEL')]

        rp <- as.data.table(reshape2::melt( re$p ))
        setnames(rp, 2:3, c('ROW_ID','P')) 
        rp <- merge(rp, DT, by='ROW_ID')
        rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
        rp <- rp[, list(Q= quantile(P, probs=ps), P=c('CL','M','CU')), by=c('LOC','SEX')]
        rp <- dcast.data.table(rp, LOC+SEX~P, value.var='Q')
        rp[, LABEL:= paste0(round(M, d=2),'% (',round(CL, d=2),'% - ',round(CU,d=2),'%)') ]
        rp <- merge(unique(subset(DT, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
        nsprev.by.sex.loc <- copy(rp)

        # extract risk ratio of unsuppressed VL female:male and male:female
        rp <- as.data.table(reshape2::melt( re$p ))
        setnames(rp, 2:3, c('ROW_ID','P')) 
        rp <- merge(rp, DT, by='ROW_ID')
        rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
        rp <- merge(unique(subset(DT, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
        rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations~SEX_LABEL, value.var='P')
        rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC','iterations')]
        rp <- melt(rp, id.vars=c('LOC','iterations'))[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC','variable')]
        rp <- dcast.data.table(rp, LOC+variable~P, value.var='Q')
        rp[, LABEL:= paste0(round(M, d=2),' (',round(CL, d=2),' - ',round(CU,d=2),')') ]
        rp <- merge(unique(subset(DT, select=c(LOC,LOC_LABEL))), rp, by=c('LOC'))
        nsprev.by.loc <- copy(rp)


        # extract risk ratio of unsuppressed VL female:male and male:female by age
        rp <- as.data.table(reshape2::melt( re$p ))
        setnames(rp, 2:3, c('ROW_ID','P')) 
        rp <- merge(rp, DT, by='ROW_ID')
        rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations+AGE+AGE_LABEL~SEX_LABEL, value.var='P')
        rp[, PR_FM:= F/M]
        rp[, PR_MF:=M/F]
        rp <- melt(rp, id.vars=c('LOC_LABEL','AGE_LABEL','iterations'), measure.vars=c('PR_FM','PR_MF'))[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','AGE_LABEL','variable')]
        rp <- dcast.data.table(rp, LOC_LABEL+AGE_LABEL+variable~P, value.var='Q')
        nsprev.ratio.by.loc.age <- copy(rp)

        # extract risk ratio of suppressed VL female:male and male:female by age
        rp <- as.data.table(reshape2::melt( re$p ))
        setnames(rp, 2:3, c('ROW_ID','P')) 
        rp <- merge(rp, DT, by='ROW_ID')
        rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations+AGE+AGE_LABEL~SEX_LABEL, value.var='P')
        rp[, PR_FM:= (1-F)/(1-M)]
        rp[, PR_MF:=(1-M)/(1-F)]
        rp <- melt(rp, id.vars=c('LOC_LABEL','AGE_LABEL','iterations'), measure.vars=c('PR_FM','PR_MF'))[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','AGE_LABEL','variable')]
        rp <- dcast.data.table(rp, LOC_LABEL+AGE_LABEL+variable~P, value.var='Q')
        sprev.ratio.by.loc.age <- copy(rp)

        # extract if difference in female:male risk ratio of unsuppressed VL is different in fishing vs inland
        rp <- as.data.table(reshape2::melt( re$p ))
        setnames(rp, 2:3, c('ROW_ID','P')) 
        rp <- merge(rp, DT, by='ROW_ID')
        rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
        rp <- merge(unique(subset(DT, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
        rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations~SEX_LABEL, value.var='P')
        rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC_LABEL','iterations')]
        rp <- dcast.data.table(rp, iterations~LOC_LABEL, value.var='PR_FM')
        rp[, PR_FM_D:= fishing-inland]	
        rp[, list(Q= quantile(PR_FM_D, probs=ps), P=c('CL','M','CU'))]

        filename <- file.path(vl.out.dir, paste0("notsuppamonginfected_220729_round",round,".rda"))
        save(DT, re, 
             sprev.ratio.by.loc.age,
             nsprev.by.age, naprev.by.age,
             nsprev.by.loc, nsprev.by.sex.loc, 
             nsprev.ratio.by.loc.age,
             file=filename)

        #	make table version unsuppressed
        nsprev.by.age[, LABEL:= paste0(round(M*100, d=1),' (',round(CL*100, d=1),' - ',round(CU*100,d=1),')') ]
        set(nsprev.by.age, NULL, 'SEX_LABEL', nsprev.by.age[, factor(as.character(SEX_LABEL), levels=c('F','M'))])	
        nsprev.ratio.by.loc.age[, LABEL2:= paste0(round(M, d=2),' (',round(CL, d=2),' - ',round(CU,d=2),')') ]
        dt <- subset(nsprev.by.age, AGE_LABEL%in%c(20,25,30,35,40,45))	
        dt <- dcast.data.table(dt, LOC_LABEL+AGE_LABEL~SEX_LABEL, value.var='LABEL')
        tmp <- subset(nsprev.ratio.by.loc.age, variable=="PR_FM" & AGE_LABEL%in%c(20,25,30,35,40,45), c(LOC_LABEL, AGE_LABEL, LABEL2))
        dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
        tmp <- subset(nsprev.ratio.by.loc.age, variable=="PR_MF" & AGE_LABEL%in%c(20,25,30,35,40,45), c(LOC_LABEL, AGE_LABEL, LABEL2))
        setnames(tmp, 'LABEL2', 'PR_MF')
        dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))

        filename <- file.path(vl.out.dir, paste0("notsuppamonginfected_round",round,"_220729.csv"))
        fwrite(dt, row.names=FALSE, file=file.path(vl.out.dir,filename))

        #	make table version suppressed
        nsprev.by.age[, LABEL:= paste0(round((1-M)*100, d=1),' (',round((1-CU)*100, d=1),' - ',round((1-CL)*100,d=1),')') ]
        set(nsprev.by.age, NULL, 'SEX_LABEL', nsprev.by.age[, factor(as.character(SEX_LABEL), levels=c('F','M'))])	
        sprev.ratio.by.loc.age[, LABEL2:= paste0(round(M, d=2),' (',round(CL, d=2),' - ',round(CU,d=2),')') ]
        dt <- subset(nsprev.by.age, AGE_LABEL%in%c(20,25,30,35,40,45))	
        dt <- dcast.data.table(dt, LOC_LABEL+AGE_LABEL~SEX_LABEL, value.var='LABEL')
        tmp <- subset(sprev.ratio.by.loc.age,
                      variable=="PR_FM" & AGE_LABEL%in%c(20,25,30,35,40,45), 
                      c(LOC_LABEL, AGE_LABEL, LABEL2))
        dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
        tmp <- subset(sprev.ratio.by.loc.age, 
                      variable=="PR_MF" & AGE_LABEL%in%c(20,25,30,35,40,45),
                      c(LOC_LABEL, AGE_LABEL, LABEL2))
        setnames(tmp, 'LABEL2', 'PR_MF')
        dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))

        filename <- file.path(vl.out.dir, paste0("suppamonginfected_round",round,"_220729.csv"))
        fwrite(dt, row.names=FALSE, file=filename)

        TRUE
    }

    foreach(
            r = vla[, unique(ROUND)],
            .combine='c'
            ) %dopar% {
        cat('Running Round', r, '\n')
        .fit.stan.and.plot.by.round(vla[ ROUND ==r, ], 
                                    iter=10e3, warmup=5e2, chains=1)
    } -> tmp

    return(tmp)
}

vl.suppofpop.by.gender.loc.age.gp<- function(DT, refit=FALSE, vl.out.dir.=vl.out.dir)
{
    cat("\n\n--- Analyse suppression among participants ---\n\n")

    # DT <- copy(dall); refit=FALSE
    DT <- .preprocess.ds.oli(DT)

    tmp <- c('N', 'HIV_N', 'VLNS_N', 'ARV_N')
    vla <- .preprocess.make.vla(DT, select=tmp)

    # Stan file locations
    file.stan <- file.path(path.stan, 'vl_binomial_gp.stan')

    .fit.stan.and.plot.by.round <- function(DT, itr=10e3, wrmp=5e2, chns=1, cntrl = list(max_treedepth= 15, adapt_delta= 0.999))
    {
        #  DT <- copy(vla[ROUND == 16]); refit=FALSE
        round <- DT[, unique(ROUND)]
        stopifnot(length(round) == 1)
        cat('Fitting stan model for round ', round, '\n')
        stan.data <- .make.stan.data.gp(DT, num.var='VLNS_N', den.var='N')

        filename <- paste0( "220729f_suppAmongPop_gp_stan_round",round,".rds") 
        filename <- file.path(vl.out.dir., filename)

        if(file.exists(filename) & refit == FALSE)
        {
            cat('Loading previously run HMC... \n')
            fit <- readRDS(filename)
        }else{
            stan.model <- stan_model(file=file.stan, model_name= 'gp_all')
            fit <- sampling(stan.model, data=stan.data, 
                            iter=itr, warmup=wrmp,
                            chains=chns, control=cntrl)
            saveRDS(fit, file=filename)
        }

        #####################
        # Analyse posterior #
        #####################

        tmp <- grep('f_tilde|alpha|rho_[0-9]+|sex', names(fit), value=TRUE)
        min( summary(fit)$summary[tmp, 'n_eff'] , na.rm=T)	
        re <- rstan::extract(fit)		

        # extract hyperparams rho
        # _______________________

        ps <- c(0.025,0.25,0.5,0.75,0.975)
        tmp <- grep('rho_[0-9]|alpha_[0-9]', names(re), value=TRUE)
        .f <- function(x) transpose(as.data.table(quantile(x, probs=ps)))
        tmp <- lapply(re[tmp], .f)
        tmp <- rbindlist(tmp, idcol='par')
        names(tmp) <- c('GP_hyper_par','CL','IL','M','IU','CU')

        .f <- function(reg, x)
            gsub('^([a-z]+)_([0-9])([0-9])',reg,x)
        tmp <- tmp[ ! GP_hyper_par %like% 'bound']
        tmp[, SEX := as.integer(.f('\\2',GP_hyper_par))]
        tmp[, LOC := as.integer(.f('\\3',GP_hyper_par))]
        tmp[, GP_hyper_par := .f('\\1',GP_hyper_par)]

        tmp1 <- unique(DT[, .(SEX,SEX_LABEL,LOC,LOC_LABEL)])
        nspop.gp.pars <- merge(tmp1, tmp, by=c('SEX','LOC'))
        nspop.gp.pars[, col:='Posterior']
        tmp <- .get.prior.ranges(stan.data, DT=tmp1,
                                 shape=re$rho_hyper_par_shape2[1],
                                 scale=re$rho_hyper_par_scale2[1])

        p <- .plot.gp.hyperparameters(nspop.gp.pars, tmp)
        filename <- paste0('220729f_notsuppAmongPop_gppars_round',round,'.pdf')
        ggsave2(p, file=filename, w=6, h=3)


        # make use stan.data for PPC.
        # ___________________________

        ppDT <- copy(DT)
        cols <- c('M', 'CU', 'CL')

        ppDT[, (cols) := binconf(VLNS_N, N, return.df=T)]

        # make prevalence plot by age
        # ___________________________

        ps <- c(0.025,0.5,0.975)
        tmp <- cbind( apply(re$p_predict_00, 2, quantile, probs=ps),
                     apply(re$p_predict_10, 2, quantile, probs=ps),
                     apply(re$p_predict_01, 2, quantile, probs=ps),
                     apply(re$p_predict_11, 2, quantile, probs=ps)
        )
        rownames(tmp) <- c('CL','M','CU')
        tmp <- as.data.table(reshape2::melt(tmp))
        nspop.by.age <- dcast.data.table(tmp, Var2~Var1, value.var='value')
        tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
        nspop.by.age <- cbind(tmp, nspop.by.age) 
        nspop.by.age <- merge(unique(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), nspop.by.age, by=c('SEX','LOC'))
        nspop.by.age[, STAT:='VLNS']

        p <- ggplot(nspop.by.age, aes(x=AGE_LABEL,colour=SEX_LABEL)) + 		
            geom_ribbon(aes(ymin=CL, ymax=CU, fill=SEX_LABEL), colour=NA, alpha=0.2) +			
            geom_line(aes(x=AGE_LABEL, y=M)) +
            scale_x_continuous( expand=c(0,0) ) + 
            scale_y_continuous(labels=scales:::percent, limits=c(0,.4), expand=c(0,0)) +
            scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
            scale_fill_manual(values=c('M'='royalblue3','F'='deeppink2')) +
            facet_wrap(~LOC_LABEL, ncol=2) +
            theme_bw() +
            theme(legend.position='bottom') +
            labs(x='\nage at visit (years)', 
                 y='population with unsuppressed viral load\n(95% credibility interval)\n', 
                 colour='gender', fill='gender', linetype='location')

        filename <- paste0('220729f_nsuppAmongPop_vs_age_by_gender_fishinland_gp_round_',round,'.pdf')
        ggsave2(p, file=filename, w=6, h=5)

        p <- p + 
            geom_point(data=ppDT, aes(y=M, pch=SEX_LABEL, size=N)) +
            # geom_linerange(data=ppDT, aes(ymin=CL, ymax=CU), linetype='dotted' ) +
            scale_size(range = c(0, 3)) +
            labs(pch='gender', linetype='gender', size='population size')

        filename <- paste0('220729f_nsuppAmongPop_vs_age_by_gender_fishinland_data_gp_round_',round,'.pdf')
        ggsave2(p, file=filename, w=6, h=5)


        # extract basic not supp estimates
        # ________________________________

        DT[, list(N=sum(VLNS_N),
                  P=sum(VLNS_N) / sum(N), 
                  N2=sum(N)-sum(VLNS_N),
                  P2= 1-sum(VLNS_N) / sum(N)
                  ), by=c('LOC_LABEL','SEX_LABEL')]		

        ps <- c(0.025,0.5,0.975)
        tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
        rp <- as.data.table(reshape2::melt( tmp ))
        setnames(rp, 1:3, c('iterations','ROW_ID','P')) 
        tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
        tmp[, ROW_ID:= seq_len(nrow(tmp))]
        tmp <- merge(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL,N)), tmp, by=c('SEX','LOC','AGE_LABEL'), all=TRUE)

        set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N]) 
        set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N])
        set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
        set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
        set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])
        set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])	
        rp <- merge(tmp, rp, by=c('ROW_ID'))	
        rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC_LABEL','SEX_LABEL','iterations')]
        rp <- rp[, list(Q= quantile(P, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','SEX_LABEL')]
        rp <- dcast.data.table(rp, LOC_LABEL+SEX_LABEL~P, value.var='Q')
        rp[, LABEL:= .write.CIs(M, CL, CU, percent=T)]
        nspop.by.sex.loc <- copy(rp)

        # extract risk ratio of suppressed VL female:male and male:female
        # _______________________________________________________________

        tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
        rp <- as.data.table(reshape2::melt( tmp ))
        setnames(rp, 1:3, c('iterations','ROW_ID','P')) 
        tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
        tmp[, ROW_ID:= seq_len(nrow(tmp))]
        tmp <- merge(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL,N)), tmp, by=c('SEX','LOC','AGE_LABEL'), all=TRUE)
        set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N]) 
        set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N])
        set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
        set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
        set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])
        set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])	
        rp <- merge(tmp, rp, by=c('ROW_ID'))	
        rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC_LABEL','SEX_LABEL','iterations')]
        rp <- dcast.data.table(rp, LOC_LABEL+iterations~SEX_LABEL, value.var='P')
        rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC_LABEL','iterations')]
        rp <- melt(rp, id.vars=c('LOC_LABEL','iterations'))
        rp <- rp[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','variable')]
        rp <- dcast.data.table(rp, LOC_LABEL+variable~P, value.var='Q')
        rp[, LABEL:= paste0(sprintf('%2.2f',M),' (',sprintf('%2.2f',CL),' - ',sprintf('%2.2f',CU),')') ]	
        nspop.by.loc <- copy(rp)	

        # extract risk ratio of unsuppressed VL F:M and M:F by age
        # ________________________________________________________

        tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
        rp <- as.data.table(reshape2::melt( tmp ))
        setnames(rp, 1:3, c('iterations','ROW_ID','P')) 
        tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
        tmp[, ROW_ID:= seq_len(nrow(tmp))]
        tmp <- merge(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL,N)), tmp, by=c('SEX','LOC','AGE_LABEL'), all=TRUE)
        set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N]) 
        set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N])
        set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
        set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
        set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])
        set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])	
        rp <- merge(tmp, rp, by=c('ROW_ID'))	
        rp <- dcast.data.table(rp, LOC_LABEL+iterations+AGE_LABEL~SEX_LABEL, value.var='P')
        rp[, PR_FM:= F/M]
        rp[, PR_MF:=M/F]
        rp <- melt(rp, id.vars=c('LOC_LABEL','AGE_LABEL','iterations'), measure.vars=c('PR_FM','PR_MF'))
        rp <- rp[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','AGE_LABEL','variable')]	
        rp <- dcast.data.table(rp, LOC_LABEL+AGE_LABEL+variable~P, value.var='Q')
        rp[, LABEL:= paste0(sprintf('%2.2f',M),' (',sprintf('%2.2f',CL),' - ',sprintf('%2.2f',CU),')') ]	
        nspop.ratio.by.loc.age <- copy(rp)


        # extract if diff in F:M riskratio of unsuppr VL is different in fish vs inland
        # _____________________________________________________________________________

        tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
        rp <- as.data.table(reshape2::melt( tmp ))
        setnames(rp, 1:3, c('iterations','ROW_ID','P')) 
        tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
        tmp[, ROW_ID:= seq_len(nrow(tmp))]
        tmp <- merge(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL,N)), tmp, by=c('SEX','LOC','AGE_LABEL'), all=TRUE)
        set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N]) 
        set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N])
        set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
        set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
        set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])
        set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])	
        rp <- merge(tmp, rp, by=c('ROW_ID'))	
        rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC_LABEL','SEX_LABEL','iterations')]
        rp <- dcast.data.table(rp, LOC_LABEL+iterations~SEX_LABEL, value.var='P')
        rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC_LABEL','iterations')]
        rp <- dcast.data.table(rp, iterations~LOC_LABEL, value.var='PR_FM')
        rp[, PR_FM_D:= fishing-inland]	
        rp[, list(Q= quantile(PR_FM_D, probs=ps), P=c('CL','M','CU'))]

        filename=paste0("220729f_suppAmongPop_round",round,".rda")
        save(DT, re, nspop.by.age,
             nspop.by.sex.loc, nspop.ratio.by.loc.age,
             file=file.path(vl.out.dir.,filename))


        # make table version suppressed
        # _____________________________

        nspop.by.age[, LABEL := .write.CIs(M, CL, CU, percent=T, d=1)]

        set(nspop.by.age, NULL, 'SEX_LABEL', nspop.by.age[, factor(as.character(SEX_LABEL), levels=c('F','M'))])	
        setnames(nspop.ratio.by.loc.age,'LABEL', 'LABEL2')
        dt <- subset(nspop.by.age, AGE_LABEL%in%c(20.5,25.5,30.5,35.5,40.5,45.5))	
        dt <- dcast.data.table(dt, LOC_LABEL+AGE_LABEL~SEX_LABEL, value.var='LABEL')
        tmp <- subset(nspop.ratio.by.loc.age, variable=="PR_FM" & AGE_LABEL%in%c(20.5,25.5,30.5,35.5,40.5,45.5), c(LOC_LABEL, AGE_LABEL, LABEL2))
        dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
        tmp <- subset(nspop.ratio.by.loc.age, variable=="PR_MF" & AGE_LABEL%in%c(20.5,25.5,30.5,35.5,40.5,45.5), c(LOC_LABEL, AGE_LABEL, LABEL2))
        setnames(tmp, 'LABEL2', 'PR_MF')
        dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))

        filename=paste0("220729f_suppAmongPop_round",round,".csv")
        fwrite(dt, row.names=FALSE, file=file.path(vl.out.dir.,filename))

        TRUE
    }

    if(parallelise)
    {
        # Use parallel backend if run locally
        foreach(
                r = vla[, unique(ROUND)],
                .combine='c'
                ) %dopar% {
            cat('Running Round', r, '\n')
            .fit.stan.and.plot.by.round(vla[ ROUND ==r, ])
        } -> tmp

    }else{

        # else do not
        for( r in args$round)
               .fit.stan.and.plot.by.round(vla[ ROUND == r, ]) 
    }

    return(tmp)

}

vl.suppofpop.by.gender.loc.age.icar<- function(DT)
{

    # DT <- copy(dall)
    vl.out.dir <- file.path(vl.out.dir)
    DT <- .preprocess.ds.oli(DT)

    tmp <- c('N', 'HIV_N', 'VLNS_N', 'ARV_N')
    vla <- .preprocess.make.vla(DT, select=tmp)

    file.stan.2 <- file.path(path.stan, 'vl_suppofpop_by_gender_loc_age_icar_1.stan')
    stan.model2 <- stan_model(file.stan.2, model_name= 'icar_age_interactions')


    .fit.stan.and.plot.by.round <- function(DT, iter=20e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999))
    {

        # DT <- copy(DT[ROUND == 16] )
        round <- DT[, unique(ROUND)]
        stopifnot(length(round) == 1)
        cat('Fitting stan model for round ', round, '\n')

        # Get Stan Data
        # _____________

        #	second order RW prior
        node1 <-  c(DT[, seq.int(1, max(AGE)-1L)], DT[, seq.int(1, max(AGE)-2L)])
        node2 <-  c(DT[, seq.int(2, max(AGE))], DT[, seq.int(3, max(AGE))])
        tmp <- sort(node1, index.return=TRUE)$ix

        stan.data <- list(
                          N = nrow(DT),
                          TOTAL = DT[,N],
                          K = DT[,VLNS_N],
                          AGE_N = DT[, max(AGE)],
                          AGE = DT[, AGE],
                          SEX = DT[, SEX],
                          LOC = DT[, LOC],
                          node1 = node1[tmp],
                          node2 = node2[tmp],
                          N_edges =  length(stan.data$node1)
        )

        fit <- sampling(stan.model2, data=stan.data, iter=iter, warmup=warmup, chains=chains, control=control)
        filename <-  paste0("notsuppAmongPop_icar_stan_220729_round_",round,".rds")
        save(fit, file=file.path(vl.out.dir,filename))

        min( summary(fit)$summary[, 'n_eff'] )	
        re <- rstan::extract(fit)	
        ps <- c(0.025,0.5,0.975)


        # make prevalence plot by age
        # ___________________________

        tmp <- apply(re$p, 2, quantile, probs=ps)
        rownames(tmp) <- c('CL','M','CU')
        tmp <- as.data.table(reshape2::melt(tmp))	
        nspop.by.age <- cbind(DT, dcast.data.table(tmp, Var2~Var1, value.var='value'))	
        p <- ggplot(nspop.by.age) + 		
            geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2) +			
            geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL)) +
            scale_x_continuous( expand=c(0,0) ) + 
            scale_y_continuous(label=scales:::percent) +
            scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
            facet_wrap(~LOC_LABEL, ncol=2) +
            theme_bw() +
            labs(x='\nage at visit (years)', 
                 y='individuals with unsuppressed viral load\n(95% credibility interval)\n', 
                 colour='gender')

        filename <- paste0('220729_notsuppAmongPop_vs_age_by_gender_fishinland_round',round,'.pdf')
        ggsave2(p, filename, w=6, h=5)


        # extract basic not supp estimates
        # ________________________________

        DT[, list(N=sum(VLNS_N), P=sum(VLNS_N) / sum(N), N2=sum(N)-sum(VLNS_N), P2= 1-sum(VLNS_N) / sum(N)), by=c('LOC_LABEL','SEX_LABEL')]		


        rp <- as.data.table(reshape2::melt( re$p ))
        setnames(rp, 2:3, c('ROW_ID','P')) 
        rp <- merge(rp, DT, by='ROW_ID')
        rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
        rp <- rp[, list(Q= quantile(P, probs=ps), P=c('CL','M','CU')), by=c('LOC','SEX')]
        rp <- dcast.data.table(rp, LOC+SEX~P, value.var='Q')
        rp[, LABEL:= paste0(round(M, d=2),'% (',round(CL, d=2),'% - ',round(CU,d=2),'%)') ]
        rp <- merge(unique(subset(DT, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
        nspop.by.sex.loc <- copy(rp)


        # extract risk ratio of unsuppressed VL F:M and M:F
        # _________________________________________________

        rp <- as.data.table(reshape2::melt( re$p ))
        setnames(rp, 2:3, c('ROW_ID','P')) 
        rp <- merge(rp, DT, by='ROW_ID')
        rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
        rp <- merge(unique(subset(DT, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
        rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations~SEX_LABEL, value.var='P')
        rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC','iterations')]
        rp <- melt(rp, id.vars=c('LOC','iterations'))[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC','variable')]
        rp <- dcast.data.table(rp, LOC+variable~P, value.var='Q')
        rp[, LABEL:= paste0(round(M, d=2),' (',round(CL, d=2),' - ',round(CU,d=2),')') ]
        rp <- merge(unique(subset(DT, select=c(LOC,LOC_LABEL))), rp, by=c('LOC'))
        nspop.ratio.by.loc <- copy(rp)


        # extract if diff in F:M risk ratio of unsupp VL is different in fish vs inland
        # _____________________________________________________________________________

        rp <- as.data.table(reshape2::melt( re$p ))
        setnames(rp, 2:3, c('ROW_ID','P')) 
        rp <- merge(rp, DT, by='ROW_ID')
        rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
        rp <- merge(unique(subset(DT, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
        rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations~SEX_LABEL, value.var='P')
        rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC_LABEL','iterations')]
        rp <- dcast.data.table(rp, iterations~LOC_LABEL, value.var='PR_FM')
        rp[, PR_FM_D:= fishing-inland]	
        rp[, list(Q= quantile(PR_FM_D, probs=ps), P=c('CL','M','CU'))]

        # extract risk ratio of unsuppressed VL female:male and male:female by age
        # ________________________________________________________________________

        rp <- as.data.table(reshape2::melt( re$p ))
        setnames(rp, 2:3, c('ROW_ID','P')) 
        rp <- merge(rp, DT, by='ROW_ID')
        rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations+AGE+AGE_LABEL~SEX_LABEL, value.var='P')
        rp[, PR_FM:= F/M]
        rp[, PR_MF:=M/F]
        rp <- melt(rp, id.vars=c('LOC_LABEL','AGE_LABEL','iterations'), measure.vars=c('PR_FM','PR_MF'))[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','AGE_LABEL','variable')]
        rp <- dcast.data.table(rp, LOC_LABEL+AGE_LABEL+variable~P, value.var='Q')
        nspop.ratio.by.loc.age <- copy(rp)


        filename <- paste0("notsuppAmongPop_220729_round", round, ".rda")
        cat('Saving', filename, '...\n')
        save(DT, re, nspop.by.age, nspop.ratio.by.loc, nspop.ratio.by.loc.age,
             file=file.path(vl.out.dir,filename))

        # make table
        # __________

        nspop.by.age[, LABEL:= paste0(round(M*100, d=1),' (',round(CL*100, d=1),' - ',round(CU*100,d=1),')') ]
        set(nspop.by.age, NULL, 'SEX_LABEL', nspop.by.age[, factor(as.character(SEX_LABEL), levels=c('F','M'))])	
        nspop.ratio.by.loc.age[, LABEL2:= paste0(round(M, d=2),' (',round(CL, d=2),' - ',round(CU,d=2),')') ]
        dt <- subset(nspop.by.age, AGE_LABEL%in%c(20,25,30,35,40,45))	
        dt <- dcast.data.table(dt, LOC_LABEL+AGE_LABEL~SEX_LABEL, value.var='LABEL')
        tmp <- subset(nspop.ratio.by.loc.age, variable=="PR_FM" & AGE_LABEL%in%c(20,25,30,35,40,45), c(LOC_LABEL, AGE_LABEL, LABEL2))
        dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
        tmp <- subset(nspop.ratio.by.loc.age, variable=="PR_MF" & AGE_LABEL%in%c(20,25,30,35,40,45), c(LOC_LABEL, AGE_LABEL, LABEL2))
        setnames(tmp, 'LABEL2', 'PR_MF')
        dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))

        filename <- paste0( "notsuppAmongPop_220729_round",round,".csv")
        fwrite(dt, row.names=FALSE, file=file.path(vl.out.dir,filename))

        TRUE
    }

    foreach(
            r = vla[, unique(ROUND)],
            .combine='c'
            ) %dopar% {
        cat('Running Round', r, '\n')
        .fit.stan.and.plot.by.round(vla[ ROUND ==r, ], iter=5e3)
    } -> tmp

    return(tmp)
}

vl.vlrunningprops.by.gender.loc.age.round <- function(DT)
{
    # DT <- copy(dall)
    vl.out.dir <- file.path(vl.out.dir)
    DT <- .preprocess.ds.oli(DT)
    ans <- .preprocess.make.vla.2(DT)

    # HIV prevalence
    # ______________

    p <- ggplot(ans) + 		
        geom_ribbon(aes(x=AGEYRS, ymin=PHIV_CL, ymax=PHIV_CU, fill=SEX), alpha=0.2) +
        geom_line(aes(x=AGEYRS, y=PHIV_MEAN, colour=SEX)) +
        scale_x_continuous( expand=c(0,0) ) + 
        scale_y_continuous(labels=scales:::percent, expand=c(0,0)) +
        scale_colour_manual(values=c('men'='royalblue3','women'='deeppink2')) +
        scale_fill_manual(values=c('men'='royalblue3','women'='deeppink2')) +
        facet_grid(ROUND~FC) +
        theme_bw() +
        theme(legend.position='bottom') +
        labs(x='\nage at visit (years)', 
             y='HIV prevalence (95% CI)\n', 
             colour='gender', fill='gender', linetype='location',
             title='Rolling average prevalence estimates'
        )

    filename <- '220729_hivprevalence_vs_age_by_gender_fishinland_rounds.pdf'
	ggsave2(p, file=filename, w=6, h=10)
	
	# HIV unsuppressed viral load
	# ___________________________

	p <- ggplot(ans) + 		
        geom_ribbon(aes(x=AGEYRS, ymin=PVLNS_CL, ymax=PVLNS_CU, fill=SEX), alpha=0.2) +			
        geom_line(aes(x=AGEYRS, y=PVLNS_MEAN, colour=SEX)) +
        scale_x_continuous( expand=c(0,0) ) + 
        scale_y_continuous(labels=scales:::percent) +
        scale_colour_manual(values=c('men'='royalblue3','women'='deeppink2')) +
        scale_fill_manual(values=c('men'='royalblue3','women'='deeppink2')) +
        facet_grid(ROUND~FC) +
        theme_bw() +
        theme(legend.position='bottom') +
        labs(title='HIV prevalence: rolling average estimates',
             x='\nage at visit (years)', 
             y='proportion unsuppressed HIV (95% CI)\n', 
             colour='gender', fill='gender', linetype='location')

    filename <- '220729_hivnotsupp_vs_age_by_gender_fishinland_rounds.pdf'
    ggsave2(p, file=filename, w=6, h=10)


    # HIV unsuppressed viral load among HIV+
    # ______________________________________

    p <- ggplot(ans) + 		
        geom_ribbon(aes(x=AGEYRS, ymin=PVLNSofHIV_CL, ymax=PVLNSofHIV_CU, fill=SEX), alpha=0.2) +
        # geom_line(aes(x=AGEYRS, y=PARVofHIV_MEAN, colour=SEX), linetype='dotted') +
        geom_line(aes(x=AGEYRS, y=PVLNSofHIV_MEAN, colour=SEX)) +
        scale_x_continuous( expand=c(0,0) ) + 
        scale_y_continuous(labels=scales:::percent, expand=c(0,0)) +
        scale_colour_manual(values=c('men'='royalblue3','women'='deeppink2')) +
        scale_fill_manual(values=c('men'='royalblue3','women'='deeppink2')) +
        facet_grid(ROUND~FC) +
        theme_bw() +
        theme(legend.position='bottom') +
        labs(title='Unsuppressed among HIV+: rolling average estimate',
             x='\nage at visit (years)', 
             y='proportion unsuppressed HIV among infected (95% CI)\n', 
             colour='gender', fill='gender',
             linetype='location')

    filename <- '220729_hivnotsuppofhiv_vs_age_by_gender_fishinland_rounds.pdf'
    ggsave2(p, file=filename, w=6, h=10)


    # write results to file
    # _____________________

    setkey(ans, ROUND, FC, SEX, AGEYRS)     
    ans[, PHIV_L:= .write.CIs(PHIV_MEAN, PHIV_CL, PHIV_CU, percent=T, d=1)]
    ans[, PVLNS_L:= .write.CIs(PVLNS_MEAN, PVLNS_CL, PVLNS_CU, percent=T, d=1)]
    ans[, PVLNSofHIV_L:= .write.CIs(PVLNSofHIV_MEAN, PVLNSofHIV_CL, PVLNSofHIV_CU, percent=T, d=1)]

    fwrite(ans, file=file.path(vl.out.dir,'220729_keystats_by_age_gender_fishinland_rounds.csv'))
}

vl.vlrunningmean.by.gender.loc.age.round <- function(DT)
{
    # DT <- copy(dall)
    vl.out.dir <- file.path(vl.out.dir)
    DT <- .preprocess.ds.oli(DT)

    tmp <- seq.int(min(DT$AGEYRS), max(DT$AGEYRS))
    tmp1 <- DT[, sort(unique(ROUND))]

    vla <- as.data.table(expand.grid(ROUND=tmp1,
                                     FC=c('fishing','inland'),
                                     SEX=c('M','F'),
                                     AGEYRS=tmp))

    ans <- vla[, {
        z <- which(DT$ROUND==ROUND &DT$FC==FC & DT$SEX==SEX & DT$AGEYRS<=(AGEYRS+2) & DT$AGEYRS>=(AGEYRS-2))
        z2 <- mean( DT$VLC[z] )
        z3 <- sd(DT$VLC[z])/sqrt(length(z))
        list(N= length(z),
             VLCM_M= z2,
             VLCM_CL= z2-1.96*z3,
             VLCM_CU= z2+1.96*z3
        )
    }, by=names(vla)]
    set(ans, NULL, 'SEX', ans[, factor(SEX, levels=c('M','F'), labels=c('men','women'))])


    p <- ggplot(ans) + 
        #geom_errorbar(aes(x=AGEYRS, ymin=VLCM_CL, ymax=VLCM_CU)) +		
        geom_ribbon(aes(x=AGEYRS, ymin=VLCM_CL, ymax=VLCM_CU, fill=SEX), alpha=0.2) +
        # geom_hline(yintercept=1e3) +
        geom_line(aes(x=AGEYRS, y=VLCM_M, colour=SEX)) +
        scale_x_continuous( expand=c(0,0) ) + 
        scale_y_continuous() +
        scale_colour_manual(values=c('men'='royalblue3','women'='deeppink2')) +
        scale_fill_manual(values=c('men'='royalblue3','women'='deeppink2')) +
        facet_grid(ROUND~FC, scales='free_y') +
        theme_bw() +
        theme(legend.position='bottom') +
        labs(x='\nage at visit (years)', 
             y='mean viral load (95% CI)\n', 
             colour='gender', fill='gender' ,linetype='location')
    p

    filename=paste0('220729_vlmean_vs_age_by_gender_fishinland_rounds.pdf')
	ggsave2(p, file=filename, w=8, h=5)
	
	
	# sqrt transformation did not work, gave too low means
	DT[, VLCS:= sqrt(VLC)]
	vlclo <- DT[, loess(VLCS ~ AGEYRS, control=loess.control(trace.hat='approximate'))]	
	ans <- subset(DT, select=c(FC, SEX, AGEYRS, VLC, VLCS))	
	ans[, VLCLO_M:= (vlclo$fitted)^2]
	p <- ggplot(ans) + 
        geom_line(aes(x=AGEYRS, y=VLCLO_M)) +
        scale_x_continuous( expand=c(0,0) )	
    p

	# loess mean below 0 for some age groups, not a good model
    ans <- DT[, {
        vlclo <- loess(VLC ~ AGEYRS, control=loess.control(trace.hat='approximate'))
        list(	VLC= VLC,
             AGEYRS= AGEYRS,
             VLCLO_M= vlclo$fitted 
        )				
    }, by=c('ROUND','FC','SEX')]	
    ans
    predict(vlclo, newdata=NULL, se=TRUE)

    # What about log10?
    # check distributions etc...
    ans <- DT[, {

        z <- which(DT$ROUND==ROUND &DT$FC==FC & DT$SEX==SEX & DT$AGEYRS<=(AGEYRS+2) & DT$AGEYRS>=(AGEYRS-2))
        z1 <- pmax(log(DT$VLC[z], 10), 0)
        z2 <- mean( z1 )
        z3 <- sd( z1 )/sqrt(length(z))

        list(N= length(z),
             z1= z1,
             VLCM_M= z2,
             VLCM_CL= z2-1.96*z3,
             VLCM_CU= z2+1.96*z3
        )
    }, by=c('AGEYRS','ROUND', 'FC', 'SEX')]
    ans
}

vl.vldistribution.by.gender.loc<- function()
{
    # TODO: add facets per round?

    # DT <- copy(dall)
    DT <- .preprocess.ds.oli(DT)

    tmp <- seq.int(min(DT$AGEYRS), max(DT$AGEYRS))
    tmp1 <- DT[, sort(unique(ROUND))]

    # plot proportion of population with viral load > x
    x <- seq(log(1),log(max(DT$VLC)), length.out=1e3)
    x <- c(0,exp(x))

    vld <- as.data.table(expand.grid(ROUND=tmp1,
                                     X=x,
                                     SEX=c('M','F'),
                                     FC=c('fishing','inland'),
                                     HIV_AND_VLD=c(0,1)))

    ans <- vld[, {
        n <- sum(DT$ROUND==ROUND, DT$SEX==SEX & DT$FC==FC & DT$HIV_AND_VLD>=HIV_AND_VLD)
        k <- sum(DT$ROUND==ROUND, DT$SEX==SEX & DT$FC==FC & DT$HIV_AND_VLD>=HIV_AND_VLD & X<DT$VLC)
        z<- as.vector( unname( binconf(k, n) ) )				
        list(N=n, K=k, P_M= z[1], P_CL=z[2], P_CU=z[3] )
    }, by=names(vld)]

    set(ans, NULL, 'HIV_AND_VLD', factor(ans[,HIV_AND_VLD], levels=c(0,1), labels=c('all study participants','infected study participants\nwith detectable viral load')))
    set(ans, NULL, 'SEX', ans[, factor(SEX, levels=c('M','F'), labels=c('men','women'))])

    ans <- subset(ans, !(HIV_AND_VLD=='infected study participants\nwith detectable viral load' & X<VL_DETECTABLE) )
    ans <- subset(ans, !(HIV_AND_VLD=='all study participants' & X<VL_DETECTABLE) )

    p <- ggplot(ans) +
        geom_line(aes(x=X, y=P_M, group=interaction(FC,SEX), colour=SEX, linetype=FC)) +
        scale_x_log10() +
        scale_y_continuous(labels=scales:::percent, expand=c(0,0)) +
        scale_colour_manual(values=c('men'='royalblue3','women'='deeppink2')) +
        geom_text(aes(x=1e3, y=P_M * 1.03, label="")) +
        theme_bw() +
        facet_wrap(~HIV_AND_VLD, scales='free', ncol=2) +
        labs(	x='\nviral load\n(copies / ml)', 
             y='proportion of individuals with larger viral load\n',
             colour='gender', linetype='location')

    # TODO: add `_round?`
    filename <- '220729_vldistribution_by_gender_fishinland.pdf'
	ggsave2(p, file=filename, w=9, h=5)
}

vl.keystats.by.gender.loc<- function(DT)
{

        #TODO: by_round as well?

        # DT <- copy(dall)
	vl.out.dir <- file.path(vl.out.dir)
        DT <- .preprocess.ds.oli(DT)
	
	# entire population: 
	# mean viral load, proportion with DVL
	DT[, mean(VLC)]
	# 2290.494
	binconf( length(which(DT$VLD==1)), nrow(DT) )
	# PointEst   Lower      Upper
	# 0.05717964 0.05393704 0.06060469

	# stratified by men/women inland/fishing
	# ______________________________________

        .f <- function(x,y)
        {
                as.vector( unname ( binconf(sum(x), length(y) )  ) )
        }

	ans <- DT[, {
			z  <- .f(HIV_STATUS==1, HIV_STATUS)
			z2 <- .f(VLNS==1, VLNS)
			z3 <- .f(VLNS==1, which(HIV_STATUS==1))

			list(N= length(HIV_STATUS),
				 PHIV_MEAN= z[1],
				 PHIV_CL= z[2],
				 PHIV_CU= z[3],				 
				 PVLNS_MEAN= z2[1],
				 PVLNS_CL= z2[2],
				 PVLNS_CU= z2[3],
				 PVLNSofHIV_MEAN= z3[1],
				 PVLNSofHIV_CL= z3[2],
				 PVLNSofHIV_CU= z3[3],				 
				 VLC_MEAN= mean(VLC))	
			}, by='SEX']

	ans[, FC:='overall']	
	tmp <- DT[, {
				z  <- .f(HIV_STATUS==1, HIV_STATUS)
				z2 <- .f(VLNS==1, VLNS)
				z3 <- .f(VLNS==1, which(HIV_STATUS==1))

				list(N= length(HIV_STATUS),
						PHIV_MEAN= z[1],
						PHIV_CL= z[2],
						PHIV_CU= z[3],				 
						PVLNS_MEAN= z2[1],
						PVLNS_CL= z2[2],
						PVLNS_CU= z2[3],
						PVLNSofHIV_MEAN= z3[1],
						PVLNSofHIV_CL= z3[2],
						PVLNSofHIV_CU= z3[3],				 						
						VLC_MEAN= mean(VLC))	
			}, by=c('FC','SEX')]

	ans <- rbind(tmp, ans)

	set(ans, NULL, 'FC', factor(ans$FC, levels=c('overall','fishing','inland')))
	set(ans, NULL, 'SEX', factor(ans$SEX, levels=c('F','M')))
	setkey(ans, FC, SEX)
	
        .f <- function(m, l, u)
        {
                .r <- function(x) round(x*100, digits=1)
                paste0( .r(m),' [', .r(l),'-', .r(u) ,']' )
        }

	ans[, PHIV_L:= .f( PHIV_MEAN, PHIV_CL, PHIV_CU) ]
	ans[, PVLNS_L:= .f( PVLNS_MEAN, PVLNS_CL, PVLNS_CU) ]
	ans[, PVLNSofHIV_L:= .f( PVLNSofHIV_MEAN, PVLNSofHIV_CL, PVLNSofHIV_CU)]
	
        # TODO: change name?
	fwrite(ans, file=file.path(vl.out.dir, '220729_keystats_by_gender_fishinland.csv'))
}

vl.vlratio.by.loc<- function(DT)
{

    # DT <- copy(dall)
	vl.out.dir <- file.path(vl.out.dir)
    DT <- .preprocess.ds.oli(DT)
	
	ans	<- as.data.table(expand.grid(BS= 1:1e3, FC=c('fishing','inland')))
	set.seed(42)

	ans <- ans[, {
				zm <- which(DT$FC==FC & DT$SEX=='M')
				zf <- which(DT$FC==FC & DT$SEX=='F')
				zm <- sample(zm, length(zm), replace=TRUE)
				zf <- sample(zf, length(zf), replace=TRUE)
				list(VLCM_M=mean(DT$VLC[zm]), VLCM_F=mean(DT$VLC[zf])) 
			}, by=c('FC','BS')]

	ans[, VLCR:= VLCM_M/VLCM_F]	

	ans <- ans[, list( 
                          V = quantile(VLCR, prob=c(0.5, 0.025, 0.975)),
                          P = c('M','CL','CU')
                          ), by=c('FC')]
	ans <- dcast.data.table(ans, FC~P, value.var='V')
	
	#
	# stratified by men/women inland/fishing
	# 
    .f <- function(x,y)
        as.vector( unname ( binconf(sum(x), length(y) )  ) )

	ans <- DT[, {
        z  <- .f(HIV_STATUS==1, HIV_STATUS) 
        z2 <- .f(VLNS==1, VLNS) 
        z3 <- .f(VLNS==1, which(HIV_STATUS==1))

        list( N =length(HIV_STATUS),
             PHIV_MEAN= z[1],
             PHIV_CL= z[2],
             PHIV_CU= z[3],				 
             PVLNS_MEAN= z2[1],
             PVLNS_CL= z2[2],
             PVLNS_CU= z2[3],
             PVLNSofHIV_MEAN= z3[1],
             PVLNSofHIV_CL= z3[2],
             VLC_MEAN= mean(VLC) )	
    }, by='SEX']

    ans[, FC:='overall']	

    tmp <- DT[, {
        z  <- .f(HIV_STATUS==1, HIV_STATUS) 
        z2 <- .f(VLNS==1, VLNS) 
        z3 <- .f(VLNS==1, which(HIV_STATUS==1)) 

        list(N= length(HIV_STATUS),
             PHIV_MEAN= z[1],
             PHIV_CL= z[2],
             PHIV_CU= z[3],				 
             PVLNS_MEAN= z2[1],
             PVLNS_CL= z2[2],
             PVLNS_CU= z2[3],
             PVLNSofHIV_MEAN= z3[1],
             PVLNSofHIV_CL= z3[2],
             PVLNSofHIV_CU= z3[3],				 		
             VLC_MEAN= mean(VLC))	
    }, by=c('FC','SEX')]

    ans <- rbind(tmp, ans)
    set(ans, NULL, 'FC', factor(ans$FC, levels=c('overall','fishing','inland')))
    set(ans, NULL, 'SEX', factor(ans$SEX, levels=c('F','M')))
    setkey(ans, FC, SEX)
}

vl.age.gender<- function()
{
    require(data.table)
    prjdir	<- '~/Box Sync/OR_Work/2018/2018_RakaiViralLoad'
    infile	<- file.path(prjdir,'data','191101_data_round17_vl_gps.rda')
    load( infile )

    # drop few infecteDT with missing VL
    DT <- subset(DT, HIV_STATUS==0 | (HIV_STATUS==1 & !is.na(VL_COPIES)) )
	# set VL for uninfected to 0, and VL with undetectable VL to 0
    set(DT, DT[, which(HIV_STATUS==0)], 'VL_COPIES', 0)
    set(DT, DT[, which(HIV_STATUS==1 & VL_UNDETECTABLE==1)], 'VL_COPIES', 0)

    # 
    # calculate proportion with VL > x among participants

    # do general by as characters
    # then determine sort index
    # then calculate empirical quantile
    DT <- DT[order(SEX,VL_COPIES),]
    DT[VL]

    DT[, sort(unique(VL_COPIES))]
    #dv <- data.table(VL:= )
}

vl.get.data.round17<- function()
{
    require(data.table)
    prjdir	<- '~/Box/OR_Work/2018/2018_RakaiViralLoad'
    infile	<- 'data_raw/ViralLoad_Data_Pangea_Ratmann.rda'
    load( file.path(prjdir, infile) )

    # subset to survey round 17
    # _________________________

    ds		<- subset(as.data.table(survey_data), visit==17)
    # reset dates from Date format to numeric
    for(x in c('visit_date','lastNegDate','firstPosDate'))
    {
        set(ds, NULL, x, date2numeric(ds[[x]]))
    }
    # make all column names upper case
    setnames(ds, colnames(ds), toupper(colnames(ds)))
    # define FISHING_COMM
    ds[, FC:= as.character(factor(COMM_NUM%in%c(770,771,774,38),levels=c(TRUE,FALSE),labels=c('fishing','inland')))]
    # define ARVMED
    set(ds, ds[, which(ARVMED==8)], 'ARVMED', NA_integer_)
    set(ds, NULL, 'ARVMED', ds[, as.integer(as.character(factor(ARVMED, levels=c(1,2), labels=c('1','0'))))])

    #
    # prepare GPS coordinates
    #
    dg	<- as.data.table(gpsdat)	
    # bring dates into common format
    setnames(dg, colnames(dg), gsub('\\.','_',toupper(colnames(dg))))
    tmp	<- which(dg[, grepl('([0-9]+)/([0-9]+)/([0-9]+)',GPS_DATE)])
    set(dg, tmp, 'GPS_DATE', dg[tmp,gsub('([0-9]+)/([0-9]+)/([0-9]+)','\\3-\\1-\\2',GPS_DATE)])
    tmp	<- which(dg[, grepl('([0-9]+)-([A-Za-z]+)-([0-9]+)',GPS_DATE)])
    set(dg, tmp, 'GPS_DATE', dg[tmp,gsub('([0-9]+)-([A-Za-z]+)-([0-9]+)','20\\3-\\2-\\1',GPS_DATE)])	
    set(dg, NULL, 'GPS_DATE', dg[,gsub('Nov','11',gsub('Oct','10',gsub('Sep','09',gsub('Aug','08',gsub('July','07',gsub('Jun','06',gsub('May','05',GPS_DATE)))))))])
    # reset dates from character format to numeric
    set(dg, NULL, 'GPS_DATE', date2numeric(dg[,GPS_DATE]))
    # make households per date unique
    dg	<- unique(dg, by=c('HHID','GPS_DATE'))

    #
    # add to surveyed individuals the GPS of their households	
    # 
    tmp	<- unique(subset(ds, select=c(RCCS_STUDYID, VISIT_DATE, HHID)))
    tmp	<- merge(tmp, dg, by='HHID', all.x=TRUE)
    # some households do not have GPS coordinates
    ch	<- subset(tmp, is.na(LATITUDE_JITTER) | is.na(LATITUDE_JITTER))
    if(nrow(ch))
    {
        cat('\nNumber of households without GPS coordinates, n=', nrow(ch))
        write.csv(ch, file=file.path(prjdir,'data/check_missing_coordinates.csv'))
        #	521 households without GPS coordinates
    }
    # for every individual, extract house closest in time
    tmp		<- subset(tmp, !is.na(LATITUDE_JITTER) & !is.na(LATITUDE_JITTER))
    tmp2	<- tmp[, list(GPS_DATE= GPS_DATE[which.min(abs(GPS_DATE-VISIT_DATE))[1]]), by=c('RCCS_STUDYID','VISIT_DATE')]
    tmp		<- merge(tmp, tmp2, by=c('RCCS_STUDYID','VISIT_DATE','GPS_DATE'))
    stopifnot(nrow(tmp)==nrow(tmp2))
    set(tmp, NULL, c('COMM','HOUSE'), NULL)	
    ds		<- merge(ds, tmp, by=c('RCCS_STUDYID','VISIT_DATE','HHID'), all.x=TRUE)

    #
    # extract viral loads from round 17
    #
    dvl		<- subset(as.data.table(viralLoads), visit==17)
    setnames(dvl, colnames(dvl), toupper(colnames(dvl)))
    setnames(dvl, c('DATE','COPIES','DONEBY'), c('VL_DATE','VL_COPIES','VL_DONEBY'))
    set(dvl, NULL, 'VL_DATE', date2numeric(dvl[,VL_DATE]))
    stopifnot( !nrow(subset(dvl, is.na(VL_COPIES))) )
    # check if one viral load measurement per person
    tmp <- dvl[, list(N_VL=length(VL_COPIES)), by='RCCS_STUDYID']
    stopifnot( !nrow(subset(tmp, N_VL>1)) )	
    # merge with main data
    set(dvl, NULL, 'VISIT', dvl[, as.numeric(VISIT)])
    set(dvl, NULL, 'VL_DONEBY', NULL)
    ds	<- merge(ds, dvl, by=c('RCCS_STUDYID','VISIT'), all.x=TRUE)

    # check if viral load for all infected
    ch	<- subset(ds, HIV_STATUS==1 & is.na(VL_COPIES))
    if(nrow(ch))
    {
        cat('\nFound infected individuals without VL measurement, n=',nrow(ch))
        write.csv(ch, file=file.path(prjdir,'data/check_missing_viralloads.csv'))
        # 13 HIV+ individuals without VL
    }
    ds[, HIV_AND_VL:= as.integer(HIV_STATUS==1 & !is.na(VL_COPIES))]


    save(ds, file=file.path(prjdir,'data','191101_data_round17_vl_gps.rda'))
}

prop.dectectable.viraemia<- function()
{
    require(data.table)
    require(rgdal)
    require(rgeos)
    library(raster)
    require(RColorBrewer) #Map colours

    # load data
    infile	<- '~/Box Sync/OR_Work/2018/2018_RakaiViralLoad/data/merged_round17_vl_gps.rda'
    load(infile)

    tmp <- ds[, list(
                     HIV_POS=  sum(HIV_STATUS==1), 
                     HIV_NEG=  sum(HIV_STATUS==0),  
                     HIV_PREV= sum(HIV_STATUS==1)/length(HIV_STATUS)
                     ), by='COMM_NUM']

    thr	<- 1e3		
    tmp2 <- ds[HIV_STATUS==1]
    tmp2 <- tmp2[, list(
                        VL_D=  sum(VL_COPIES>thr), 
                        VL_U=  sum(VL_COPIES<=thr),  
                        VL_DP= sum(VL_COPIES>thr)/length(VL_COPIES)
                        ), by='COMM_NUM']
    tmp	<- merge(tmp, tmp2, by='COMM_NUM')
    tmp[, POP_VL_DP:= HIV_PREV*VL_DP]
    ggplot(tmp, aes(y=COMM_NUM, x=POP_VL_DP)) + geom_point()

    tmp3	<- subset(ds, HIV_STATUS==1 & COMM_NUM==38)
    ggplot(tmp3, aes(x=VL_COPIES)) + geom_histogram() + facet_grid(~ARVMED)

    tmp3	<- subset(ds, HIV_STATUS==1 & COMM_NUM==38 & ARVMED==2)
    tmp3[, VL_COPIES_C:= cut(VL_COPIES, breaks=c(0,1,10,100,1000,1e4,1e5,1e6,1e7,1e10), right=FALSE)]
    tmp3[, table(VL_COPIES_C)]

    tmp3	<- subset(ds, HIV_STATUS==1 & COMM_NUM==38 & ARVMED==1)
    tmp3[, VL_COPIES_C:= cut(VL_COPIES, breaks=c(0,1,10,100,1000,1e4,1e5,1e6,1e7,1e10), right=FALSE)]
    tmp3[, table(VL_COPIES_C)]

    tmp3	<- subset(ds, HIV_STATUS==1 & COMM_NUM==38)
    tmp3[, table(VL_COPIES>1)]
}

make.map.190129	<- function(DT)
{
    require(data.table)
    require(rgdal)
    require(rgeos)
    library(raster)
    require(RColorBrewer) #Map colours

    # DT <- copy(dall)
    DT <- .preprocess.ds.oli(DT)
    infile <- file.path(indir.deepsequence.data, 'RCCS_R15_R18', 'Rakai_community_geography_R15.rda')

    # Get longitude and latitude
    tmp <- new.env()
    load(infile, envir=tmp)
    ds <- as.data.table(tmp$comgps)
    cols <- c('latitude', 'longitude')
    ds[, (cols):=lapply(.SD, unlist), .SDcols=cols]
    ds[, COMM_NUM := as.integer(COMM_NUM)]
    setnames(ds, cols, paste0(toupper(cols), '_JITTER'))

    # merge with data
    DT <- merge(DT, ds, all.x=T, by='COMM_NUM')
    DT[is.na(LATITUDE_JITTER), uniqueN(COMM_NUM)] -> tmp
    DT[, cat('missing geoloc for', tmp, 'communities\n')]

    #convert the data into a data table
    setnames(DT, 'VLU', 'VL_UNDETECTABLE')
    DT <- DT[,.(STUDY_ID, SEX, AGEYRS, HIV_STATUS, LATITUDE_JITTER, LONGITUDE_JITTER, VL_COPIES, VL_UNDETECTABLE)]
    #set the NA VL to 0
    DT[is.na(VL_COPIES), VL_COPIES:=0]
    DT[,VL_DETECTABLE := as.numeric(VL_COPIES>=1000)]
    DT[,RCCS_STUDYID2:= seq_len(.N) ]

    ##############################
    # Load in Uganda Shape files #
    ##############################

    uganda1<-raster::getData('GADM',country="UGA",level=1)# Admin unit 1
    uganda3<- raster::getData('GADM', country='UGA', level=3)

    rakai1<-subset(uganda1, NAME_1=="Rakai")
    rakai3<- subset(uganda3, NAME_1=="Rakai")
    masaka1<-subset(uganda1, NAME_1=="Masaka")
    # Create a smaller Rakai for plotting (not current Rakai region no longer includes kabula subdistrict 3)
    #minirak<-rakai3[which(rakai3$NAME_2!="Kabula" | rakai3$NAME_3=="Lyantonde Tc" | rakai3$NAME_3=="Lyantonde"),]
    minirak<-rakai3[which(rakai3$NAME_2!="Kabula" | rakai3$NAME_3=="Lyantonde Tc"),]

    ##################################
    #set the coordinates of the data #
    ##################################

    DT <- DT[! is.na(LATITUDE_JITTER) & ! is.na(LONGITUDE_JITTER)]
    coordinates(DT) <- ~ LONGITUDE_JITTER+LATITUDE_JITTER
    #set coordinate system to match uganda files
    proj4string(DT) <- proj4string(uganda1)

    #convert to m in order to build a 30x30m grid
    newcrs <- CRS("+proj=robin +datum=WGS84")
    DTnew  <- spTransform(DT, newcrs)
    rakai1trans  <- spTransform(rakai1, newcrs)
    miniraktrans <- spTransform(minirak, newcrs)
    masaka1trans <- spTransform(masaka1, newcrs)

    #######################################
    #Combine rakai1trans and masaka1trans #
    #######################################

    outline<- union(rakai1trans, masaka1trans)
    #find the extent of the data
    exnew<- extent(DTnew)
    #extent of the maps
    exmap<- extent(outline)

    #chose extent to cover all the data and rakai district

    #With a 30m grid, I think the same individuals are usually entering calculations for a large number of grid points
    #Do we really need a 30m grid? Why not 100m?

    grid<- raster(xmn=min(exnew[1], exmap[1]), xmx= exnew[2], ymn=exmap[3], ymx=exnew[4], res=100 )
    grid[]<- 1:ncell(grid) #No longer needed

    # set the coordinate reference system to match
    proj4string(grid)<- proj4string(DTnew) 

    #restrict grid to map
    gridmask<- mask(grid, outline) #Restrict the map after
    # plot(gridmask)

    #consider the grid points in a data frame
    id<- as.data.table(1:ncell(gridmask))
    setnames(id, "V1", "ID")
    griddf<- as.data.table(SpatialPoints(grid))
    griddf<- data.table(id, griddf)
    setnames(griddf, gsub('y','LAT_GRID',gsub('x','LONG_GRID',colnames(griddf))))

    bw      <- 3000
    bw2	<- bw*bw
    #require(mvtnorm)
    #dmvnorm( c(3.84,0) )	# ~ 9.996634e-05 
    threshold <- bw*3.84 	# cut if density is < 1e-4
    threshold <- threshold*threshold	# square the threshold, to avoid sqrt calculations in loop 		
    norm.const	<- 1/(2*pi*bw2)

    tmp		<- griddf[1:1e4,]
    anst <- system.time({
        ans	<- tmp[, {
            z1	<- LONG_GRID - DTnew@coords[,'LONGITUDE_JITTER']
            z2	<- LAT_GRID - DTnew@coords[,'LATITUDE_JITTER']
            z1	<- z1*z1 + z2*z2 		# square distance
            z2	<- which(z1<threshold)	# avoid sqrt on 2e4 entries
            w	<- norm.const*exp(-0.5*z1/bw2)	#now with correct normalising constant
            #	Xiayue
            #z3 <-  z1*z1 + z2*z2
            #z4 <- which(z3<threshold)
            #z <- cbind(matrix(z1[z4],ncol=1),matrix(z2[z4],ncol=1))
            #OR: the source code in Boom seems quite slow, with Cholesky decomposition etc. DIY faster?
            #w <- dmvn(z,mu=c(0,0),bw^2*diag(2))	
            #z2 <- z4
            #	olli
            # z1	<- z1*z1 + z2*z2 		# square distance
            # z2	<- which(z1<threshold)	# avoid sqrt on 2e4 entries
            # # to avoid very large output data, calculate directly all smooths here
            # z1	<- sqrt(z1[z2])			# sqrt on few entries					
            # w	<- dnorm(z1, mean=0, sd=bw) # OR: I agree the normalising constant is not right
            # code assumes @coords and @data has same order. 
            list( 	HIV_STATUS_MEAN=mean( DTnew@data$HIV_STATUS[z2] ),				#no weighting by distance
                 HIV_STATUS_KERNEL=sum( DTnew@data$HIV_STATUS[z2]*w )/sum(w),		#Gaussian kernel
                 VL_COPIES_KERNEL_GEOMMEAN = exp(sum(w*log(DTnew@data$VL_COPIES[z2]+1))/sum(w))-1, #Geometric Mean Kernel
                 VL_DETECTABLE_KERNEL = sum( DTnew@data$VL_DETECTABLE[z2]*w)/sum(w) #Detectable Prevelance
            )
        }, by=c('ID','LONG_GRID','LAT_GRID')]
    })

    grid[]<- ans[, VL_DETECTABLE_KERNEL]
    gridmask<- mask(grid, outline)

    #Breaks chosen by looking at data - need refining
    plot(gridmask,
         breaks = c(0, 0.025, 0.05, 0.075, 0.1, 0.5), 
         col=brewer.pal(11, "RdYlGn")[c(10,9,5,4,3)],
         axes=FALSE, box=FALSE, ylim= c(exmap[3],-6000),legend=FALSE)

    plot(outline, add=TRUE)
    par(xpd=TRUE)
    legend("right", 
           legend=c("0-2.5","2.5-5","5-7.5","7.5-10", ">10"),
           fill=brewer.pal(11, "RdYlGn")[c(10,9,5,4,3)],
           horiz = FALSE, inset=-0.175,
           title= "Prevelence of \n Detectable \n Viremia (%)",
           cex=0.8, box.lty = 0)

    grid[]<- ans[, VL_COPIES_KERNEL_GEOMMEAN]
    gridmask<- mask(grid, outline)
    plot(gridmask, breaks = c(0, 0.8, 1.5, 2.5, 3, 145), col=brewer.pal(11, "RdYlGn")[c(10,9,5,4,3)] ,  axes=FALSE, box=FALSE, ylim= c(exmap[3],-6000), legend=FALSE)
    plot(outline, add=TRUE)
    par(xpd=TRUE)

    legend("right",
           legend=c("0-0.8","0.8-1.5","1.5-2.5","2.5-3", ">3"),
           fill=brewer.pal(11, "RdYlGn")[c(10,9,5,4,3)],
           horiz = FALSE, inset=-0.175,
           title= "Geometric Mean \n VL (Copies/ml)",
           cex=0.8, box.lty = 0)
}

make.map.220810	<- function(DT, plotDT=NA, plotcols=NA)
{
    # DT <- copy(dall); plotDT <- copy(vlc); plotcols=NA
    # plotcols <- 'PVLNSofHIV_MEAN'
    library(rnaturalearth)
    library(osmdata)
    require(ggpubr)
    library(sf)

    if( ! is.data.table(plotDT))
    {
        bool.plot.loc <- TRUE
    }
    if( is.data.table(plotDT) & !is.na(plotcols))
    {
        bool.plot.loc <- FALSE
        stopifnot( all(plotcols %in% names(plotDT)))
        cols <- c('COMM_NUM', 'ROUND', 'SEX', plotcols)
        plotDT <- plotDT[, ..cols]
    }

    # theme
    mytheme <- theme(
                     panel.background = element_rect(fill = 'lightblue', color='blue'),
                     panel.grid.major = element_line(colour = "transparent")
    )
    noticks <- theme(axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank()) 

    # Get longitude and latitude
    infile <- file.path(indir.deepsequence.data, 'RCCS_R15_R18', 'Rakai_community_geography_R15.rda')

    # Fix factors stored as integers (COMM_NUM)
    # TODO: merge plotDT with rest
    tmp <- new.env()
    load(infile, envir=tmp)
    ds <- as.data.table(tmp$comgps)
    cols <- c('latitude', 'longitude')
    ds[, (cols):=lapply(.SD, unlist), .SDcols=cols]
    .f <- function(x) as.integer(as.character(x))
    ds[, COMM_NUM := .f(COMM_NUM)]
    tmp1 <- unique(DT[, .(COMM_NUM, FC)])
    ds <- merge(ds, tmp1)
    ds[is.na(longitude)]
    setnames(ds, cols, paste0(toupper(cols), '_JITTER'))
    ds_sf <- st_as_sf(ds, coords=c("LONGITUDE_JITTER", "LATITUDE_JITTER"), crs=4326)
    ds_sf_t <- st_transform(ds_sf, crs=2163)

    ds_sf_t <- merge(ds_sf_t, plotDT, by='COMM_NUM')

    plotDT[]

    coord2 <- data.table(
                         x=range(ds$LONGITUDE_JITTER)+c(-.1, .1), 
                         y=range(ds$LATITUDE_JITTER) +c(-.1, .1)
    )


    # Latitude from -1.28538 to 3.66088 
    # and longitude from 29.65 to 34.66659.

    long_bounds=c(22, 43) # x
    lati_bounds=c(-10, 10)  # y

    # get uganda borders + lakes
    # __________________________

    uganda_borders <- c('uganda', 'kenya', "Tanzania", 'rwanda', 'democratic republic of the congo')

    africa_sf <- ne_countries(scale=50, continent='africa', returnclass='sf')
    uganda_sf <- ne_countries(scale=50, country=c('uganda', 'tanzania'), returnclass='sf')

    # Define bounding box:
    bb <- c(long_bounds[1], lati_bounds[1],
            long_bounds[2], lati_bounds[2])
    lakes <- opq(bbox=bb, timeout=150) |>
    add_osm_feature(key='natural', value='water') |>
    osmdata_sf()

    # Download data of interest
    lakes10 <- ne_download(scale = 10, type = 'lakes', category = 'physical', returnclass='sf')
    roads10 <- ne_download(scale = 10, type = 'roads', category = 'cultural', returnclass='sf')
    rivers10 <- ne_download(scale=10, type='rivers_lake_centerlines', category='physical', returnclass='sf')

    p1 <- ggplot() +
        geom_sf(data=africa_sf, fill='antiquewhite') + 
        geom_sf(data=lakes10, color="blue", fill="lightblue") +
        geom_sf_text(data=africa_sf, aes(label=name), color = "darkred", fontface = "bold", check_overlap = FALSE) +
        geom_rect(data=coord2, aes(xmin = x[1], xmax = x[2], ymin = y[1], ymax = y[2]), color = "red", fill = NA)  +
        coord_sf(xlim=long_bounds, ylim=lati_bounds) +
        labs(x='', y='') + 
        mytheme

    if('FC.x' %in% names(ds_sf_t))
    {
        setnames(ds_sf_t, 'FC.x', 'FC')
        ds_sf_t$FC.y <- NULL
    }
    p2 <- ggplot() +
        geom_sf(data=africa_sf, fill='antiquewhite') + 
        geom_sf(data=rivers10, color="blue", fill="lightblue") +
        geom_sf(data=roads10, color="grey80") +
        geom_sf(data=lakes10, color="blue", fill="lightblue") +
        geom_sf(data=ds_sf_t, aes(fill=FC), colour='black', pch=21, alpha=.8, size=3)+
        coord_sf(xlim=coord2$x, ylim=coord2$y) +
        scale_fill_manual(values=palettes$comm) + 
        mytheme +
        theme(legend.position=c(0,0),
              legend.background = element_rect(fill = "white", color = "black"),
              legend.justification = c("left", "bottom"),
              panel.border = element_rect(color = "red", size = 2)) +
        labs(fill='Community type')


    # CAN I PUT THE LEGEND ON NORTH-EAST CORNER to save space?
    p <- ggarrange(p1, p2, ncol=1, heights=c(1.1,1))

    if( ! bool.plot.loc )
    {
        filename <- 'uganda_communities_map_220808.pdf'
        ggsave2(p, file=filename, w=5.5, h=7.5)
    }else{

        tmp <- list(
                    list(sex='women', round=16),
                    list(sex='women', round=19),
                    list(sex='men', round=16),
                    list(sex='men', round=19)
        )

        .f <- function(lst)
        {
            idx <- which(ds_sf_t$ROUND == lst$round & ds_sf_t$SEX == lst$sex )
            DS <- ds_sf_t[idx,] 
            ggplot() +
                geom_sf(data=africa_sf, fill='antiquewhite') + 
                geom_sf(data=rivers10, color="blue", fill="lightblue") +
                geom_sf(data=roads10, color="grey80") +
                geom_sf(data=lakes10, color="blue", fill="lightblue") +
                geom_sf(data=DS, aes(fill=PVLNSofHIV_MEAN), colour='black', pch=21, alpha=.8, size=3) + 
                scale_fill_gradient2(low='green', high='red', midpoint=median(DS$PVLNSofHIV_MEAN)) +
                facet_grid( SEX ~ ROUND ) + 
                coord_sf(xlim=coord2$x, ylim=coord2$y) +
                mytheme +
                theme(legend.position=c(0,0),
                      legend.background = element_rect(fill = "white", color = "black"),
                      legend.justification = c("left", "bottom"),
                      panel.border = element_rect(color = "red", size = 2)
                      ) +
                labs(fill='')
        }

        tmp1 <- lapply(tmp, .f)
        p <- ggarrange(tmp1[[1]], tmp1[[2]], 
                       tmp1[[3]], tmp1[[4]],
                       align='hv',
                       nrow=2, ncol=2)
        p <- annotate_figure(p,
                             top=text_grob('Prevalence of unsuppressed among HIV+', 
                                           face='bold', size='14'))

        return(p)
    }
}

glm_random_effects <- function(formula_glm)
{
    glm_reffs <- glmer(data=dglm,
                       formula=formula_glm,
                       weights=N_HIV,
                       control=glmerControl(optimizer='Nelder_Mead'),
                       family='binomial')
    ilink <- family(glm_reffs)$linkinv

    dpreds <- copy(dglm)
    dpreds[, ROUND:=as.factor(ROUND)]
    dpreds[, PRED:=ilink(predict(glm_reffs))]
    dpreds[, `:=` (CL=qbinom(p=.25, size=N_HIV, prob=PRED)/N_HIV,
                   CU=qbinom(p=.975, size=N_HIV, prob=PRED)/N_HIV), by=PRED]

    cat("The proportion of observations falling within prediction intervals are:\n",
        dpreds[, mean(PVLNSofHIV_MEAN <= CU & PVLNSofHIV_MEAN >= CL)], "\n")

    g <- ggplot(dpreds, aes(x=ordered(COMM_NUM,levels=comm_lvls), color=ROUND, pch=FC2) ) +
        geom_linerange(aes(ymin=CL, ymax=CU), position=position_dodge(width=.4)) +
        geom_point(aes(y=PVLNSofHIV_MEAN, size=N_HIV), position=position_dodge(width=.4)) +
        geom_vline(xintercept=6.5, linetype='dotted') + 
        geom_vline(xintercept=36.5, linetype='dotted') + 
        facet_grid(SEX~.) + 
        theme_bw()
    print(g)

    glm_reffs
}

glm_no_random_effects <- function(formula_glm)
{
    glm_base <- glm(data=dglm,
                    formula=formula_glm,
                    weights=N_HIV,
                    family='binomial')
    ilink <- family(glm_base)$linkinv

    dpreds <- copy(dglm)
    dpreds[, ROUND:=as.factor(ROUND)]
    dpreds[, PRED:=ilink(predict(glm_base))]
    dpreds[, `:=` (CL=qbinom(p=.25, size=N_HIV, prob=PRED)/N_HIV,
                   CU=qbinom(p=.975, size=N_HIV, prob=PRED)/N_HIV), by=PRED]

    cat("The proportion of observations falling within prediction intervals are:\n",
        dpreds[, mean(PVLNSofHIV_MEAN <= CU & PVLNSofHIV_MEAN >= CL)], "\n")

    g <- ggplot(dpreds, aes(x=ordered(COMM_NUM, levels=comm_lvls), pch=FC2, color=ROUND) ) +
        geom_linerange(aes(ymin=CL, ymax=CU), position=position_dodge(width=.4)) +
        geom_point(aes(y=PVLNSofHIV_MEAN, size=N_HIV), position=position_dodge(width=.4)) +
        geom_vline(xintercept=6.5, linetype='dotted') + 
        geom_vline(xintercept=36.5, linetype='dotted') + 
        facet_grid(SEX~.) + 
        theme_bw()
    print(g)

    glm_base
}

glm_random_effects_stan <- function(formula_glm)
{
    glm_reffs <- stan_glmer(data=dglm,
                            formula=formula_glm,
                            weights=N_HIV,
                            control=glmerControl(optimizer='Nelder_Mead'),
                            family='binomial')
    ilink <- family(glm_reffs)$linkinv

    dpreds <- copy(dglm)
    dpreds[, ROUND:=as.factor(ROUND)]
    dpreds[, PRED:=ilink(predict(glm_reffs))]
    dpreds[, `:=` (CL=qbinom(p=.25, size=N_HIV, prob=PRED)/N_HIV,
                   CU=qbinom(p=.975, size=N_HIV, prob=PRED)/N_HIV), by=PRED]

    cat("The proportion of observations falling within prediction intervals are:\n",
        dpreds[, mean(PVLNSofHIV_MEAN <= CU & PVLNSofHIV_MEAN >= CL)], "\n")

    g <- ggplot(dpreds, aes(x=ordered(COMM_NUM,levels=comm_lvls), color=ROUND, pch=FC2) ) +
        geom_linerange(aes(ymin=CL, ymax=CU), position=position_dodge(width=.4)) +
        geom_point(aes(y=PVLNSofHIV_MEAN, size=N_HIV), position=position_dodge(width=.4)) +
        geom_vline(xintercept=6.5, linetype='dotted') + 
        geom_vline(xintercept=36.5, linetype='dotted') + 
        facet_grid(SEX~.) + 
        theme_bw()
    print(g)

    glm_reffs
}


mcmc_intervals_2 <- function(
                             x,
                             pars = character(),
                             regex_pars = character(),
                             transformations = list(),
                             ...,
                             prob = 0.5,
                             prob_outer = 0.9,
                             point_est = c("median", "mean", "none"),
                             outer_size = 0.5,
                             inner_size = 2,
                             point_size = 4,
                             rhat = numeric(),
                             color = NA
                             )
{
    .gs <- function(x) as.integer(gsub('[A-z]|\\)|\\(|:','',x))
    p_mint$data$COMM_NUM <- .gs(p$data$parameter)
    tmp <- unique(dglm[, .(COMM_NUM, FC2),])
    p_mint$data <- merge(p$data, tmp)

    out <- ggplot(p_mint$data,
                  aes(x=m, y=parameter, xmin=ll, xmax=hh, col=FC2)) +
        geom_point(size=1.5) +
        geom_errorbarh(height=0, size=0.5) +
        geom_errorbarh(aes(xmin=l, xmax=h), height=0, size=1)

    out
}
