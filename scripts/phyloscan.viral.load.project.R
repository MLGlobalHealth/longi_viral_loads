ggsave2 <- function(p, file, w, h)
{
        filename <- file
        filename2 <- gsub('pdf$','png',filename)
        cat('Saving', filename, '...\n')
        ggsave(p, file=file.path(outdir,filename2), width=w, height=h)
        ggsave(p, file=file.path(outdir, filename), width=w, height=h)	
}


date2numeric<- function( x )
{
	if(!class(x)%in%c('Date','character'))	return( x )
	x	<- as.POSIXlt(x)
	tmp	<- x$year + 1900
	x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
	x	
}

.preprocess.ds.oli <- function(DT)
{
        # DT <- copy(dall)
        DT <- subset(DT, AGEYRS <= 50)

        # remove HIV+ individuals with missing VLs and 
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
                z <- which(DT$ROUND==ROUND, DT$FC==FC & DT$SEX==SEX & DT$AGEYRS==AGEYRS)	
                list(N          = length(z),
                     HIV_N      = sum(DT$HIV_STATUS[z]==1),
                     VLNS_N     = sum(DT$VLNS[z]==1),
                     ARV_N      = sum(DT$ARVMED[z]==0 & DT$HIV_STATUS[z]==1 & !is.na(DT$ARVMED[z])),
                     VL_MEAN    = mean(DT$VLC[z]),
                     VL_SD      = sd(DT$VLC[z]),
                     VL_MEAN_SD = sd(DT$VLC[z]) / sqrt(length(z)),
                     HIV_N      = sum(DT$HIV_STATUS[z]==1),
                )
        }, by=names(vla)]

	setnames(vla, c('FC','SEX','AGEYRS'), c('LOC_LABEL','SEX_LABEL','AGE_LABEL'))
	vla[, LOC:= as.integer(LOC_LABEL=='fishing')]
	vla[, SEX:= as.integer(SEX_LABEL=='M')]
	vla[, AGE:= AGE_LABEL-14L]
	vla[, ROW_ID:= seq_len(nrow(vla))]

        # Extract selected fields
        vla[select]
}


vl.get.eligible.round17<- function()
{
	require(data.table)
	
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/rakai_elibility.rda"
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/"
	load(infile)
	
	#	subset to data of interest
	de	<- as.data.table(eldat)	
	de	<- subset(de, status%in%c('_Participated','Away','Blood refusal','Missing data','Other','Refused','urine sample'))
	de	<- subset(de, visit==17)
	
}

vl.vlprops.by.comm.gender.loc<- function(DT, write.csv=FALSE)
{

	# remove HIV+ individuals with missing VLs
	DT <- subset(DT, HIV_STATUS==0 | HIV_AND_VL==1)
	
	# define VL_COPIES for uninfected
	set(DT, NULL, 'VLC', DT$VL_COPIES)
	set(DT, DT[,which(HIV_STATUS==0)], 'VLC', 0)
	
	# define undetectable VL (machine-undetectable)
	# define suppressed VL (according to WHO criteria)	
	set(DT, NULL, 'VLU', DT[, as.integer(VLC < VL_DETECTABLE)])
	set(DT, NULL, 'VLS', DT[, as.integer(VLC < VIREMIC_VIRAL_LOAD )])
	set(DT, NULL, 'VLD', DT[, as.integer(VLC >= VL_DETECTABLE)])
	set(DT, NULL, 'VLNS', DT[, as.integer(VLC>= VIREMIC_VIRAL_LOAD )])
	set(DT, NULL, 'HIV_AND_VLD', DT[, as.integer(VLD==1 & HIV_AND_VL==1)])
        
        # cols <- c('VLU', 'VLS', 'VLD', 'VLNS', 'HIV_AND_VLD')
        # DT[, lapply(.SD, mean, na.rm=T) , .SDcols=cols]
	
	# reset VLC below machine detectable to 0
	set(DT, DT[, which(HIV_AND_VL==1 & VLU==1)], 'VLC', 0)
	setkey(DT, FC, SEX, AGEYRS)
	
	# merge two communities that fully overlap, so we have 40 communities in the end 
	set(DT, DT[, which(COMM_NUM==22)], 'COMM_NUM', 1)

	# calculate HIV prevalence and proportion not suppressed of HIV+ by community and gender
        .f <- function(x,y)
        {
                as.vector( unname ( binconf(sum(x), length(y) )  ) )
        }
        

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
        vlc[, factor(SEX, levels = c('M', 'F'), labels=c('men', 'women'))]
	set(vlc, NULL, 'SEX', vlc[, factor(SEX, levels=c('M','F'), labels=c('men','women'))])

	
	p <- ggplot(vlc) +
                scale_x_continuous(labels=scales:::percent) +
                scale_y_continuous(labels=scales:::percent) +
                geom_errorbar(aes(x=PHIV_MEAN, ymin=PVLNSofHIV_CL, ymax=PVLNSofHIV_CU), alpha=0.2) +
                geom_errorbarh(aes(y=PVLNSofHIV_MEAN, xmin=PHIV_CL, xmax=PHIV_CU), alpha=0.2) +
                geom_point(aes(x=PHIV_MEAN, y=PVLNSofHIV_MEAN, colour=FC)) +
                geom_text(aes(x=PHIV_MEAN, y=PVLNSofHIV_MEAN, label=COMM_NUM), size=2) +
                facet_grid(ROUND~SEX) +
                scale_colour_manual(values=palettes$comm) + 
                theme_bw() +
                theme(legend.position='bottom') + 
                labs(x='\nHIV prevalence', 
                     y='proportion unsuppressed HIV among infected\n', 
                     colour='location')

        filename <- file.path(out.dir,'220829_hivnotsuppofhiv_vs_hivprev_by_round_gender_fishinland.pdf')
        filename2 <- gsub('pdf$','png',filename)
        ggsave(p, filename=filename2, w=9, h=12)
        ggsave(p, filename=filename, w=9, h=12)
		
        if(write.csv)
        {
                #	write results to file
                filename <- file.path(out.dir,'220829_hivnotsuppofhiv_vs_hivprev_by_round_gender_fishinland.csv')
                write.csv(vlc, file=filename)
        }

        list(DT=vlc, p=p)
}

vl.prevalence.by.gender.loc.age.gp <- function(DT) 
{
        # DT <- copy(dall)
	outdir <- file.path(out.dir)
        DT <- .preprocess.ds.oli(DT)
	
        tmp <- c('N', 'HIV_N', 'VLNS_N', 'ARV_N')
        vla <- .preprocess.make.vla(DT, select=tmp)

        # Stan file locations
        file.stan.1 <- file.path(path.stan, 'vl_prevalence_by_gender_loc_age_gp_1.stan')
        file.stan.2 <- file.path(path.stan, 'vl_prevalence_by_gender_loc_age_gp_2.stan')

	#stan.model <- stan_model(file=file.stan.1, model_name= 'gp_one')
	stan.model <- stan_model(file=file.stan.2, model_name= 'gp_all')
	vla2 <- subset(vla, SEX==0 & LOC==0)

        .fit.stan.and.plot.by.round <- function(DT, iter=2e3, warmup=5e2, chains=1, control=list(max_treedepth= 15, adapt_delta= 0.999))
        {
                #  DT <- copy(vla[ROUND == 16])
                round <- DT[, unique(ROUND)]
                stopifnot(length(round) == 1)
                cat('Fitting stan model for round ', round, '\n')

                # make stan data
                # ______________

                tmp <- seq(DT[, min(AGE_LABEL)], DT[, max(AGE_LABEL)+1], 0.5)
                tmp1 <- which(tmp%%1==0.5)

                stan.data <- list(
                        x_predict = tmp,
                        N_predict = length(tmp),
                        observed_idx = tmp1,
                        N_observed = length(tmp1),
                        y_observed_00 = DT[SEX==0 & LOC==0, HIV_N],
                        y_observed_10 = DT[SEX==1 & LOC==0, HIV_N],
                        y_observed_01 = DT[SEX==0 & LOC==1, HIV_N],
                        y_observed_11 = DT[SEX==1 & LOC==1, HIV_N],
                        total_observed_00 = DT[SEX==0 & LOC==0, N],
                        total_observed_10 = DT[SEX==1 & LOC==0, N],
                        total_observed_01 = DT[SEX==0 & LOC==1, N],
                        total_observed_11 = DT[SEX==1 & LOC==1, N],
                        alpha_hyper_par_00 = 2,
                        alpha_hyper_par_10 = 2,
                        alpha_hyper_par_01 = 2,
                        alpha_hyper_par_11 = 2,
                        rho_hyper_par_00 = diff(range(tmp))/3,
                        rho_hyper_par_10 = diff(range(tmp))/3,
                        rho_hyper_par_01 = diff(range(tmp))/3,
                        rho_hyper_par_11 = diff(range(tmp))/3
                )

                fit <- sampling(stan.model,
                                data=stan.data,
                                iter=iter, warmup=warmup, chains=chains, 
                                control = control)

                filename <- file.path(outdir, paste0("hivprevalence_gp_stanfit_round",round,"_220829.rds"))
                saveRDS(fit, file=filename)

                # Analyse posterior
                # _________________

                cat('The minimum effective sample size is:\n')
                min( summary(fit)$summary[, 'n_eff'] )	
                re <- rstan::extract(fit)
                ps <- c(0.025,0.25,0.5,0.75,0.975)

                # extract hyperparams rho		
                .f <- function(x) quantile(x, probs=ps)
                tmp <- cbind( .f(re$rho_00),   .f(re$rho_10),   .f(re$rho_01),   .f(re$rho_11),
                              .f(re$alpha_00), .f(re$alpha_10), .f(re$alpha_01), .f(re$alpha_11) )			
                colnames(tmp) <- c('rho_00','rho_10','rho_01','rho_11','alpha_00','alpha_10','alpha_01','alpha_11')
                rownames(tmp) <- c('CL','IL','M','IU','CU')
                tmp <- as.data.table(reshape2::melt(tmp))
                setnames(tmp, 'Var2', 'GP_hyper_par')
                tmp[, SEX := as.integer(gsub('^([a-z]+)_([0-9])([0-9])','\\2',GP_hyper_par))]
                tmp[, LOC := as.integer(gsub('^([a-z]+)_([0-9])([0-9])','\\3',GP_hyper_par))]
                tmp[, GP_hyper_par := gsub('^([a-z]+)_([0-9])([0-9])','\\1',GP_hyper_par)]
                tmp <- dcast.data.table(tmp, LOC+SEX+GP_hyper_par~Var1, value.var='value')
                prev.hiv.gp.pars <- merge(unique(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), tmp, by=c('SEX','LOC'))
                p <- ggplot(prev.hiv.gp.pars, aes(x= paste0(GP_hyper_par, ' ', LOC_LABEL, ' ', SEX_LABEL))) +
                                geom_point(aes(y=M)) +
                                geom_errorbar(aes(ymin=CL, ymax=CU)) +
                                coord_flip() +
                                theme_bw() +
                                labs(x='GP hyperparameter\n', y='')

                filename=paste0('220829_hivprevalence_gppars_round',round,'.pdf')
                ggsave2(p, file=filename, w=6, h=3)

                # make prevalence plot by age
                .f <- function(x) apply(x, 2, quantile, probs=ps)
                tmp <- cbind( .f(re$p_predict_00), .f(re$p_predict_10),
                              .f(re$p_predict_01), .f(re$p_predict_11))
                rownames(tmp) <- c('CL','IL','M','IU','CU')
                tmp <- as.data.table(reshape2::melt(tmp))
                prev.hiv.by.age <- dcast.data.table(tmp, Var2~Var1, value.var='value')
                tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
                prev.hiv.by.age <- cbind(tmp, prev.hiv.by.age) 
                prev.hiv.by.age <- merge(unique(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), prev.hiv.by.age, by=c('SEX','LOC'))

                p <- ggplot(prev.hiv.by.age) + 		
                        geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2) +
                        geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL)) +
                        scale_x_continuous( expand=c(0,0) ) + 
                        scale_y_continuous(labels=scales:::percent) +
                        scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
                        facet_wrap(~LOC_LABEL, ncol=2) +
                        theme_bw() +
                        labs(x='\nage at visit (years)', 
                             y='HIV prevalence (95% credibility interval)\n', 
                             colour='gender', 
                             linetype='location')
                filename=paste0('220829_hivprevalence_vs_age_by_gender_fishinland_stan_round',round,'.pdf')
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

                # plot prevalence ratio female:male and male:female by age
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
                        scale_y_log10() +
                        facet_wrap(~LOC_LABEL, ncol=2) +
                        theme_bw() +
                        labs(x='\nage at visit (years)', 
                             y='female to male HIV prevalence ratio\n(95% credibility interval)\n')

                filename=paste0('220829_hivprevalenceratio_vs_age_by_fishinland_stan_round',round,'.pdf')
                ggsave2(p, file=filename, w=8, h=5)

                #	extract if difference in female:male prevalence risk ratio in fishing vs inland
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
                

                filename=file.path(outdir, "220929f_hivprevalence.rda")
                save(DT, re, prev.hiv.by.age, prevratio.hiv.by.loc,
                     prev.hiv.by.sex.loc, prevratio.hiv.by.loc.age,
                     file=filename)


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

                filename = file.path(outdir, paste0("220829f_hivprevalence_round",round,".csv"))
                fwrite(dt, row.names=FALSE, file=filename)
	
                TRUE
        }
        
        foreach(
                r = vla[, unique(ROUND)],
                .combine='c'
        ) %dopar% {
                cat('Running Round', r, '\n')
                .fit.stan.and.plot.by.round(vla[ ROUND ==r, ], iter=2e3)
        } -> tmp

        return(tmp)

}

vl.prevalence.by.gender.loc.age.icar<- function(DT)
{
        # DT <- copy(dall)
	outdir <- file.path(out.dir)
        DT <- .preprocess.ds.oli(DT)
	
	# subset(DT, FC=='inland' & SEX=='F' & AGEYRS==25)
        tmp <- seq.int(min(DT$AGEYRS), max(DT$AGEYRS))
        tmp1 <- DT[, sort(unique(ROUND))]
	vla <- as.data.table(expand.grid(ROUND=tmp1,
                                         FC=c('fishing','inland'),
                                         SEX=c('M','F'),
                                         AGEYRS=tmp))
	vla <- vla[, { 
                z<- which(DT$ROUND==ROUND, DT$FC==FC & DT$SEX==SEX & DT$AGEYRS==AGEYRS)	
                list(N      = length(z),
                     HIV_N  = sum(DT$HIV_STATUS[z]==1),
                     VLNS_N = sum(DT$VLNS[z]==1),
                     ARV_N  = sum(DT$ARVMED[z]==0 & DT$HIV_STATUS[z]==1 & !is.na(DT$ARVMED[z]))
                     )},
            by=names(vla)]

	setnames(vla, c('FC','SEX','AGEYRS'), c('LOC_LABEL','SEX_LABEL','AGE_LABEL'))
	vla[, LOC:= as.integer(LOC_LABEL=='fishing')]
	vla[, SEX:= as.integer(SEX_LABEL=='M')]
	vla[, AGE:= AGE_LABEL-14L]
	vla[, ROW_ID:= seq_len(nrow(vla))]

        # Stan file locations
        file.stan.1 <- file.path(path.stan, 'vl_prevalence_by_gender_loc_age_icar_1.stan')
        file.stan.2 <- file.path(path.stan, 'vl_prevalence_by_gender_loc_age_icar_2.stan')

        .fit.stan.and.plot.by.round <- function(DT, iter=2e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999))
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

                #save(fit, file=file.path(outdir, "hivprevalence_icar_stanfit_200428.rda"))		# trends by age quite rough, using Cauchy prior on sigma
                #save(fit, file=file.path(outdir, "hivprevalence_icar_stanfit_200428c.rda"))	# trends by age still quite rough, using N(0,0.1) prior on sigma
                filename <- paste0("hivprevalence_icar_stanfit_round",round,"_220829.rds")
                saveRDS(fit, file=file.path(outdir, filename))

                min( summary(fit)$summary[, 'n_eff'] )
                re <- rstan::extract(fit)
                ps <- c(0.025,0.5,0.975)
                
                quantile(re$sigma_loc0, probs=ps)
                #	in 200428b: 0.2230766 0.2943186 0.3861041 
                #	in 200428c: 0.2022420 0.2578030 0.3270586  
                #	in 200428d: 0.4380249 0.5558544 0.7091595
                
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

                filename=paste0('220829_hivprevalence_vs_age_by_gender_fishinland_icar_round',round,'.pdf')
                filename=file.path(outdir, filename)
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
                #	   LOC SEX LOC_LABEL SEX_LABEL         CL         CU         M                 LABEL
                #1:   0   0    inland         F 0.15553741 0.17130775 0.1634250 0.16% (0.16% - 0.17%)
                #2:   0   1    inland         M 0.08610792 0.09970375 0.0927843  0.09% (0.09% - 0.1%)
                #3:   1   0   fishing         F 0.42138562 0.46424074 0.4427624 0.44% (0.42% - 0.46%)
                #4:   1   1   fishing         M 0.30626567 0.34474245 0.3254705 0.33% (0.31% - 0.34%)


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
                #	   LOC LOC_LABEL variable        CL        CU         M              LABEL
                #1:   0    inland    PR_FM 1.6122846 1.9240272 1.7615063 1.76 (1.61 - 1.92)
                #2:   0    inland    PR_MF 0.5197432 0.6202379 0.5676960 0.57 (0.52 - 0.62)
                #3:   1   fishing    PR_FM 1.2608839 1.4687378 1.3605713 1.36 (1.26 - 1.47)
                #4:   1   fishing    PR_MF 0.6808567 0.7930944 0.7349854 0.73 (0.68 - 0.79)
                

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
                filename=paste0('220829_hivprevalenceration_vs_age_by_fishinland_icar_round',round,'.pdf') 
                filename=file.path(outdir, filename)
                cat('Saving', filename, '...\n')
                filename2 <- gsub('pdf$','png',filename)
                ggsave(p, filename=filename2, w=10, h=5)
                ggsave(p, filename=filename, w=10, h=5)

                # extract if difference in female:male prevalence risk ratio in fishing vs inland
                #________________________________________________________________________________
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
                #1: 0.2182593 CL
                #2: 0.4008081  M
                #3: 0.5894657 CU
                
                filename=file.path(outdir, paste0('hivprevalence_round',round,'220829.rda'))
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
                .fit.stan.and.plot.by.round(vla[ ROUND == r, ], 
                                            iter=2e3, warmup=5e2, chains=1,
                                            control = list(max_treedepth= 15, adapt_delta= 0.999)
                                            )

        } -> tmp
        tmp

        return(tmp)
}

vl.meanviralload.by.gender.loc.age.icar<- function(DT)
{
	
        # DT <- copy(dall)
	outdir <- file.path(out.dir)
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
        filename = file.path(outdir, filename)
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
	        saveRDS(fit, file=file.path(outdir, filename) )
	
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
                filename=file.path(outdir, filename)
                filename2 <- gsub('pdf$','png',filename)
                ggsave(p, filename=filename2, w=6, h=5)
	        ggsave(p, filename=filename, w=6, h=5)

                TRUE
        }

        foreach(
                r = vla[, unique(ROUND)],
                .combine='c'
        ) %dopar% {
                cat('Running Round', r, '\n')
                .fit.stan.and.plot.by.round(vla[ ROUND ==r, ], 
                                            iter=2e3, warmup=5e2, chains=1)
        } -> tmp

        return(tmp)
}

vl.suppofinfected.by.gender.loc.age.icar<- function()
{
        # DT <- copy(dall)
	outdir <- file.path(out.dir)
        DT <- .preprocess.ds.oli(DT)
	
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
                     ARV_N      = sum(DT$ARVMED[z]==0 & DT$HIV_STATUS[z]==1 & !is.na(DT$ARVMED[z]))
                )				
        }, by=names(vla)]

	setnames(vla, c('FC','SEX','AGEYRS'), c('LOC_LABEL','SEX_LABEL','AGE_LABEL'))
	vla[, LOC:= as.integer(LOC_LABEL=='fishing')]
	vla[, SEX:= as.integer(SEX_LABEL=='M')]
	vla[, AGE:= AGE_LABEL-14L]
	vla[, ROW_ID:= seq_len(nrow(vla))]
	
        # Stan file locations
        file.stan.1 <- file.path(path.stan, 'vl_suppofinfected_by_gender_loc_age_icar_1.stan')
	
        list.files(path.stan, pattern='suppofinfected')
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
                saveRDS(fit, file=file.path(outdir,filename))

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
                fit2 <- sampling(stan.model, data=stan.data, iter=iter, warmup=warmup, chains=chains, control=control)	

                filename=paste0('notARVAmongInfected_icar_stan_round',round,'_220729.rds')
                save(fit2, file=file.path(outdir,filename))

                
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
                
                filename <- paste0('220729d_notsuppAmongInfected_vs_age_by_gender_fishinland_v1_round',round,'.pdf')
                cat('Saving', filename, '...\n')
                filename2 <- gsub('pdf$','png',filename)
                ggsave(p, file=file.path(out.dir, filename2), w=6, h=5)
                ggsave(p, file=file.path(out.dir, filename), w=6, h=5)	

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
                
                filename <- paste0('220729d_notsuppAmongInfected_vs_age_by_gender_fishinland_v2_round',round,'.pdf')
                cat('Saving', filename, '...\n')
                filename2 <- gsub('pdf$','png',filename)
                ggsave(p, file=file.path(out.dir, filename2), w=6, h=5)
                ggsave(p, file=file.path(out.dir,filename), w=6, h=5)	

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
                filename <- paste0('220729d_suppAmongInfected_vs_age_by_gender_fishinland_stan_v2_round',round,'.pdf')
                cat('Saving', filename, '...\n')
                filename2 <- gsub('pdf$','png',filename)
                ggsave(p, file=file.path(out.dir, filename2), w=6, h=5)
                ggsave(p, file=file.path(out.dir,filename), w=6, h=5)

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
                
                filename <- paste0('220729d_notsuppAmongInfected_vs_age_by_gender_fishinland_v3_round',round,'.pdf')
                cat('Saving', filename, '...\n')
                filename2 <- gsub('pdf$','png',filename)
                ggsave(p, file=file.path(out.dir, filename2), w=9, h=8)
                ggsave(p, file=file.path(out.dir,filename), w=9, h=8)

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
                cat('Saving', filename, '...\n')
                filename2 <- gsub('pdf$','png',filename)
                ggsave(p, file=file.path(out.dir, filename2), w=9, h=8)
                ggsave(p, file=file.path(out.dir,filename), w=9, h=8)
                

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

                filename <- file.path(out.dir, paste0("notsuppamonginfected_220729_round",round,".rda"))
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

                filename <- file.path(outdir, paste0("notsuppamonginfected_round",round,"_220729.csv"))
                fwrite(dt, row.names=FALSE, file=file.path(outdir,filename))
                        
                #	make table version suppressed
                nsprev.by.age[, LABEL:= paste0(round((1-M)*100, d=1),' (',round((1-CU)*100, d=1),' - ',round((1-CL)*100,d=1),')') ]
                set(nsprev.by.age, NULL, 'SEX_LABEL', nsprev.by.age[, factor(as.character(SEX_LABEL), levels=c('F','M'))])	
                sprev.ratio.by.loc.age[, LABEL2:= paste0(round(M, d=2),' (',round(CL, d=2),' - ',round(CU,d=2),')') ]
                dt <- subset(nsprev.by.age, AGE_LABEL%in%c(20,25,30,35,40,45))	
                dt <- dcast.data.table(dt, LOC_LABEL+AGE_LABEL~SEX_LABEL, value.var='LABEL')
                tmp <- subset(sprev.ratio.by.loc.age, variable=="PR_FM" & AGE_LABEL%in%c(20,25,30,35,40,45), c(LOC_LABEL, AGE_LABEL, LABEL2))
                dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
                tmp <- subset(sprev.ratio.by.loc.age, variable=="PR_MF" & AGE_LABEL%in%c(20,25,30,35,40,45), c(LOC_LABEL, AGE_LABEL, LABEL2))
                setnames(tmp, 'LABEL2', 'PR_MF')
                dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))

                filename <- file.path(outdir, paste0("suppamonginfected_round",round,"_220729.csv"))
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

vl.suppofinfected.by.gender.loc.age.gp<- function(DT)
{
	
        # DT <- copy(dall)
	outdir <- file.path(out.dir)
        DT <- .preprocess.ds.oli(DT)
	require(data.table)
	
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
                                     ARV_N      = sum(DT$ARVMED[z]==0 & DT$HIV_STATUS[z]==1 & !is.na(DT$ARVMED[z]))
                                )				
			}, by=names(vla)]

	setnames(vla, c('FC','SEX','AGEYRS'), c('LOC_LABEL','SEX_LABEL','AGE_LABEL'))
	vla[, LOC:= as.integer(LOC_LABEL=='fishing')]
	vla[, SEX:= as.integer(SEX_LABEL=='M')]
	vla[, AGE:= AGE_LABEL-14L]
	vla[, ROW_ID:= seq_len(nrow(vla))]
	

        file.stan.1 <- file.path(path.stan, 'vl_suppofinfected_by_gender_loc_age_gp.stan')
	stan.model <- stan_model(file.stan.1, model_name='gp_all')	
		
	.fit.stan.and.plot.by.round <- function(DT , iter=10e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999))
        {
                # DT <- copy(DT[ROUND == 16] )
                round <- DT[, unique(ROUND)]
                stopifnot(length(round) == 1)
                cat('Fitting stan model for round ', round, '\n')

                stan.data <- list(
                        x_predict = seq(DT[, min(AGE_LABEL)], DT[, max(AGE_LABEL)+1], 0.5),
                        y_observed_00 = DT[SEX==0 & LOC==0, HIV_N-VLNS_N],
                        y_observed_10 = DT[SEX==1 & LOC==0, HIV_N-VLNS_N],
                        y_observed_01 = DT[SEX==0 & LOC==1, HIV_N-VLNS_N],
                        y_observed_11 = DT[SEX==1 & LOC==1, HIV_N-VLNS_N],
                        total_observed_00 = DT[SEX==0 & LOC==0, HIV_N],
                        total_observed_10 = DT[SEX==1 & LOC==0, HIV_N],
                        total_observed_01 = DT[SEX==0 & LOC==1, HIV_N],
                        total_observed_11 = DT[SEX==1 & LOC==1, HIV_N],
                        alpha_hyper_par_00 = 2,
                        alpha_hyper_par_10 = 2,
                        alpha_hyper_par_01 = 2,
                        alpha_hyper_par_11 = 2
                )
                stan.data$N_predict <- length(stan.data$x_predict)
                stan.data$observed_idx <- which(stan.data$x_predict%%1==0.5)
                stan.data$N_observed <- length(stan.data$observed_idx)
                stan.data$rho_hyper_par_00 <- diff(range(stan.data$x_predict))/3
                stan.data$rho_hyper_par_10 <- diff(range(stan.data$x_predict))/3
                stan.data$rho_hyper_par_01 <- diff(range(stan.data$x_predict))/3
                stan.data$rho_hyper_par_11 <- diff(range(stan.data$x_predict))/3
                
                fit <- sampling(stan.model, data=stan.data, iter=iter, warmup=warmup, chains=chains, control=control)
                filename <- paste0( '220729f_notsuppAmongInfected_gp_stan_round',round,'.rds')
                saveRDS(fit, file=file.path(outdir,filename))
                
                
                # compare to self-report
                # ______________________

                stan.data <- list(
                        x_predict = seq(DT[, min(AGE_LABEL)], DT[, max(AGE_LABEL)+1], 0.5),
                        y_observed_00 = DT[SEX==0 & LOC==0, HIV_N-ARV_N],
                        y_observed_10 = DT[SEX==1 & LOC==0, HIV_N-ARV_N],
                        y_observed_01 = DT[SEX==0 & LOC==1, HIV_N-ARV_N],
                        y_observed_11 = DT[SEX==1 & LOC==1, HIV_N-ARV_N],
                        total_observed_00 = DT[SEX==0 & LOC==0, HIV_N],
                        total_observed_10 = DT[SEX==1 & LOC==0, HIV_N],
                        total_observed_01 = DT[SEX==0 & LOC==1, HIV_N],
                        total_observed_11 = DT[SEX==1 & LOC==1, HIV_N],
                        alpha_hyper_par_00 = 2,
                        alpha_hyper_par_10 = 2,
                        alpha_hyper_par_01 = 2,
                        alpha_hyper_par_11 = 2
                )
                stan.data$N_predict <- length(stan.data$x_predict)
                stan.data$observed_idx <- which(stan.data$x_predict%%1==0.5)
                stan.data$N_observed <- length(stan.data$observed_idx)
                stan.data$rho_hyper_par_00 <- diff(range(stan.data$x_predict))/3
                stan.data$rho_hyper_par_10 <- diff(range(stan.data$x_predict))/3
                stan.data$rho_hyper_par_01 <- diff(range(stan.data$x_predict))/3
                stan.data$rho_hyper_par_11 <- diff(range(stan.data$x_predict))/3

                fit2 <- sampling(stan.model, data=stan.data, iter=iter, warmup=warmup, chains=chains, control=control)
                filename <- paste0('200428f_notARVAmongInfected_icar_stan_round',round,'.rds')
                saveRDS(fit2, file=file.path(outdir,filename))
                
                # Analyse posterior
                # _________________

                re <- rstan::extract(fit)
                re2 <- rstan::extract(fit2)
                ps <- c(0.025,0.5,0.975)
                tmp <- summary(fit)$summary
                tmp[grepl('^p_predict_',rownames(tmp)),]
                
                #
                #	extract hyperparams rho
                ps <- c(0.025,0.25,0.5,0.75,0.975)
                tmp <- cbind(quantile(re$rho_00, probs=ps),
                             quantile(re$rho_10, probs=ps),
                             quantile(re$rho_01, probs=ps),
                             quantile(re$rho_11, probs=ps),
                             quantile(re$alpha_00, probs=ps),
                             quantile(re$alpha_10, probs=ps),
                             quantile(re$alpha_01, probs=ps),
                             quantile(re$alpha_11, probs=ps) )			

                colnames(tmp) <- c('rho_00','rho_10','rho_01','rho_11','alpha_00','alpha_10','alpha_01','alpha_11')
                rownames(tmp) <- c('CL','IL','M','IU','CU')
                tmp <- as.data.table(reshape2::melt(tmp))
                setnames(tmp, 'Var2', 'GP_hyper_par')
                tmp[, SEX:= as.integer(gsub('^([a-z]+)_([0-9])([0-9])','\\2',GP_hyper_par))]
                tmp[, LOC:= as.integer(gsub('^([a-z]+)_([0-9])([0-9])','\\3',GP_hyper_par))]
                tmp[, GP_hyper_par:= gsub('^([a-z]+)_([0-9])([0-9])','\\1',GP_hyper_par)]
                tmp <- dcast.data.table(tmp, LOC+SEX+GP_hyper_par~Var1, value.var='value')
                nsinf.gp.pars <- merge(unique(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), tmp, by=c('SEX','LOC'))
                # Make plots
                # __________

                p <- ggplot(nsinf.gp.pars, aes(x= paste0(GP_hyper_par, ' ', LOC_LABEL, ' ', SEX_LABEL))) +
                        geom_point(aes(y=M)) +
                        geom_errorbar(aes(ymin=CL, ymax=CU)) +
                        coord_flip() +
                        theme_bw() +
                        labs(x='GP hyperparameter\n', y='')

                filename <- paste0('220729f_notsuppAmongInfected_gppars_round',round,'.pdf')
                cat('Saving', filename, '...\n')
                filename2 <- gsub('pdf$','png',filename)
                ggsave(p, file=file.path(out.dir, filename2), w=6, h=3)
                ggsave(p, file=file.path(outdir, filename), w=6, h=3)
                
                
                #	make prevalence plot by age
                tmp <- cbind(apply(re$p_predict_00, 2, quantile, probs=ps),
                             apply(re$p_predict_10, 2, quantile, probs=ps),
                             apply(re$p_predict_01, 2, quantile, probs=ps),
                             apply(re$p_predict_11, 2, quantile, probs=ps))

                rownames(tmp) <- c('CL','IL','M','IU','CU')
                tmp <- as.data.table(reshape2::melt(tmp))
                nsinf.by.age <- dcast.data.table(tmp, Var2~Var1, value.var='value')
                tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
                nsinf.by.age <- cbind(tmp, nsinf.by.age) 
                nsinf.by.age <- merge(unique(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), nsinf.by.age, by=c('SEX','LOC'))
                nsinf.by.age[, STAT:='VLNS']
                
                tmp <- cbind( apply(re2$p_predict_00, 2, quantile, probs=ps),
                                apply(re2$p_predict_10, 2, quantile, probs=ps),
                                apply(re2$p_predict_01, 2, quantile, probs=ps),
                                apply(re2$p_predict_11, 2, quantile, probs=ps)
                                )
                rownames(tmp) <- c('CL','IL','M','IU','CU')
                tmp <- as.data.table(reshape2::melt(tmp))
                nainf.by.age <- dcast.data.table(tmp, Var2~Var1, value.var='value')
                tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
                nainf.by.age <- cbind(tmp, nainf.by.age) 
                nainf.by.age <- merge(unique(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), nainf.by.age, by=c('SEX','LOC'))
                nainf.by.age[, STAT:='VLNA']
                tmp <- subset(nainf.by.age, select=c(SEX, LOC, AGE_LABEL, M, CL, CU))
                setnames(tmp, c('M','CL','CU'), c('M2','CL2','CU2'))
                tmp <- merge(nsinf.by.age, tmp, by=c('SEX','LOC','AGE_LABEL'))			
                p <- ggplot(tmp) + 		
                        geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2) +
                        geom_hline(yintercept=c(0.9^3, 0.95^3)) +
                        geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL)) +
                        scale_x_continuous( expand=c(0,0) ) + 
                        scale_y_continuous(label=scales:::percent) +
                        scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
                        facet_wrap(~LOC_LABEL, ncol=2) +
                        theme_bw() +
                        labs(x='\nage at visit (years)', 
                             y='HIV+ individuals with suppressed viral load\n(95% credibility interval)\n', 
                             colour='gender', 
                             linetype='location')

                filename <- paste0('220729f_suppAmongInfected_vs_age_by_gender_fishinland_stan_round',round,'.pdf')
                cat('Saving', filename, '...\n')
                filename2 <- gsub('pdf$','png',filename)
                ggsave(p, file=file.path(out.dir, filename2), w=6, h=5)
                ggsave(p, file=file.path(out.dir,filename), w=6, h=5)		

                p <- ggplot(tmp) + 		
                        geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2) +
                        geom_hline(yintercept=c(0.9^3, 0.95^3)) +
                        geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL)) +
                        geom_line(aes(x=AGE_LABEL, y=M2, colour=SEX_LABEL), linetype=2) +
                        scale_x_continuous( expand=c(0,0) ) + 
                        scale_y_continuous(label=scales:::percent) +
                        scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
                        facet_wrap(~LOC_LABEL, ncol=2) +
                        theme_bw() +
                        labs(x='\nage at visit (years)', 
                             y='HIV+ individuals with suppressed viral load\n(95% credibility interval)\n', 
                             colour='gender')

                filename <- paste0('220729f_suppAmongInfected_vs_age_by_gender_fishinland_stan_v2_round',round,'.pdf')
                cat('Saving', filename, '...\n')
                filename2 <- gsub('pdf$','png',filename)
                ggsave(p, file=file.path(out.dir, filename2), w=9, h=8)
                ggsave(p, file=file.path(outdir,filename), w=6, h=5)	

                tmp <- rbind(nsinf.by.age, nainf.by.age, fill=TRUE)
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

                filename <- paste0('220729f_suppAmongInfected_vs_age_by_gender_fishinland_stan_v3_round',round,'.pdf')
                filename2 <- gsub('pdf$','png',filename)
                cat('Saving', filename, '...\n')
                ggsave(p, file=file.path(outdir,filename2), w=9, h=8)
                ggsave(p, file=file.path(outdir,filename), w=9, h=8)
                
                
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
                rp[, LABEL:= paste0(round(M, d=2),' (',round(CL, d=2),' - ',round(CU,d=2),')') ]	
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
                rp[, PR_MF:=M/F]
                rp <- melt(rp, id.vars=c('LOC_LABEL','AGE_LABEL','iterations'), measure.vars=c('PR_FM','PR_MF'))
                rp <- rp[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','AGE_LABEL','variable')]	
                rp <- dcast.data.table(rp, LOC_LABEL+AGE_LABEL+variable~P, value.var='Q')
                rp[, LABEL:= paste0(round(M, d=2),' (',round(CL, d=2),' - ',round(CU,d=2),')') ]	
                nsinf.ratio.by.loc.age <- copy(rp)
                
                
                # extract if difference in female:male risk ratio of unsuppressed VL is different in fishing vs inland
                # ____________________________________________________________________________________________________

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
                save(DT, re, re2, nainf.by.age,
                     nsinf.by.age, nsinf.by.sex.loc,
                     nsinf.by.loc, nsinf.ratio.by.loc.age,
                     file=file.path(outdir,filename))
                
                
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

                filename <- paste0('220729f_suppamonginfected_round',round,'csv')
                fwrite(dt, row.names=FALSE, file=file.path(outdir,filename))	

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

vl.suppofpop.by.gender.loc.age.gp<- function(DT)
{
	
        # DT <- copy(dall)
	outdir <- file.path(out.dir)
        DT <- .preprocess.ds.oli(DT)

        tmp <- c('N', 'HIV_N', 'VLNS_N', 'ARV_N')
        vla <- .preprocess.make.vla(DT, select=tmp)
	
        # Stan file locations
        file.stan.1 <- file.path(path.stan, 'vl_suppofpop_by_gender_loc_age_gp_1.stan')
	stan.model <- stan_model(file=file.stan.1, model_name= 'gp_all')

        .fit.stan.and.plot.by.round <- function(DT, iter=10e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999))
        {
                #  DT <- copy(vla[ROUND == 16])
                round <- DT[, unique(ROUND)]
                stopifnot(length(round) == 1)
                cat('Fitting stan model for round ', round, '\n')

                # Define Stan Data
                tmp <- seq(DT[, min(AGE_LABEL)], DT[, max(AGE_LABEL)+1], 0.5)
                tmp1 <- length(which(tmp%%1==0.5))
                stan.data <- list(
                        x_predict = tmp,
                        N_predict = length(tmp),
                        observed_idx = which(tmp%%1==0.5),
                        N_observed = tmp1,
                        y_observed_00 = DT[SEX==0 & LOC==0, VLNS_N],
                        y_observed_10 = DT[SEX==1 & LOC==0, VLNS_N],
                        y_observed_01 = DT[SEX==0 & LOC==1, VLNS_N],
                        y_observed_11 = DT[SEX==1 & LOC==1, VLNS_N],
                        total_observed_00 = DT[SEX==0 & LOC==0, N],
                        total_observed_10 = DT[SEX==1 & LOC==0, N],
                        total_observed_01 = DT[SEX==0 & LOC==1, N],
                        total_observed_11 = DT[SEX==1 & LOC==1, N],
                        rho_hyper_par_00 = diff(range(tmp))/3,
                        rho_hyper_par_10 = diff(range(tmp))/3,
                        rho_hyper_par_01 = diff(range(tmp))/3,
                        rho_hyper_par_11 = diff(range(tmp))/3,
                        alpha_hyper_par_00 = 2,
                        alpha_hyper_par_10 = 2,
                        alpha_hyper_par_01 = 2,
                        alpha_hyper_par_11 = 2
                )

                fit <- sampling(stan.model, data=stan.data, iter=iter, warmup= warmup, chains=chains, control=control)

                filename <- paste0( "220729f_suppAmongPop_gp_stan_round",round,".rds") 
                saveRDS(fit, file=file.path(outdir, filename))

                #####################
                # Analyse posterior #
                #####################

                min( summary(fit)$summary[, 'n_eff'] )	
                re <- rstan::extract(fit)		
                
                # extract hyperparams rho
                # _______________________

                ps <- c(0.025,0.25,0.5,0.75,0.975)
                .f <- function(x) quantile(x, probs = ps)
                tmp <- cbind( .f(re$rho_00),
                              .f(re$rho_10),
                              .f(re$rho_01),
                              .f(re$rho_11),
                              .f(re$alpha_00),
                              .f(re$alpha_10),
                              .f(re$alpha_01),
                              .f(re$alpha_11) )			
                colnames(tmp) <- c('rho_00','rho_10','rho_01','rho_11','alpha_00','alpha_10','alpha_01','alpha_11')
                rownames(tmp) <- c('CL','IL','M','IU','CU')

                tmp <- as.data.table(reshape2::melt(tmp))
                setnames(tmp, 'Var2', 'GP_hyper_par')
                tmp[, SEX:= as.integer(gsub('^([a-z]+)_([0-9])([0-9])','\\2',GP_hyper_par))]
                tmp[, LOC:= as.integer(gsub('^([a-z]+)_([0-9])([0-9])','\\3',GP_hyper_par))]
                tmp[, GP_hyper_par:= gsub('^([a-z]+)_([0-9])([0-9])','\\1',GP_hyper_par)]
                tmp <- dcast.data.table(tmp, LOC+SEX+GP_hyper_par~Var1, value.var='value')
                nspop.gp.pars <- merge(unique(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), tmp, by=c('SEX','LOC'))

                p <- ggplot(nspop.gp.pars, aes(x= paste0(GP_hyper_par, ' ', LOC_LABEL, ' ', SEX_LABEL))) +
                        geom_point(aes(y=M)) +
                        geom_errorbar(aes(ymin=CL, ymax=CU)) +
                        coord_flip() +
                        theme_bw() +
                        labs(x='GP hyperparameter\n', y='')

                filename <- paste0('220729f_notsuppAmongPop_gppars_round',round,'.pdf')
                filename2 <- gsub('pdf$','png',filename)
                cat('Saving', filename, '...\n')
                ggsave(p, file=file.path(outdir,filename), w=6, h=3)
                ggsave(p, file=file.path(outdir,filename2), w=6, h=3)
                
                
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
                
                p <- ggplot(nspop.by.age) + 		
                        geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2) +			
                        geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL)) +
                        scale_x_continuous( expand=c(0,0) ) + 
                        scale_y_continuous(label=scales:::percent) +
                        scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
                        facet_wrap(~LOC_LABEL, ncol=2) +
                        theme_bw() +
                        labs(x='\nage at visit (years)', 
                             y='population with unsuppressed viral load\n(95% credibility interval)\n', 
                             colour='gender', 
                             linetype='location')

                filename <- paste0('220729f_nsuppAmongPop_vs_age_by_gender_fishinland_round_',round,'.pdf')
                filename2 <- gsub('pdf$','png',filename)
                cat('Saving', filename, '...\n')
                ggsave(p, file=file.path(outdir,filename), w=6, h=5)
                ggsave(p, file=file.path(outdir,filename2), w=6, h=5)
                
                
                # extract basic not supp estimates
                # ________________________________

                DT[, list(N=sum(VLNS_N), P=sum(VLNS_N) / sum(N), N2=sum(N)-sum(VLNS_N), P2= 1-sum(VLNS_N) / sum(N)), by=c('LOC_LABEL','SEX_LABEL')]		
                
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
                     file=file.path(outdir,filename))


                # make table version suppressed
                # _____________________________

                nspop.by.age[, LABEL:= paste0(sprintf('%2.1f',M*100),' (',sprintf('%2.1f',CL*100),' - ',sprintf('%2.1f',CU*100),')') ]
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
                fwrite(dt, row.names=FALSE, file=file.path(outdir,filename))

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

vl.suppofpop.by.gender.loc.age.icar<- function(DT)
{
	
        # DT <- copy(dall)
	outdir <- file.path(out.dir)
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
                save(fit, file=file.path(outdir,filename))
                
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
                     file=file.path(outdir,filename))
                
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
                fwrite(dt, row.names=FALSE, file=file.path(outdir,filename))

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

vl.vlrunningprops.by.gender.loc.age<- function(DT)
{
        # TODO: 
        # Check by round

        # DT <- copy(dall)
	outdir <- file.path(out.dir)
        DT <- .preprocess.ds.oli(DT)

	tmp <- seq.int(min(DT$AGEYRS), max(DT$AGEYRS))
        tmp1 <- DT[, sort(unique(ROUND))]
	vla <- as.data.table(expand.grid(ROUND=tmp1, 
                                         FC=c('fishing','inland'),
                                         SEX=c('M','F'),
                                         AGEYRS=tmp))

        .f <- function(x,y) as.vector( unname ( binconf( sum(x), length(y) )))
	ans <- vla[, {		
				z <- which(DT$ROUND == ROUND, DT$FC==FC & DT$SEX==SEX & DT$AGEYRS<=(AGEYRS+2) & DT$AGEYRS>=(AGEYRS-1))				
                                z2 <- .f(DT$HIV_STATUS[z]==1, z)
                                z3 <- .f(DT$VLNS[z]==1, z)
                                z4 <- .f(DT$VLNS[z]==1, which(DT$HIV_STATUS[z]==1))
                                z5 <- .f(DT$ARVMED[z]==0 & DT$HIV_STATUS[z] & !is.na(DT$ARVMED), which(DT$hiv_status[z]==1 & !is.na(DT$arvmed[z]))) 

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
                                     PVLNSofHIV_CU= z4[3],
                                     PARVofHIV_MEAN= z5[1],
                                     PARVofHIV_CL= z5[2],
                                     PARVofHIV_CU= z5[3]
                                )				
        }, by=names(vla)]

	set(ans, NULL, 'SEX', ans[, factor(SEX, levels=c('M','F'), labels=c('men','women'))])
	
	# HIV prevalence
	# ______________

	p <- ggplot(ans) + 		
                geom_ribbon(aes(x=AGEYRS, ymin=PHIV_CL, ymax=PHIV_CU, group=interaction(SEX,FC)), alpha=0.2) +
                geom_line(aes(x=AGEYRS, y=PHIV_MEAN, colour=SEX)) +
                scale_x_continuous( expand=c(0,0) ) + 
                scale_y_continuous(label=scales:::percent) +
                scale_colour_manual(values=c('men'='royalblue3','women'='deeppink2')) +
                facet_wrap(~FC, ncol=2) +
                theme_bw() +
                labs(x='\nage at visit (years)', 
                     y='HIV prevalence (95% CI)\n', 
                     colour='gender', linetype='location')

	ggsave(p, file=file.path(outdir,'220729_hivprevalence_vs_age_by_gender_fishinland.pdf'), w=6, h=5)
	
	# HIV unsuppressed viral load
	# ___________________________

	p <- ggplot(ans) + 		
                geom_ribbon(aes(x=AGEYRS, ymin=PVLNS_CL, ymax=PVLNS_CU, group=interaction(SEX,FC)), alpha=0.2) +			
                geom_line(aes(x=AGEYRS, y=PVLNS_MEAN, colour=SEX)) +
                scale_x_continuous( expand=c(0,0) ) + 
                scale_y_continuous(label=scales:::percent) +
                scale_colour_manual(values=c('men'='royalblue3','women'='deeppink2')) +
                facet_wrap(~FC, ncol=2) +
                theme_bw() +
                labs(x='\nage at visit (years)', 
                     y='proportion unsuppressed HIV (95% CI)\n', 
                     colour='gender', 
                     linetype='location')

        ggsave(p, file=file.path(outdir,'220729_hivnotsupp_vs_age_by_gender_fishinland.pdf'), w=6, h=5)

	
	# HIV unsuppressed viral load among HIV+
	# ______________________________________

	p <- ggplot(ans) + 		
                geom_ribbon(aes(x=AGEYRS, ymin=PVLNSofHIV_CL, ymax=PVLNSofHIV_CU, group=interaction(SEX,FC)), alpha=0.2) +
                geom_line(aes(x=AGEYRS, y=PARVofHIV_MEAN, colour=SEX), linetype='dotted') +
                geom_line(aes(x=AGEYRS, y=PVLNSofHIV_MEAN, colour=SEX)) +
                scale_x_continuous( expand=c(0,0) ) + 
                scale_y_continuous(label=scales:::percent) +
                scale_colour_manual(values=c('men'='royalblue3','women'='deeppink2')) +
                facet_wrap(~FC, ncol=2) +
                theme_bw() +
                labs(x='\nage at visit (years)', 
                     y='proportion unsuppressed HIV among infected (95% CI)\n', 
                     colour='gender', 
                     linetype='location')

        ggsave(p, file=file.path(outdir,'220729_hivnotsuppofhiv_vs_age_by_gender_fishinland.pdf'), w=6, h=5)
	
	
	# write results to file
	# _____________________

	setkey(ans, FC, SEX, AGEYRS)
	ans[, PHIV_L:= paste0( round(PHIV_MEAN*100, d=1),' [', round(PHIV_CL*100, d=1),'-', round(PHIV_CU*100, d=1),']' )]
	ans[, PVLNS_L:= paste0( round(PVLNS_MEAN*100, d=1),' [', round(PVLNS_CL*100, d=1),'-', round(PVLNS_CU*100, d=1),']' )]
	ans[, PVLNSofHIV_L:= paste0( round(PVLNSofHIV_MEAN*100, d=1),' [', round(PVLNSofHIV_CL*100, d=1),'-', round(PVLNSofHIV_CU*100, d=1),']' )]
	write.csv(ans, file=file.path(outdir,'220729_keystats_by_age_gender_fishinland.csv'))
}

vl.vlrunningmean.by.gender.loc.age<- function()
{
	require(Hmisc)
	require(data.table)
	require(ggplot2)
	VL_DETECTABLE <- 4e2
	VIREMIC_VIRAL_LOAD <- 1e3
	
	
	require(data.table)
	prjdir	<- '~/Box/OR_Work/2018/2018_RakaiViralLoad'
	infile	<- file.path(prjdir,'data','191101_data_round17_vl_gps.rda')
	load( infile )
	setalloccol(ds)
	
	# remove HIV+ individuals with missing VLs
	ds <- subset(ds, HIV_STATUS==0 | HIV_AND_VL==1)
	
	# define VL_COPIES for uninfected
	set(ds, NULL, 'VLC', ds$VL_COPIES)
	set(ds, ds[,which(HIV_STATUS==0)], 'VLC', 0)
	
	# define undetectable VL (machine-undetectable)
	# define suppressed VL (according to WHO criteria)	
	set(ds, NULL, 'VLU', ds[, as.integer(VLC<VL_DETECTABLE)])
	set(ds, NULL, 'VLS', ds[, as.integer(VLC<VIREMIC_VIRAL_LOAD)])
	set(ds, NULL, 'VLD', ds[, as.integer(VLC>=VL_DETECTABLE)])
	set(ds, NULL, 'VLNS', ds[, as.integer(VLC>=VIREMIC_VIRAL_LOAD)])
	set(ds, NULL, 'HIV_AND_VLD', ds[, as.integer(VLD==1 & HIV_AND_VL==1)])
	
	# reset VLC below machine detectable to 0
	set(ds, ds[, which(HIV_AND_VL==1 & VLU==1)], 'VLC', 0)
	setkey(ds, FC, SEX, AGEYRS)
	
	tmp <- seq.int(min(ds$AGEYRS), max(ds$AGEYRS))
	vla <- as.data.table(expand.grid(FC=c('fishing','inland'), SEX=c('M','F'), AGEYRS=tmp))
	ans <- vla[, {
				z <- which(ds$FC==FC & ds$SEX==SEX & ds$AGEYRS<=(AGEYRS+2) & ds$AGEYRS>=(AGEYRS-2))
				z2 <- mean( ds$VLC[z] )
				z3 <- sd(ds$VLC[z])/sqrt(length(z))
				list(N= length(z),
					VLCM_M= z2,
					VLCM_CL= z2-1.96*z3,
					VLCM_CU= z2+1.96*z3
					)
			}, by=c('FC','SEX','AGEYRS')]
	set(ans, NULL, 'SEX', ans[, factor(SEX, levels=c('M','F'), labels=c('men','women'))])
	
	p <- ggplot(ans) + 
                #geom_errorbar(aes(x=AGEYRS, ymin=VLCM_CL, ymax=VLCM_CU)) +		
                geom_ribbon(aes(x=AGEYRS, ymin=VLCM_CL, ymax=VLCM_CU, group=interaction(SEX,FC)), alpha=0.2) +
                geom_hline(yintercept=1e3) +
                geom_line(aes(x=AGEYRS, y=VLCM_M, colour=SEX)) +
                scale_x_continuous( expand=c(0,0) ) + 
                scale_y_continuous() +
                scale_colour_manual(values=c('men'='royalblue3','women'='deeppink2')) +
                facet_wrap(~FC, ncol=2, scales='free_y') +
                theme_bw() +
                labs(x='\nage at visit (years)', 
                     y='mean viral load (95% CI)\n', 
                     colour='gender', linetype='location')

	ggsave(p, file=file.path(prjdir,'results_200220','200220_vlmean_vs_age_by_gender_fishinland.pdf'), w=8, h=5)
	
	
	#	loess mean below 0 for some age groups, not a good model
	#	sqrt transformation did not work, gave too low means
	ds[, VLCS:= sqrt(VLC)]
	vlclo <- ds[, loess(VLCS ~ AGEYRS, control=loess.control(trace.hat='approximate'))]	
	ans <- subset(ds, select=c(FC, SEX, AGEYRS, VLC, VLCS))	
	ans[, VLCLO_M:= (vlclo$fitted)^2]
	ggplot(ans) + 
			geom_line(aes(x=AGEYRS, y=VLCLO_M)) +
			scale_x_continuous( expand=c(0,0) )	
	ans <- ds[, {
				vlclo <- loess(VLC ~ AGEYRS, control=loess.control(trace.hat='approximate'))
				list(	VLC= VLC,
						AGEYRS= AGEYRS,
						VLCLO_M= vlclo$fitted 
				)				
			}, by=c('FC','SEX')]	
	predict(vlclo, newdata=NULL, se=TRUE)
}

vl.vldistribution.by.gender.loc<- function()
{
	require(Hmisc)
	require(data.table)
	VL_DETECTABLE <- 4e2
	VIREMIC_VIRAL_LOAD <- 1e3
	
	
	require(data.table)
	prjdir	<- '~/Box/OR_Work/2018/2018_RakaiViralLoad'
	infile	<- file.path(prjdir,'data','191101_data_round17_vl_gps.rda')
	load( infile )
	setalloccol(ds)
	
	# remove HIV+ individuals with missing VLs
	ds <- subset(ds, HIV_STATUS==0 | HIV_AND_VL==1)
	
	# define VL_COPIES for uninfected
	set(ds, NULL, 'VLC', ds$VL_COPIES)
	set(ds, ds[,which(HIV_STATUS==0)], 'VLC', 0)
	
	# define undetectable VL (machine-undetectable)
	# define suppressed VL (according to WHO criteria)	
	set(ds, NULL, 'VLU', ds[, as.integer(VLC<VL_DETECTABLE)])
	set(ds, NULL, 'VLS', ds[, as.integer(VLC<VIREMIC_VIRAL_LOAD)])
	set(ds, NULL, 'VLD', ds[, as.integer(VLC>=VL_DETECTABLE)])
	set(ds, NULL, 'VLNS', ds[, as.integer(VLC>=VIREMIC_VIRAL_LOAD)])
	set(ds, NULL, 'HIV_AND_VLD', ds[, as.integer(VLD==1 & HIV_AND_VL==1)])
	
	# reset VLC below machine detectable to 1e-6 (for plotting)
	set(ds, ds[, which(HIV_AND_VL==1 & VLU==1)], 'VLC', 1e-6)
	
	# plot proportion of population with viral load > x
	x <- seq(log(1),log(max(ds$VLC)), length.out=1e3)
	x <- c(0,exp(x))
	vld <- as.data.table(expand.grid(X=x, SEX=c('M','F'), FC=c('fishing','inland'), HIV_AND_VLD=c(0,1)))
	
	ans <- vld[, {
				n <- length(which(ds$SEX==SEX & ds$FC==FC & ds$HIV_AND_VLD>=HIV_AND_VLD))
				k <- length(which(ds$SEX==SEX & ds$FC==FC & ds$HIV_AND_VLD>=HIV_AND_VLD & X<ds$VLC))
				z<- as.vector( unname( binconf(k, n) ) )				
				list(N=n, K=k, P_M= z[1], P_CL=z[2], P_CU=z[3] )
			}, by=c('HIV_AND_VLD','FC','SEX','X')]
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

	ggsave(p, file=file.path(prjdir,'results_200220','200220_vldistribution_by_gender_fishinland.pdf'), w=9, h=5)
}

vl.keystats.by.gender.loc<- function()
{
	require(Hmisc)
	require(data.table)
	VL_DETECTABLE <- 4e2
	VIREMIC_VIRAL_LOAD <- 1e3
	
	
	require(data.table)
	prjdir	<- '~/Box/OR_Work/2018/2018_RakaiViralLoad'
	infile	<- file.path(prjdir,'data','191101_data_round17_vl_gps.rda')
	load( infile )
	setalloccol(ds)
	
	# remove HIV+ individuals with missing VLs
	ds <- subset(ds, HIV_STATUS==0 | HIV_AND_VL==1)
	
	# define VL_COPIES for uninfected
	set(ds, NULL, 'VLC', ds$VL_COPIES)
	set(ds, ds[,which(HIV_STATUS==0)], 'VLC', 0)
	
	# define undetectable VL (machine-undetectable)
	# define suppressed VL (according to WHO criteria)	
	set(ds, NULL, 'VLU', ds[, as.integer(VLC<VL_DETECTABLE)])
	set(ds, NULL, 'VLS', ds[, as.integer(VLC<VIREMIC_VIRAL_LOAD)])
	set(ds, NULL, 'VLD', ds[, as.integer(VLC>=VL_DETECTABLE)])
	set(ds, NULL, 'VLNS', ds[, as.integer(VLC>=VIREMIC_VIRAL_LOAD)])
	set(ds, NULL, 'HIV_AND_VLD', ds[, as.integer(VLD==1 & HIV_AND_VL==1)])
	ds[, table(VLD, VLNS)]
	#	   VLNS
	#VLD     0     1
  	#0 17577     0
  	#1    90   976
	#	--> there are only 90 individuals with VL in 4e2-1e3, so setting 4e2 or 1e3 is essentially the same
	
	# reset VLC below machine detectable to 0
	set(ds, ds[, which(HIV_AND_VL==1 & VLU==1)], 'VLC', 0)
	
	# entire population: 
	# mean viral load, proportion with DVL
	ds[, mean(VLC)]
	# 2290.494
	binconf( length(which(ds$VLD==1)), nrow(ds) )
	# PointEst   Lower      Upper
	# 0.05717964 0.05393704 0.06060469

	#
	# stratified by men/women inland/fishing
	# 
	ans <- ds[, {
			z<- as.vector( unname( binconf( length(which(HIV_STATUS==1)), length(HIV_STATUS) ) ) )
			z2<- as.vector( unname( binconf( length(which(VLNS==1)), length(VLNS) ) ) )
			z3<- as.vector( unname( binconf( length(which(VLNS==1)), length(which(HIV_STATUS==1)) ) ) )
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
	tmp <- ds[, {
				z<- as.vector( unname( binconf( length(which(HIV_STATUS==1)), length(HIV_STATUS) ) ) )
				z2<- as.vector( unname( binconf( length(which(VLNS==1)), length(VLNS) ) ) )
				z3<- as.vector( unname( binconf( length(which(VLNS==1)), length(which(HIV_STATUS==1)) ) ) )
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
	
	ans[, PHIV_L:= paste0( round(PHIV_MEAN*100, d=1),' [', round(PHIV_CL*100, d=1),'-', round(PHIV_CU*100, d=1),']' )]
	ans[, PVLNS_L:= paste0( round(PVLNS_MEAN*100, d=1),' [', round(PVLNS_CL*100, d=1),'-', round(PVLNS_CU*100, d=1),']' )]
	ans[, PVLNSofHIV_L:= paste0( round(PVLNSofHIV_MEAN*100, d=1),' [', round(PVLNSofHIV_CL*100, d=1),'-', round(PVLNSofHIV_CU*100, d=1),']' )]
	
	#FC SEX    N  PHIV_MEAN    PHIV_CL   PHIV_CU  PVLD_MEAN    PVLD_CL    PVLD_CU PVLDofHIV_MEAN PVLDofHIV_CL PVLDofHIV_CU  VLC_MEAN           PHIV_L           PVLD_L      PVLDofHIV_L
	#1: overall   F 9984 0.21764824 0.20966345 0.2258502 0.05348558 0.04924138 0.05807325      0.2457432    0.2281006    0.2642832 1376.0877   21.8 [21-22.6]    5.3 [4.9-5.8] 24.6 [22.8-26.4]
	#2: overall   M 8659 0.14943989 0.14208609 0.1571046 0.06143897 0.05657296 0.06669393      0.4111283    0.3846208    0.4381619 3344.8235 14.9 [14.2-15.7]    6.1 [5.7-6.7] 41.1 [38.5-43.8]
	#3: fishing   F 1938 0.44272446 0.42074507 0.4649305 0.11455108 0.10112790 0.12949930      0.2587413    0.2305585    0.2890747 3771.5119 44.3 [42.1-46.5] 11.5 [10.1-12.9] 25.9 [23.1-28.9]
	#4: fishing   M 2108 0.32542694 0.30575906 0.3457299 0.14753321 0.13303557 0.16331313      0.4533528    0.4164628    0.4907623 9052.8373 32.5 [30.6-34.6] 14.8 [13.3-16.3] 45.3 [41.6-49.1]
	#5:  inland   F 8046 0.16343525 0.15551676 0.1716750 0.03877703 0.03477391 0.04322036      0.2372624    0.2150559    0.2609994  799.1138 16.3 [15.6-17.2]    3.9 [3.5-4.3] 23.7 [21.5-26.1]
	#6:  inland   M 6551 0.09281026 0.08602036 0.1000774 0.03373531 0.02962926 0.03838786      0.3634868    0.3262210    0.4024669 1508.0821     9.3 [8.6-10]      3.4 [3-3.8] 36.3 [32.6-40.2]
	
	write.csv(ans, file=file.path(prjdir,'results_200220','200220_keystats_by_gender_fishinland.csv'))
}

vl.vlratio.by.loc<- function()
{
	require(Hmisc)
	require(data.table)
	VL_DETECTABLE <- 4e2
	VIREMIC_VIRAL_LOAD <- 1e3
	
	
	require(data.table)
	prjdir	<- '~/Box/OR_Work/2018/2018_RakaiViralLoad'
	infile	<- file.path(prjdir,'data','191101_data_round17_vl_gps.rda')
	load( infile )
	setalloccol(ds)
	
	# remove HIV+ individuals with missing VLs
	ds <- subset(ds, HIV_STATUS==0 | HIV_AND_VL==1)
	
	# define VL_COPIES for uninfected
	set(ds, NULL, 'VLC', ds$VL_COPIES)
	set(ds, ds[,which(HIV_STATUS==0)], 'VLC', 0)
	
	# define undetectable VL (machine-undetectable)
	# define suppressed VL (according to WHO criteria)	
	set(ds, NULL, 'VLU', ds[, as.integer(VLC<VL_DETECTABLE)])
	set(ds, NULL, 'VLS', ds[, as.integer(VLC<VIREMIC_VIRAL_LOAD)])
	set(ds, NULL, 'VLD', ds[, as.integer(VLC>=VL_DETECTABLE)])
	set(ds, NULL, 'VLNS', ds[, as.integer(VLC>=VIREMIC_VIRAL_LOAD)])
	set(ds, NULL, 'HIV_AND_VLD', ds[, as.integer(VLD==1 & HIV_AND_VL==1)])
	ds[, table(VLD, VLNS)]
	
	# reset VLC below machine detectable to 0
	set(ds, ds[, which(HIV_AND_VL==1 & VLU==1)], 'VLC', 0)
	
	
	ans	<- as.data.table(expand.grid(BS= 1:1e3, FC=c('fishing','inland')))
	set.seed(42)
	ans <- ans[, {
				zm <- which(ds$FC==FC & ds$SEX=='M')
				zf <- which(ds$FC==FC & ds$SEX=='F')
				zm <- sample(zm, length(zm), replace=TRUE)
				zf <- sample(zf, length(zf), replace=TRUE)
				list(VLCM_M=mean(ds$VLC[zm]), VLCM_F=mean(ds$VLC[zf])) 
			}, by=c('FC','BS')]
	ans[, VLCR:= VLCM_M/VLCM_F]	
	ans <- ans[, list( V=quantile(VLCR, prob=c(0.5, 0.025, 0.975)),
				P= c('M','CL','CU')
				), by=c('FC')]
	ans <- dcast.data.table(ans, FC~P, value.var='V')
	
	#
	# stratified by men/women inland/fishing
	# 
	ans <- ds[, {
				z<- as.vector( unname( binconf( length(which(HIV_STATUS==1)), length(HIV_STATUS) ) ) )
				z2<- as.vector( unname( binconf( length(which(VLNS==1)), length(VLNS) ) ) )
				z3<- as.vector( unname( binconf( length(which(VLNS==1)), length(which(HIV_STATUS==1)) ) ) )
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
	tmp <- ds[, {
				z<- as.vector( unname( binconf( length(which(HIV_STATUS==1)), length(HIV_STATUS) ) ) )
				z2<- as.vector( unname( binconf( length(which(VLNS==1)), length(VLNS) ) ) )
				z3<- as.vector( unname( binconf( length(which(VLNS==1)), length(which(HIV_STATUS==1)) ) ) )
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
	
	# drop few infecteds with missing VL
	ds <- subset(ds, HIV_STATUS==0 | (HIV_STATUS==1 & !is.na(VL_COPIES)) )
	# set VL for uninfected to 0, and VL with undetectable VL to 0
	set(ds, ds[, which(HIV_STATUS==0)], 'VL_COPIES', 0)
	set(ds, ds[, which(HIV_STATUS==1 & VL_UNDETECTABLE==1)], 'VL_COPIES', 0)
	
	# 
	# calculate proportion with VL > x among participants
	
	# do general by as characters
	# then determine sort index
	# then calculate empirical quantile
	ds <- ds[order(SEX,VL_COPIES),]
	ds[VL]
	
	ds[, sort(unique(VL_COPIES))]
	#dv <- data.table(VL:= )
}

vl.get.data.round17<- function()
{
	require(data.table)
	prjdir	<- '~/Box/OR_Work/2018/2018_RakaiViralLoad'
	infile	<- 'data_raw/ViralLoad_Data_Pangea_Ratmann.rda'
	load( file.path(prjdir, infile) )
	
	#
	# subset to survey round 17
	#
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
	
	tmp		<- ds[, list(		HIV_POS= length(which(HIV_STATUS==1)), 
					HIV_NEG= length(which(HIV_STATUS==0)),  
					HIV_PREV= length(which(HIV_STATUS==1))/length(HIV_STATUS)
					), by='COMM_NUM']
			
	thr	<- 1e3		
	tmp2	<- subset(ds, HIV_STATUS==1)[, list(		VL_D= length(which(VL_COPIES>thr)), 
					VL_U= length(which(VL_COPIES<=thr)),  
					VL_DP= length(which(VL_COPIES>thr))/length(VL_COPIES)
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

make.map.190129	<- function()
{
	require(data.table)
	require(rgdal)
	require(rgeos)
	library(raster)
	require(RColorBrewer) #Map colours
	
	# load data
	infile	<- '~/Box Sync/OR_Work/2018/2018_RakaiViralLoad/data/merged_round17_vl_gps.rda'
	load(infile)
		
	#convert the data into a data table
	dt<- as.data.table(ds)
	dt<- dt[,.(RCCS_STUDYID, SEX, AGEYRS, HIV_STATUS, LATITUDE_JITTER, LONGITUDE_JITTER, VL_COPIES, VL_UNDETECTABLE)]
	#set the NA VL to 0
	dt[is.na(VL_COPIES), VL_COPIES:=0]
	dt[,VL_DETECTABLE := as.numeric(VL_COPIES>=1000)]
	dt[,RCCS_STUDYID2:= seq_len(nrow(dt)) ]
		
	#################################################### load in maps
	# Load in Uganda Shape files 
	uganda1<-raster::getData('GADM',country="UGA",level=1)# Admin unit 1
	uganda3<- raster::getData('GADM', country='UGA', level=3)
	rakai1<-subset(uganda1, NAME_1=="Rakai")
	rakai3<- subset(uganda3, NAME_1=="Rakai")
	masaka1<-subset(uganda1, NAME_1=="Masaka")
	# Create a smaller Rakai for plotting (not current Rakai region no longer includes kabula subdistrict 3)
	#minirak<-rakai3[which(rakai3$NAME_2!="Kabula" | rakai3$NAME_3=="Lyantonde Tc" | rakai3$NAME_3=="Lyantonde"),]
	minirak<-rakai3[which(rakai3$NAME_2!="Kabula" | rakai3$NAME_3=="Lyantonde Tc"),]
	
	####################################################### Convert the data to meters
	#set the coordinates of the data
	coordinates(dt)<- ~ LONGITUDE_JITTER+LATITUDE_JITTER
	#set coordinate system to match uganda files
	proj4string(dt) <- proj4string(uganda1)
	
	#convert to m in order to build a 30x30m grid
	newcrs <- CRS("+proj=robin +datum=WGS84")
	dtnew<- spTransform(dt, newcrs)
	rakai1trans<- spTransform(rakai1, newcrs)
	miniraktrans<- spTransform(minirak, newcrs)
	masaka1trans<- spTransform(masaka1, newcrs)
	
	###################################################### Build Grid
	#Combine rakai1trans and masaka1trans
	outline<- union(rakai1trans, masaka1trans)
	#find the extent of the data
	exnew<- extent(dtnew)
	#extent of the maps
	exmap<- extent(outline)
	
	#chose extent to cover all the data and rakai district
	
	#With a 30m grid, I think the same individuals are usually entering calculations for a large number of grid points
	#Do we really need a 30m grid? Why not 100m?

	grid<- raster(xmn=min(exnew[1], exmap[1]), xmx= exnew[2], ymn=exmap[3], ymx=exnew[4], res=100 )
	#grid[]<- 1:ncell(grid) #No longer needed
	
	# set the coordinate reference system to match
	proj4string(grid)<- proj4string(dtnew) 
	
	#restrict grid to map
	#gridmask<- mask(grid, outline) #Restrict the map after
	#plot(gridmask)
	
	#consider the grid points in a data frame
	id<- as.data.table(1:ncell(gridmask))
	setnames(id, "V1", "ID")
	griddf<- as.data.table(SpatialPoints(grid))
	griddf<- data.table(id, griddf)
	setnames(griddf, gsub('y','LAT_GRID',gsub('x','LONG_GRID',colnames(griddf))))
	
	bw			<- 3000
	bw2			<- bw*bw
	#require(mvtnorm)
	#dmvnorm( c(3.84,0) )	# ~ 9.996634e-05 
	threshold	<- bw*3.84 	# cut if density is < 1e-4
	threshold	<- threshold*threshold	# square the threshold, to avoid sqrt calculations in loop 		
	norm.const	<- 1/(2*pi*bw2)
	
	tmp			<- griddf[1:1e4,]
	anst<- system.time({
		ans	<- tmp[, {
					z1	<- LONG_GRID - dtnew@coords[,'LONGITUDE_JITTER']
					z2	<- LAT_GRID - dtnew@coords[,'LATITUDE_JITTER']
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
					list( 	HIV_STATUS_MEAN=mean( dtnew@data$HIV_STATUS[z2] ),				#no weighting by distance
						HIV_STATUS_KERNEL=sum( dtnew@data$HIV_STATUS[z2]*w )/sum(w),		#Gaussian kernel
						VL_COPIES_KERNEL_GEOMMEAN = exp(sum(w*log(dtnew@data$VL_COPIES[z2]+1))/sum(w))-1, #Geometric Mean Kernel
                                                VL_DETECTABLE_KERNEL = sum( dtnew@data$VL_DETECTABLE[z2]*w)/sum(w) #Detectable Prevelance
      	)
				}, by=c('ID','LONG_GRID','LAT_GRID')]
	})
	 grid[]<- ans[, VL_DETECTABLE_KERNEL]
  	gridmask<- mask(grid, outline)
  	#Breaks chosen by looking at data - need refining
  	plot(gridmask, breaks = c(0, 0.025, 0.05, 0.075, 0.1, 0.5), 
       		col=brewer.pal(11, "RdYlGn")[c(10,9,5,4,3)] ,  axes=FALSE, box=FALSE, ylim= c(exmap[3],-6000), legend=FALSE)
  	plot(outline, add=TRUE)
  	par(xpd=TRUE)
  	legend("right", legend=c("0-2.5","2.5-5","5-7.5","7.5-10", ">10"),fill=brewer.pal(11, "RdYlGn")[c(10,9,5,4,3)],horiz = FALSE, inset=-0.175, title= "Prevelence of \n Detectable \n Viremia (%)",  cex=0.8, box.lty = 0)
  	grid[]<- ans[, VL_COPIES_KERNEL_GEOMMEAN]
  	gridmask<- mask(grid, outline)
  	plot(gridmask, breaks = c(0, 0.8, 1.5, 2.5, 3, 145), col=brewer.pal(11, "RdYlGn")[c(10,9,5,4,3)] ,  axes=FALSE, box=FALSE, ylim= c(exmap[3],-6000), legend=FALSE)
  	plot(outline, add=TRUE)
  	par(xpd=TRUE)
  	legend("right", legend=c("0-0.8","0.8-1.5","1.5-2.5","2.5-3", ">3"),fill=brewer.pal(11, "RdYlGn")[c(10,9,5,4,3)],horiz = FALSE, inset=-0.175, title= "Geometric Mean \n VL (Copies/ml)",  cex=0.8, box.lty = 0)
}


