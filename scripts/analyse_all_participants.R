# AIMS:
# - adapt older code from Oli to multi-round, longitudinal settings
# TODO: add ARVMED to dall in preprocess_data.R

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
n.cores <- min(4, parallel::detectCores()-1 )
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
vl.meanviralload.by.gender.loc.age.icar(dall)

# Estimate suppofinfected (whatever that means)
# ________________________
vl.suppofinfected.by.gender.loc.age.icar(dall)


vl.suppofinfected.by.gender.loc.age.gp<- function()
{
	
        # DT <- copy(dall)
	outdir <- file.path(out.dir)
        DT <- .preprocess.ds.oli(DT)
	require(data.table)
	
	tmp <- seq.int(min(DT$AGEYRS), max(DT$AGEYRS))

	tmp <- seq.int(min(DT$AGEYRS), max(DT$AGEYRS))
        tmp1 <- DT[, sort(unique(ROUND))]
	vla <- as.data.table(expand.grid(ROUND=tmp1,
                                         FC=c('fishing','inland'),
                                         SEX=c('M','F'),
                                         AGEYRS=tmp))
	vla <- vla[, {		
				z <- which(DT$FC==FC & DT$SEX==SEX & DT$AGEYRS==AGEYRS)	
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
	
	stan.code2 <- ""
		
	stan.model <- stan_model(model_name= 'gp_all',model_code = gsub('\t',' ',stan.code2))	
	stan.data <- list()
	stan.data$x_predict <- seq(vla[, min(AGE_LABEL)], vla[, max(AGE_LABEL)+1], 0.5)
	stan.data$N_predict <- length(stan.data$x_predict)
	stan.data$observed_idx <- which(stan.data$x_predict%%1==0.5)
	stan.data$N_observed <- length(stan.data$observed_idx)
	stan.data$y_observed_00 <- vla[SEX==0 & LOC==0, HIV_N-VLNS_N]
	stan.data$y_observed_10 <- vla[SEX==1 & LOC==0, HIV_N-VLNS_N]
	stan.data$y_observed_01 <- vla[SEX==0 & LOC==1, HIV_N-VLNS_N]
	stan.data$y_observed_11 <- vla[SEX==1 & LOC==1, HIV_N-VLNS_N]
	stan.data$total_observed_00 <- vla[SEX==0 & LOC==0, HIV_N]
	stan.data$total_observed_10 <- vla[SEX==1 & LOC==0, HIV_N]
	stan.data$total_observed_01 <- vla[SEX==0 & LOC==1, HIV_N]
	stan.data$total_observed_11 <- vla[SEX==1 & LOC==1, HIV_N]
	stan.data$rho_hyper_par_00 <- diff(range(stan.data$x_predict))/3
	stan.data$rho_hyper_par_10 <- diff(range(stan.data$x_predict))/3
	stan.data$rho_hyper_par_01 <- diff(range(stan.data$x_predict))/3
	stan.data$rho_hyper_par_11 <- diff(range(stan.data$x_predict))/3
	stan.data$alpha_hyper_par_00 <- 2
	stan.data$alpha_hyper_par_10 <- 2
	stan.data$alpha_hyper_par_01 <- 2
	stan.data$alpha_hyper_par_11 <- 2
        
	fit <- sampling(stan.model, data=stan.data, iter=10e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999))
	save(fit, file=file.path(outdir, "200428f_notsuppAmongInfected_gp_stanfit.rda"))
	
	
	#
	#	compare to self-report
	#
	stan.data <- list()
	stan.data$x_predict <- seq(vla[, min(AGE_LABEL)], vla[, max(AGE_LABEL)+1], 0.5)
	stan.data$N_predict <- length(stan.data$x_predict)
	stan.data$observed_idx <- which(stan.data$x_predict%%1==0.5)
	stan.data$N_observed <- length(stan.data$observed_idx)
	stan.data$y_observed_00 <- vla[SEX==0 & LOC==0, HIV_N-ARV_N]
	stan.data$y_observed_10 <- vla[SEX==1 & LOC==0, HIV_N-ARV_N]
	stan.data$y_observed_01 <- vla[SEX==0 & LOC==1, HIV_N-ARV_N]
	stan.data$y_observed_11 <- vla[SEX==1 & LOC==1, HIV_N-ARV_N]
	stan.data$total_observed_00 <- vla[SEX==0 & LOC==0, HIV_N]
	stan.data$total_observed_10 <- vla[SEX==1 & LOC==0, HIV_N]
	stan.data$total_observed_01 <- vla[SEX==0 & LOC==1, HIV_N]
	stan.data$total_observed_11 <- vla[SEX==1 & LOC==1, HIV_N]
	stan.data$rho_hyper_par_00 <- diff(range(stan.data$x_predict))/3
	stan.data$rho_hyper_par_10 <- diff(range(stan.data$x_predict))/3
	stan.data$rho_hyper_par_01 <- diff(range(stan.data$x_predict))/3
	stan.data$rho_hyper_par_11 <- diff(range(stan.data$x_predict))/3
	stan.data$alpha_hyper_par_00 <- 2
	stan.data$alpha_hyper_par_10 <- 2
	stan.data$alpha_hyper_par_01 <- 2
	stan.data$alpha_hyper_par_11 <- 2
	fit2 <- sampling(stan.model, data=stan.data, iter=10e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999))
	save(fit2, file=file.path(outdir, "200428f_notARVAmongInfected_icar_stanfit.rda"))
	
	
	re <- rstan::extract(fit)
	re2 <- rstan::extract(fit2)
	ps <- c(0.025,0.5,0.975)
	tmp <- summary(fit)$summary
	tmp[grepl('^p_predict_',rownames(tmp)),]
	
	#
	#	extract hyperparams rho
	ps <- c(0.025,0.25,0.5,0.75,0.975)
	tmp <- cbind( quantile(re$rho_00, probs=ps),
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
	nsinf.gp.pars <- merge(unique(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), tmp, by=c('SEX','LOC'))
	ggplot(nsinf.gp.pars, aes(x= paste0(GP_hyper_par, ' ', LOC_LABEL, ' ', SEX_LABEL))) +
			geom_point(aes(y=M)) +
			geom_errorbar(aes(ymin=CL, ymax=CU)) +
			coord_flip() +
			theme_bw() +
			labs(x='GP hyperparameter\n', y='')
	ggsave(file=file.path(prjdir,'results_200220','200428f_notsuppAmongInfected_gppars.pdf'), w=6, h=3)
	
	
	#
	#	make prevalence plot by age
	tmp <- cbind( apply(re$p_predict_00, 2, quantile, probs=ps),
			apply(re$p_predict_10, 2, quantile, probs=ps),
			apply(re$p_predict_01, 2, quantile, probs=ps),
			apply(re$p_predict_11, 2, quantile, probs=ps)
			)
	rownames(tmp) <- c('CL','IL','M','IU','CU')
	tmp <- as.data.table(reshape2::melt(tmp))
	nsinf.by.age <- dcast.data.table(tmp, Var2~Var1, value.var='value')
	tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
	nsinf.by.age <- cbind(tmp, nsinf.by.age) 
	nsinf.by.age <- merge(unique(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), nsinf.by.age, by=c('SEX','LOC'))
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
	nainf.by.age <- merge(unique(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), nainf.by.age, by=c('SEX','LOC'))
	nainf.by.age[, STAT:='VLNA']
	tmp <- subset(nainf.by.age, select=c(SEX, LOC, AGE_LABEL, M, CL, CU))
	setnames(tmp, c('M','CL','CU'), c('M2','CL2','CU2'))
	tmp <- merge(nsinf.by.age, tmp, by=c('SEX','LOC','AGE_LABEL'))			
	ggplot(tmp) + 		
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
	ggsave(file=file.path(prjdir,'results_200220','200428f_suppAmongInfected_vs_age_by_gender_fishinland_stan.pdf'), w=6, h=5)		
	ggplot(tmp) + 		
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
	ggsave(file=file.path(prjdir,'results_200220','200428f_suppAmongInfected_vs_age_by_gender_fishinland_stan_v2.pdf'), w=6, h=5)	
	tmp <- rbind(nsinf.by.age, nainf.by.age, fill=TRUE)
	ggplot(tmp) + 		
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
	ggsave(file=file.path(prjdir,'results_200220','200428f_suppAmongInfected_vs_age_by_gender_fishinland_stan_v3.pdf'), w=9, h=8)
	
	
	#	extract basic not supp estimates
	vla[, list(N=sum(VLNS_N), P=sum(VLNS_N) / sum(HIV_N), N2=sum(HIV_N)-sum(VLNS_N), P2= 1-sum(VLNS_N) / sum(HIV_N)), by=c('LOC_LABEL','SEX_LABEL')]
	#   LOC_LABEL SEX_LABEL   N         P   N2        P2
	#1:   fishing         M 296 0.4314869  390 0.5685131
	#2:    inland         M 203 0.3338816  405 0.6661184
	#3:   fishing         F 198 0.2307692  660 0.7692308
	#4:    inland         F 279 0.2121673 1036 0.7878327
	

	ps <- c(0.025,0.5,0.975)
	tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
	rp <- as.data.table(reshape2::melt( tmp ))
	setnames(rp, 1:3, c('iterations','ROW_ID','P')) 
	tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
	tmp[, ROW_ID:= seq_len(nrow(tmp))]
	tmp <- merge(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL,N)), tmp, by=c('SEX','LOC','AGE_LABEL'), all=TRUE)
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
	# LOC_LABEL SEX_LABEL        CL        CU         M                    LABEL
	#	fishing         M 0.4562322 0.5536721 0.5018376 50.18% (45.62% - 55.37%)
	#2: fishing         F 0.7134543 0.7781516 0.7469707  74.7% (71.35% - 77.82%)
	#3: inland         M 0.4910456 0.6555350 0.5703211  57.03% (49.1% - 65.55%)
	#4: inland         F 0.6939593 0.7646014 0.7308165  73.08% (69.4% - 76.46%)



	#	extract risk ratio of suppressed VL female:male and male:female
	tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
	rp <- as.data.table(reshape2::melt( tmp ))
	setnames(rp, 1:3, c('iterations','ROW_ID','P')) 
	tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
	tmp[, ROW_ID:= seq_len(nrow(tmp))]
	tmp <- merge(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL,N)), tmp, by=c('SEX','LOC','AGE_LABEL'), all=TRUE)
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
	#   LOC_LABEL variable        CL        CU         M              LABEL
	#1:   fishing    PR_FM 1.3365856 1.6508029 1.4873708 1.49 (1.34 - 1.65)
	#2:   fishing    PR_MF 0.6057659 0.7481750 0.6723273 0.67 (0.61 - 0.75)
	#3:    inland    PR_FM 1.1015949 1.4994209 1.2811266   1.28 (1.1 - 1.5)
	#4:    inland    PR_MF 0.6669241 0.9077748 0.7805630 0.78 (0.67 - 0.91)
	
	
	#	extract risk ratio of unsuppressed VL female:male and male:female by age
	tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
	rp <- as.data.table(reshape2::melt( tmp ))
	setnames(rp, 1:3, c('iterations','ROW_ID','P')) 
	tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
	tmp[, ROW_ID:= seq_len(nrow(tmp))]
	tmp <- merge(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL,N)), tmp, by=c('SEX','LOC','AGE_LABEL'), all=TRUE)
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
	
	
	#	extract if difference in female:male risk ratio of unsuppressed VL is different in fishing vs inland
	tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
	rp <- as.data.table(reshape2::melt( tmp ))
	setnames(rp, 1:3, c('iterations','ROW_ID','P')) 
	tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
	tmp[, ROW_ID:= seq_len(nrow(tmp))]
	tmp <- merge(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL,N)), tmp, by=c('SEX','LOC','AGE_LABEL'), all=TRUE)
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
	#	            Q  P
	#1: 1: -0.05655754 CL
	#2:  0.20449222  M
	#3:  0.45032189 CU
	
	
	save(vla, re, re2, nainf.by.age, nsinf.by.age, nsinf.by.sex.loc, nsinf.by.loc, nsinf.ratio.by.loc.age, file=file.path(outdir, "200428f_suppAmongInfected.rda"))
	
	
	#	make table version suppressed
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
	write.csv(dt, row.names=FALSE, file=file.path(outdir, "200428f_suppamonginfected.csv"))	
}
