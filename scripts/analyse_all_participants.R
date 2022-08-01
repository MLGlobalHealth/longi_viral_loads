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


vl.meanviralload.by.gender.loc.age.gp<- function()
{
	require(Hmisc)
	require(data.table)
	require(ggplot2)
	require(rstan)
	VL_DETECTABLE <- 4e2
	VIREMIC_VIRAL_LOAD <- 1e3
		
	prjdir	<- '~/Box/OR_Work/2018/2018_RakaiViralLoad'
	infile	<- file.path(prjdir,'data','191101_data_round17_vl_gps.rda')
	outdir <- file.path(prjdir,'results_200220')
	load( infile )
	setalloccol(ds)
	
	# remove HIV+ individuals with missing VLs
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
	setkey(DT, FC, SEX, AGEYRS)	
	
	tmp <- seq.int(min(DT$AGEYRS), max(DT$AGEYRS))
	vla <- as.data.table(expand.grid(FC=c('fishing','inland'), SEX=c('M','F'), AGEYRS=tmp))
	vla <- vla[, {		
				z <- which(DT$FC==FC & DT$SEX==SEX & DT$AGEYRS==AGEYRS)
				z2 <- which(DT$FC==FC & DT$SEX==SEX & DT$AGEYRS==AGEYRS & DT$VLC>0)
				list(	N= length(z),
						VL_MEAN= mean(DT$VLC[z]),
						VL_MEAN_SD= sd(DT$VLC[z]) / sqrt(length(z)),
						VL_SD= sd(DT$VLC[z]),						
						VLNZ_N= length(z2),
						VLNZ_MEAN= mean(log(DT$VLC[z2])),
						VLNZ_MEAN_SD= sd(log(DT$VLC[z2])) / sqrt(length(z2)),
						VLNZ_S2= var(log(DT$VLC[z2])),
						HIV_N= length(which(DT$HIV_STATUS[z]==1)),
						VLNS_N= length(which(DT$VLNS[z]==1)),
						ARV_N= length(which(DT$ARVMED[z]==0 & DT$HIV_STATUS[z]==1 & !is.na(DT$ARVMED[z])))
				)				
			}, by=c('FC','SEX','AGEYRS')]
	setnames(vla, c('FC','SEX','AGEYRS'), c('LOC_LABEL','SEX_LABEL','AGE_LABEL'))
	vla[, LOC:= as.integer(LOC_LABEL=='fishing')]
	vla[, SEX:= as.integer(SEX_LABEL=='M')]
	vla[, AGE:= AGE_LABEL-14L]
	vla[, ROW_ID:= seq_len(nrow(vla))]
	
	
	if(0)
	{			
		stan.code <- "
data{	
	int<lower=1> N_predict;
	real x_predict[N_predict];
	int<lower=1> N_observed;
	int<lower=1, upper=N_predict> observed_idx[N_observed];
	int zero_observed[N_observed];
	int total_observed[N_observed];
	vector<lower=0>[N_observed] mean_observed;
	vector<lower=0>[N_observed] meansd_observed;
	real<lower=0> zero_rho_hyper_par;
	real<lower=0> zero_alpha_hyper_par;
	real<lower=0> mean_rho_hyper_par;
	real<lower=0> mean_alpha_hyper_par;
}
			
parameters {
	real<lower=0> zero_rho;
	real<lower=0> zero_alpha;
	real zero_base;
	vector[N_predict] zero_f_tilde;
	real<lower=0> mean_rho;
	real<lower=0> mean_alpha;			
	real mean_base;
	vector[N_predict] mean_f_tilde;
}
			
transformed parameters {
	matrix[N_predict, N_predict] L_cov;
	vector[N_predict] logit_zero_p_predict;
	vector[N_predict] zero_p_predict;
	vector[N_predict] mean_predict;
	// GP on zeros
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, zero_alpha, zero_rho) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_zero_p_predict = zero_base + L_cov * zero_f_tilde;
	zero_p_predict = inv_logit(logit_zero_p_predict);
	// GP on means
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, mean_alpha, mean_rho) + diag_matrix(rep_vector(1e-10, N_predict)));
	mean_predict = mean_base + L_cov * mean_f_tilde;

}
			
model {
	zero_rho ~ normal(0, zero_rho_hyper_par);
	zero_alpha ~ normal(0, zero_alpha_hyper_par);
	zero_base ~ normal( 0 , 100 );
	zero_f_tilde ~ normal(0, 1);
	
	mean_rho ~ normal(0, mean_rho_hyper_par);
	mean_alpha ~ normal(0, mean_alpha_hyper_par);
	mean_base ~ normal( 0 , 100 );
	mean_f_tilde ~ normal(0, 1);
	
	zero_observed ~ binomial_logit( total_observed, logit_zero_p_predict[observed_idx] );
	mean_observed ~ normal( (1-zero_p_predict[observed_idx]) .* mean_predict[observed_idx], meansd_observed + rep_vector(1e-10, N_observed));
}
			
generated quantities {
	vector[N_predict] zero_inflated_mean_predict = (1-zero_p_predict) .* mean_predict;  
}			
"
		stan.model <- stan_model(model_name= 'gp_one_zero_inflated',model_code = gsub('\t',' ',stan.code))
		vla2 <- subset(vla, SEX==0 & LOC==0)
		stan.data <- list()	
		stan.data$x_predict <- seq(vla2[, min(AGE_LABEL)], vla[, max(AGE_LABEL)+1], 0.5)
		stan.data$N_predict <- length(stan.data$x_predict)
		stan.data$observed_idx <- which(stan.data$x_predict%%1==0.5)
		stan.data$N_observed <- length(stan.data$observed_idx)
		stan.data$zero_observed <- vla2[,VLC_ZERO] 
		stan.data$total_observed <- vla2[,N]
		stan.data$mean_observed <- vla2[,VL_MEAN]
		stan.data$meansd_observed <- vla2[,VL_MEAN_SD]	
		stan.data$zero_rho_hyper_par <- diff(range(stan.data$x_predict))/3
		stan.data$zero_alpha_hyper_par <- 2
		stan.data$mean_rho_hyper_par <- diff(range(stan.data$x_predict))/3
		stan.data$mean_alpha_hyper_par <- 2
		fit <- sampling(stan.model, data=stan.data, iter=2e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999) )	
	}
	if(0)
	{
		stan.code2 <- "
data{	
	int<lower=1> N_predict;
	real x_predict[N_predict];
	int<lower=1> N_observed;
	int<lower=1, upper=N_predict> observed_idx[N_observed];
	vector<lower=0>[N_observed] y_observed;
	vector<lower=0>[N_observed] sd_observed;
	real<lower=0> rho_hyper_par;
	real<lower=0> alpha_hyper_par;
}
			
parameters {
	real<lower=0> rho;
	real<lower=0> alpha;
	real base;
	vector[N_predict] f_tilde;	
}
			
transformed parameters {
	matrix[N_predict, N_predict] L_cov;
	vector[N_predict] mean_predict;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha, rho) + diag_matrix(rep_vector(1e-10, N_predict)));
	mean_predict = base + L_cov * f_tilde;
}
			
model {
	rho ~ normal(0, rho_hyper_par);
	alpha ~ normal(0, alpha_hyper_par);	
	base ~ normal( 0 , 100 );
	f_tilde ~ normal(0, 1);
		
	y_observed ~ normal( mean_predict[observed_idx], sd_observed );
}				
"
		stan.model <- stan_model(model_name= 'gp_one_logscale',model_code = gsub('\t',' ',stan.code2))
		vla2 <- subset(vla, SEX==0 & LOC==0)
		stan.data <- list()	
		stan.data$x_predict <- seq(vla2[, min(AGE_LABEL)], vla[, max(AGE_LABEL)+1], 0.5)
		stan.data$N_predict <- length(stan.data$x_predict)
		stan.data$observed_idx <- which(stan.data$x_predict%%1==0.5)
		stan.data$observed_idx <- stan.data$observed_idx[ !is.na(vla2$VLNZ_MEAN_SD) ]
		stan.data$N_observed <- length(stan.data$observed_idx)
		stan.data$y_observed <- vla2[!is.na(VLNZ_MEAN_SD),VLNZ_MEAN] 
		stan.data$sd_observed <- vla2[!is.na(VLNZ_MEAN_SD),VLNZ_MEAN_SD]
		stan.data$rho_hyper_par <- diff(range(stan.data$x_predict))/3
		stan.data$alpha_hyper_par <- 2
		fit <- sampling(stan.model, data=stan.data, iter=2e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999) )
	}
	
	
	
	
	
	
	stan.code3 <- "
data{	
	int<lower=1> N_predict;
	real x_predict[N_predict];
	int<lower=1> N_observed;
	int<lower=1, upper=N_predict> observed_idx[N_observed];
	vector<lower=0>[N_observed] y_observed;
	vector<lower=0>[N_observed] s2_observed; // with denominator (n-1)
	vector<lower=1>[N_observed] n_observed;
	real<lower=0> m_rho_hyper_par_alpha;
	real<lower=0> m_rho_hyper_par_beta;
	real<lower=0> m_alpha_hyper_par;
	real<lower=0> s_rho_hyper_par_alpha;
	real<lower=0> s_rho_hyper_par_beta;
	real<lower=0> s_alpha_hyper_par;
}

transformed data{
	vector[N_observed] n_m1_observed;
	vector[N_observed] inv_sqrt_n_observed;
	vector[N_observed] inv_scaled_s2_observed;

	n_m1_observed = n_observed - rep_vector(1, N_observed);
	inv_sqrt_n_observed = inv(sqrt(n_observed));
	inv_scaled_s2_observed = inv(s2_observed .* n_m1_observed);
}
			
parameters {
	real<lower=0> m_rho;
	real<lower=0> m_alpha;
	real<lower=0> s_rho;
	real<lower=0> s_alpha;
	real m_base;
	real s_base;
	vector[N_predict] m_f_tilde;
	vector[N_predict] s_f_tilde;	
}
			
transformed parameters {
	matrix[N_predict, N_predict] L_cov;
	vector[N_predict] m_predict;
	vector[N_predict] s_predict;
	// GP on mean of lognormal
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, m_alpha, m_rho) + diag_matrix(rep_vector(1e-10, N_predict)));
	m_predict = m_base + L_cov * m_f_tilde;
	// GP on log sigma of lognormal
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, s_alpha, s_rho) + diag_matrix(rep_vector(1e-10, N_predict)));
	s_predict = exp( s_base + L_cov * s_f_tilde);
}
			
model {
	m_rho ~ inv_gamma(m_rho_hyper_par_alpha, m_rho_hyper_par_beta);
	m_alpha ~ normal(0, m_alpha_hyper_par);	
	m_base ~ normal( 0 , 10 );
	s_rho ~ inv_gamma(s_rho_hyper_par_alpha, s_rho_hyper_par_beta);
	s_alpha ~ normal(0, s_alpha_hyper_par);	
	s_base ~ normal( 0 , 10 );	
	m_f_tilde ~ normal(0, 1);
	s_f_tilde ~ normal(0, 1);

	target+= normal_lpdf( y_observed | m_predict[observed_idx], s_predict[observed_idx] .* inv_sqrt_n_observed);
	target+= inv_chi_square_lpdf( s_predict[observed_idx] .* s_predict[observed_idx] .* inv_scaled_s2_observed | n_m1_observed );  	
}

generated quantities {
	vector[N_predict] exp_m_predict = exp( m_predict + s_predict .* s_predict .* rep_vector(0.5, N_predict) );  
}			
"

	
	stan.model <- stan_model(model_name= 'gp_one_logscale',model_code = gsub('\t',' ',stan.code3))	
	vla.select <- vla[, which(SEX==1 & LOC==1)]
	vla.select <- vla[, which(SEX==1 & LOC==0)]
	vla.select <- vla[, which(SEX==0 & LOC==1)]
	stan.data <- list()	
	stan.data$x_predict <- seq(vla[vla.select, min(AGE_LABEL)], vla[vla.select, max(AGE_LABEL)+1], 0.5)
	stan.data$N_predict <- length(stan.data$x_predict)
	stan.data$observed_idx <- which(stan.data$x_predict%%1==0.5)
	stan.data$observed_idx <- stan.data$observed_idx[ vla[vla.select,][,!is.na(VLNZ_S2)] ]
	stan.data$N_observed <- length(stan.data$observed_idx)	
	stan.data$y_observed <- vla[vla.select,][!is.na(VLNZ_S2),VLNZ_MEAN]
	stopifnot( length(stan.data$y_observed)==stan.data$N_observed )
	stan.data$n_observed <- vla[vla.select,][!is.na(VLNZ_S2),VLNZ_N]
	stopifnot( length(stan.data$n_observed)==stan.data$N_observed )
	stan.data$s2_observed <- vla[vla.select,][!is.na(VLNZ_S2),VLNZ_S2]
	stopifnot( length(stan.data$s2_observed)==stan.data$N_observed )	
	cl.target <- diff(range(stan.data$x_predict))/10
	cu.target <- diff(range(stan.data$x_predict))
	inv.gamma.err <- function(pars){ abs(sum( qinvgamma(c(0.025,0.975), exp(pars[1]), exp(pars[2])) - c(cl.target,cu.target) )) }	
	tmp <- optim(log(c(8,30)), inv.gamma.err)
	inv.gamma.pars <- exp(tmp$par)
	stan.data$m_rho_hyper_par_alpha <- inv.gamma.pars[1]
	stan.data$m_rho_hyper_par_beta <- inv.gamma.pars[2]
	stan.data$m_alpha_hyper_par <- 2
	stan.data$s_rho_hyper_par_alpha <- inv.gamma.pars[1]
	stan.data$s_rho_hyper_par_beta <- inv.gamma.pars[2]
	stan.data$s_alpha_hyper_par <- 2
	fit <- sampling(stan.model, data=stan.data, iter=20e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999) )
	#save(fit, file=file.path(outdir, "200428f_hivmeans_gp_stanfit_10.rda"))
	save(fit, file=file.path(outdir, "200428f_hivmeans_gp_stanfit_01.rda"))

	stan.code4 <- "
data{	
	int<lower=1> N_predict;
	real x_predict[N_predict];
	int<lower=1> N_observed_00;
	int<lower=1> N_observed_10;
	int<lower=1> N_observed_01;
	int<lower=1> N_observed_11;
	int<lower=1, upper=N_predict> observed_idx_00[N_observed_00];
	int<lower=1, upper=N_predict> observed_idx_10[N_observed_10];
	int<lower=1, upper=N_predict> observed_idx_01[N_observed_01];
	int<lower=1, upper=N_predict> observed_idx_11[N_observed_11];
	vector<lower=0>[N_observed_00] y_observed_00;
	vector<lower=0>[N_observed_10] y_observed_10;
	vector<lower=0>[N_observed_01] y_observed_01;
	vector<lower=0>[N_observed_11] y_observed_11;
	vector<lower=0>[N_observed_00] s2_observed_00; // with denominator (n-1)
	vector<lower=0>[N_observed_10] s2_observed_10; // with denominator (n-1)
	vector<lower=0>[N_observed_01] s2_observed_01; // with denominator (n-1)
	vector<lower=0>[N_observed_11] s2_observed_11; // with denominator (n-1)
	vector<lower=1>[N_observed_00] n_observed_00;
	vector<lower=1>[N_observed_10] n_observed_10;
	vector<lower=1>[N_observed_01] n_observed_01;
	vector<lower=1>[N_observed_11] n_observed_11;
	real<lower=0> m_rho_hyper_par_00;
	real<lower=0> m_alpha_hyper_par_00;
	real<lower=0> s_rho_hyper_par_00;
	real<lower=0> s_alpha_hyper_par_00;
	real<lower=0> m_rho_hyper_par_10;
	real<lower=0> m_alpha_hyper_par_10;
	real<lower=0> s_rho_hyper_par_10;
	real<lower=0> s_alpha_hyper_par_10;
	real<lower=0> m_rho_hyper_par_01;
	real<lower=0> m_alpha_hyper_par_01;
	real<lower=0> s_rho_hyper_par_01;
	real<lower=0> s_alpha_hyper_par_01;
	real<lower=0> m_rho_hyper_par_11;
	real<lower=0> m_alpha_hyper_par_11;
	real<lower=0> s_rho_hyper_par_11;
	real<lower=0> s_alpha_hyper_par_11;
}

transformed data{
	vector[N_observed_00] n_m1_observed_00;
	vector[N_observed_00] inv_sqrt_n_observed_00;
	vector[N_observed_00] inv_scaled_s2_observed_00;
	vector[N_observed_10] n_m1_observed_10;
	vector[N_observed_10] inv_sqrt_n_observed_10;
	vector[N_observed_10] inv_scaled_s2_observed_10;
	vector[N_observed_01] n_m1_observed_01;
	vector[N_observed_01] inv_sqrt_n_observed_01;
	vector[N_observed_01] inv_scaled_s2_observed_01;
	vector[N_observed_11] n_m1_observed_11;
	vector[N_observed_11] inv_sqrt_n_observed_11;
	vector[N_observed_11] inv_scaled_s2_observed_11;

	n_m1_observed_00 = n_observed_00 - rep_vector(1, N_observed_00);
	inv_sqrt_n_observed_00 = inv(sqrt(n_observed_00));
	inv_scaled_s2_observed_00 = inv(s2_observed_00 .* n_m1_observed_00);
	n_m1_observed_10 = n_observed_10 - rep_vector(1, N_observed_10);
	inv_sqrt_n_observed_10 = inv(sqrt(n_observed_10));
	inv_scaled_s2_observed_10 = inv(s2_observed_10 .* n_m1_observed_10);
	n_m1_observed_01 = n_observed_01 - rep_vector(1, N_observed_01);
	inv_sqrt_n_observed_01 = inv(sqrt(n_observed_01));
	inv_scaled_s2_observed_01 = inv(s2_observed_01 .* n_m1_observed_01);
	n_m1_observed_11 = n_observed_11 - rep_vector(1, N_observed_11);
	inv_sqrt_n_observed_11 = inv(sqrt(n_observed_11));
	inv_scaled_s2_observed_11 = inv(s2_observed_11 .* n_m1_observed_11);
}
			
parameters {
	real<lower=0> m_rho_00;
	real<lower=0> m_alpha_00;
	real<lower=0> s_rho_00;
	real<lower=0> s_alpha_00;
	real m_base_00;
	real s_base_00;
	vector[N_predict] m_f_tilde_00;
	vector[N_predict] s_f_tilde_00;	
	real<lower=0> m_rho_10;
	real<lower=0> m_alpha_10;
	real<lower=0> s_rho_10;
	real<lower=0> s_alpha_10;
	real m_base_10;
	real s_base_10;
	vector[N_predict] m_f_tilde_10;
	vector[N_predict] s_f_tilde_10;	
	real<lower=0> m_rho_01;
	real<lower=0> m_alpha_01;
	real<lower=0> s_rho_01;
	real<lower=0> s_alpha_01;
	real m_base_01;
	real s_base_01;
	vector[N_predict] m_f_tilde_01;
	vector[N_predict] s_f_tilde_01;	
	real<lower=0> m_rho_11;
	real<lower=0> m_alpha_11;
	real<lower=0> s_rho_11;
	real<lower=0> s_alpha_11;
	real m_base_11;
	real s_base_11;
	vector[N_predict] m_f_tilde_11;
	vector[N_predict] s_f_tilde_11;	
}
			
transformed parameters {
	matrix[N_predict, N_predict] L_cov;
	vector[N_predict] m_predict_00;
	vector[N_predict] s_predict_00;
	vector[N_predict] m_predict_10;
	vector[N_predict] s_predict_10;
	vector[N_predict] m_predict_01;
	vector[N_predict] s_predict_01;
	vector[N_predict] m_predict_11;
	vector[N_predict] s_predict_11;
	// GP on mean of lognormal
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, m_alpha_00, m_rho_00) + diag_matrix(rep_vector(1e-10, N_predict)));
	m_predict_00 = m_base_00 + L_cov * m_f_tilde_00;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, m_alpha_10, m_rho_10) + diag_matrix(rep_vector(1e-10, N_predict)));
	m_predict_10 = m_base_10 + L_cov * m_f_tilde_10;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, m_alpha_01, m_rho_01) + diag_matrix(rep_vector(1e-10, N_predict)));
	m_predict_01 = m_base_01 + L_cov * m_f_tilde_01;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, m_alpha_11, m_rho_11) + diag_matrix(rep_vector(1e-10, N_predict)));
	m_predict_11 = m_base_11 + L_cov * m_f_tilde_11;
	// GP on log sigma of lognormal
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, s_alpha_00, s_rho_00) + diag_matrix(rep_vector(1e-10, N_predict)));
	s_predict_00 = exp( s_base_00 + L_cov * s_f_tilde_00);
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, s_alpha_10, s_rho_10) + diag_matrix(rep_vector(1e-10, N_predict)));
	s_predict_10 = exp( s_base_10 + L_cov * s_f_tilde_10);
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, s_alpha_01, s_rho_01) + diag_matrix(rep_vector(1e-10, N_predict)));
	s_predict_01 = exp( s_base_01 + L_cov * s_f_tilde_01);
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, s_alpha_11, s_rho_11) + diag_matrix(rep_vector(1e-10, N_predict)));
	s_predict_11 = exp( s_base_11 + L_cov * s_f_tilde_11);
}
			
model {
	m_rho_00 ~ normal(0, m_rho_hyper_par_00);
	m_alpha_00 ~ normal(0, m_alpha_hyper_par_00);	
	m_base_00 ~ normal( 0 , 10 );
	s_rho_00 ~ normal(0, s_rho_hyper_par_00);
	s_alpha_00 ~ normal(0, s_alpha_hyper_par_00);	
	s_base_00 ~ normal( 0 , 10 );	
	m_f_tilde_00 ~ normal(0, 1);
	s_f_tilde_00 ~ normal(0, 1);

	m_rho_10 ~ normal(0, m_rho_hyper_par_10);
	m_alpha_10 ~ normal(0, m_alpha_hyper_par_10);	
	m_base_10 ~ normal( 0 , 10 );
	s_rho_10 ~ normal(0, s_rho_hyper_par_10);
	s_alpha_10 ~ normal(0, s_alpha_hyper_par_10);	
	s_base_10 ~ normal( 0 , 10 );	
	m_f_tilde_10 ~ normal(0, 1);
	s_f_tilde_10 ~ normal(0, 1);

	m_rho_01 ~ normal(0, m_rho_hyper_par_01);
	m_alpha_01 ~ normal(0, m_alpha_hyper_par_01);	
	m_base_01 ~ normal( 0 , 10 );
	s_rho_01 ~ normal(0, s_rho_hyper_par_01);
	s_alpha_01 ~ normal(0, s_alpha_hyper_par_01);	
	s_base_01 ~ normal( 0 , 10 );	
	m_f_tilde_01 ~ normal(0, 1);
	s_f_tilde_01 ~ normal(0, 1);

	m_rho_11 ~ normal(0, m_rho_hyper_par_11);
	m_alpha_11 ~ normal(0, m_alpha_hyper_par_11);	
	m_base_11 ~ normal( 0 , 10 );
	s_rho_11 ~ normal(0, s_rho_hyper_par_11);
	s_alpha_11 ~ normal(0, s_alpha_hyper_par_11);	
	s_base_11 ~ normal( 0 , 10 );	
	m_f_tilde_11 ~ normal(0, 1);
	s_f_tilde_11 ~ normal(0, 1);

	target+= normal_lpdf( y_observed_00 | m_predict_00[observed_idx_00], s_predict_00[observed_idx_00] .* inv_sqrt_n_observed_00);
	target+= inv_chi_square_lpdf( s_predict_00[observed_idx_00] .* s_predict_00[observed_idx_00] .* inv_scaled_s2_observed_00 | n_m1_observed_00 );
	target+= normal_lpdf( y_observed_10 | m_predict_10[observed_idx_10], s_predict_10[observed_idx_10] .* inv_sqrt_n_observed_10);
	target+= inv_chi_square_lpdf( s_predict_10[observed_idx_10] .* s_predict_10[observed_idx_10] .* inv_scaled_s2_observed_10 | n_m1_observed_10 );
	target+= normal_lpdf( y_observed_01 | m_predict_01[observed_idx_01], s_predict_01[observed_idx_01] .* inv_sqrt_n_observed_01);
	target+= inv_chi_square_lpdf( s_predict_01[observed_idx_01] .* s_predict_01[observed_idx_01] .* inv_scaled_s2_observed_01 | n_m1_observed_01 );
	target+= normal_lpdf( y_observed_11 | m_predict_11[observed_idx_11], s_predict_11[observed_idx_11] .* inv_sqrt_n_observed_11);
	target+= inv_chi_square_lpdf( s_predict_11[observed_idx_11] .* s_predict_11[observed_idx_11] .* inv_scaled_s2_observed_11 | n_m1_observed_11 );  	
}

generated quantities {
	vector[N_predict] exp_m_predict_00;
	vector[N_predict] exp_m_predict_10;
	vector[N_predict] exp_m_predict_01;
	vector[N_predict] exp_m_predict_11; 
	exp_m_predict_00 = exp( m_predict_00 + s_predict_00 .* s_predict_00 .* rep_vector(0.5, N_predict) );
	exp_m_predict_10 = exp( m_predict_10 + s_predict_10 .* s_predict_10 .* rep_vector(0.5, N_predict) );
	exp_m_predict_01 = exp( m_predict_01 + s_predict_01 .* s_predict_01 .* rep_vector(0.5, N_predict) );
	exp_m_predict_11 = exp( m_predict_11 + s_predict_11 .* s_predict_11 .* rep_vector(0.5, N_predict) );  
}	
"
	
	stan.model <- stan_model(model_name= 'gp_all_logscale',model_code = gsub('\t',' ',stan.code4))	
	stan.data <- list()	
	stan.data$x_predict <- seq(vla[, min(AGE_LABEL)], vla[, max(AGE_LABEL)+1], 0.5)
	stan.data$N_predict <- length(stan.data$x_predict)
	stan.data$observed_idx_00 <- which(stan.data$x_predict%%1==0.5)
	stan.data$observed_idx_00 <- stan.data$observed_idx_00[ subset(vla,SEX==0 & LOC==0)[, !is.na(VLNZ_S2)] ]
	stan.data$N_observed_00 <- length(stan.data$observed_idx_00)	
	stan.data$observed_idx_10 <- which(stan.data$x_predict%%1==0.5)
	stan.data$observed_idx_10 <- stan.data$observed_idx_10[ subset(vla,SEX==1 & LOC==0)[, !is.na(VLNZ_S2)] ]
	stan.data$N_observed_10 <- length(stan.data$observed_idx_10)	
	stan.data$observed_idx_01 <- which(stan.data$x_predict%%1==0.5)
	stan.data$observed_idx_01 <- stan.data$observed_idx_01[ subset(vla,SEX==0 & LOC==1)[, !is.na(VLNZ_S2)] ]
	stan.data$N_observed_01 <- length(stan.data$observed_idx_01)	
	stan.data$observed_idx_11 <- which(stan.data$x_predict%%1==0.5)
	stan.data$observed_idx_11 <- stan.data$observed_idx_11[ subset(vla,SEX==1 & LOC==1)[, !is.na(VLNZ_S2)] ]
	stan.data$N_observed_11 <- length(stan.data$observed_idx_11)
	stan.data$y_observed_00 <- vla[SEX==0 & LOC==0 & !is.na(VLNZ_S2),VLNZ_MEAN]
	stan.data$y_observed_10 <- vla[SEX==1 & LOC==0 & !is.na(VLNZ_S2),VLNZ_MEAN]
	stan.data$y_observed_01 <- vla[SEX==0 & LOC==1 & !is.na(VLNZ_S2),VLNZ_MEAN]
	stan.data$y_observed_11 <- vla[SEX==1 & LOC==1 & !is.na(VLNZ_S2),VLNZ_MEAN]
	stopifnot( length(stan.data$y_observed_00)==stan.data$N_observed_00 )
	stopifnot( length(stan.data$y_observed_10)==stan.data$N_observed_10 )
	stopifnot( length(stan.data$y_observed_01)==stan.data$N_observed_01 )
	stopifnot( length(stan.data$y_observed_11)==stan.data$N_observed_11 )
	stan.data$n_observed_00 <- vla[SEX==0 & LOC==0 & !is.na(VLNZ_S2),VLNZ_N]
	stan.data$n_observed_10 <- vla[SEX==1 & LOC==0 & !is.na(VLNZ_S2),VLNZ_N]
	stan.data$n_observed_01 <- vla[SEX==0 & LOC==1 & !is.na(VLNZ_S2),VLNZ_N]
	stan.data$n_observed_11 <- vla[SEX==1 & LOC==1 & !is.na(VLNZ_S2),VLNZ_N]
	stopifnot( length(stan.data$n_observed_00)==stan.data$N_observed_00 )
	stopifnot( length(stan.data$n_observed_10)==stan.data$N_observed_10 )
	stopifnot( length(stan.data$n_observed_01)==stan.data$N_observed_01 )
	stopifnot( length(stan.data$n_observed_11)==stan.data$N_observed_11 )
	stan.data$s2_observed_00 <- vla[SEX==0 & LOC==0 & !is.na(VLNZ_S2),VLNZ_S2]
	stan.data$s2_observed_10 <- vla[SEX==1 & LOC==0 & !is.na(VLNZ_S2),VLNZ_S2]
	stan.data$s2_observed_01 <- vla[SEX==0 & LOC==1 & !is.na(VLNZ_S2),VLNZ_S2]
	stan.data$s2_observed_11 <- vla[SEX==1 & LOC==1 & !is.na(VLNZ_S2),VLNZ_S2]
	stopifnot( length(stan.data$s2_observed_00)==stan.data$N_observed_00 )
	stopifnot( length(stan.data$s2_observed_10)==stan.data$N_observed_10 )
	stopifnot( length(stan.data$s2_observed_01)==stan.data$N_observed_01 )
	stopifnot( length(stan.data$s2_observed_11)==stan.data$N_observed_11 )
	stan.data$m_rho_hyper_par_00 <- diff(range(stan.data$x_predict))/3
	stan.data$m_alpha_hyper_par_00 <- 2
	stan.data$s_rho_hyper_par_00 <- diff(range(stan.data$x_predict))/3
	stan.data$s_alpha_hyper_par_00 <- 2
	stan.data$m_rho_hyper_par_10 <- diff(range(stan.data$x_predict))/3
	stan.data$m_alpha_hyper_par_10 <- 2
	stan.data$s_rho_hyper_par_10 <- diff(range(stan.data$x_predict))/3
	stan.data$s_alpha_hyper_par_10 <- 2
	stan.data$m_rho_hyper_par_01 <- diff(range(stan.data$x_predict))/3
	stan.data$m_alpha_hyper_par_01 <- 2
	stan.data$s_rho_hyper_par_01 <- diff(range(stan.data$x_predict))/3
	stan.data$s_alpha_hyper_par_01 <- 2
	stan.data$m_rho_hyper_par_11 <- diff(range(stan.data$x_predict))/3
	stan.data$m_alpha_hyper_par_11 <- 2
	stan.data$s_rho_hyper_par_11 <- diff(range(stan.data$x_predict))/3
	stan.data$s_alpha_hyper_par_11 <- 2	
	fit <- sampling(stan.model, data=stan.data, iter=20e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999) )
	save(fit, file=file.path(outdir, "200428f_hivmeans_gp_stanfit.rda"))
	
	stan.code2b <- "
data{	
	int<lower=1> N_predict;
  	real x_predict[N_predict];
  	int<lower=1> N_observed;
  	int<lower=1, upper=N_predict> observed_idx[N_observed];
  	int y_observed_00[N_observed];
	int y_observed_10[N_observed];
	int y_observed_01[N_observed];	
	int y_observed_11[N_observed];	
	int total_observed_00[N_observed];
	int total_observed_10[N_observed];
	int total_observed_01[N_observed];
	int total_observed_11[N_observed];
  	real<lower=0> rho_hyper_par_00;
	real<lower=0> rho_hyper_par_10;
	real<lower=0> rho_hyper_par_01;
	real<lower=0> rho_hyper_par_11;  	
  	real<lower=0> alpha_hyper_par_00;
	real<lower=0> alpha_hyper_par_10;
	real<lower=0> alpha_hyper_par_01;
	real<lower=0> alpha_hyper_par_11;
}

parameters {
	real<lower=0> rho_00;
	real<lower=0> rho_10;
	real<lower=0> rho_01;
	real<lower=0> rho_11;
	real<lower=0> alpha_00;
	real<lower=0> alpha_10;
	real<lower=0> alpha_01;
	real<lower=0> alpha_11;
	real sex0_loc0;
	real sex1_loc0;
	real sex0_loc1;
	real sex1_loc1;  	
  	vector[N_predict] f_tilde_00;
	vector[N_predict] f_tilde_10;
	vector[N_predict] f_tilde_01;
	vector[N_predict] f_tilde_11;
}

transformed parameters {
  	matrix[N_predict, N_predict] L_cov;
  	vector[N_predict] logit_p_predict_00;
	vector[N_predict] logit_p_predict_10;
	vector[N_predict] logit_p_predict_01;
	vector[N_predict] logit_p_predict_11;
	// GP for 00 and 01 (women)
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_00, rho_00) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_00 = sex0_loc0 + L_cov * f_tilde_00;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_01, rho_01) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_01 = sex0_loc1 + L_cov * f_tilde_01;
	// GP for 10 and 10 (men)
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_10, rho_10) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_10 = sex1_loc0 + L_cov * f_tilde_10;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_11, rho_11) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_11 = sex1_loc1 + L_cov * f_tilde_11;
}

model {
	rho_00 ~ normal(0, rho_hyper_par_00);
	rho_10 ~ normal(0, rho_hyper_par_10);
	rho_01 ~ normal(0, rho_hyper_par_01);
	rho_11 ~ normal(0, rho_hyper_par_11);  	
  	alpha_00 ~ normal(0, alpha_hyper_par_00);
	alpha_10 ~ normal(0, alpha_hyper_par_10);
	alpha_01 ~ normal(0, alpha_hyper_par_01);
	alpha_11 ~ normal(0, alpha_hyper_par_11);
  	sex0_loc0 ~ normal( 0 , 10 );
	sex0_loc1 ~ normal( 0 , 10 );
	sex1_loc0 ~ normal( 0 , 10 );
	sex1_loc1 ~ normal( 0 , 10 );
  	f_tilde_00 ~ normal(0, 1);
	f_tilde_01 ~ normal(0, 1);
	f_tilde_10 ~ normal(0, 1);
	f_tilde_11 ~ normal(0, 1);
  	y_observed_00 ~ binomial_logit(total_observed_00, logit_p_predict_00[observed_idx] );
	y_observed_01 ~ binomial_logit(total_observed_01, logit_p_predict_01[observed_idx] );
	y_observed_10 ~ binomial_logit(total_observed_10, logit_p_predict_10[observed_idx] );
	y_observed_11 ~ binomial_logit(total_observed_11, logit_p_predict_11[observed_idx] );
}

generated quantities {
  	vector[N_predict] p_predict_00;
	vector[N_predict] p_predict_01;
	vector[N_predict] p_predict_10;
	vector[N_predict] p_predict_11;
	p_predict_00 = inv_logit(logit_p_predict_00);
	p_predict_01 = inv_logit(logit_p_predict_01);  
	p_predict_10 = inv_logit(logit_p_predict_10);
	p_predict_11 = inv_logit(logit_p_predict_11);
}
"
		
	stan.modelB <- stan_model(model_name= 'gp_all',model_code = gsub('\t',' ',stan.code2b))
	stan.dataB <- list()
	stan.dataB$x_predict <- seq(vla[, min(AGE_LABEL)], vla[, max(AGE_LABEL)+1], 0.5)
	stan.dataB$N_predict <- length(stan.dataB$x_predict)
	stan.dataB$observed_idx <- which(stan.dataB$x_predict%%1==0.5)
	stan.dataB$N_observed <- length(stan.dataB$observed_idx)
	stan.dataB$y_observed_00 <- vla[SEX==0 & LOC==0, N-VLNZ_N]
	stan.dataB$y_observed_10 <- vla[SEX==1 & LOC==0, N-VLNZ_N]
	stan.dataB$y_observed_01 <- vla[SEX==0 & LOC==1, N-VLNZ_N]
	stan.dataB$y_observed_11 <- vla[SEX==1 & LOC==1, N-VLNZ_N]
	stan.dataB$total_observed_00 <- vla[SEX==0 & LOC==0, N]
	stan.dataB$total_observed_10 <- vla[SEX==1 & LOC==0, N]
	stan.dataB$total_observed_01 <- vla[SEX==0 & LOC==1, N]
	stan.dataB$total_observed_11 <- vla[SEX==1 & LOC==1, N]
	stan.dataB$rho_hyper_par_00 <- diff(range(stan.dataB$x_predict))/3
	stan.dataB$rho_hyper_par_10 <- diff(range(stan.dataB$x_predict))/3
	stan.dataB$rho_hyper_par_01 <- diff(range(stan.dataB$x_predict))/3
	stan.dataB$rho_hyper_par_11 <- diff(range(stan.dataB$x_predict))/3
	stan.dataB$alpha_hyper_par_00 <- 2
	stan.dataB$alpha_hyper_par_10 <- 2
	stan.dataB$alpha_hyper_par_01 <- 2
	stan.dataB$alpha_hyper_par_11 <- 2
	fitB <- sampling(stan.modelB, data=stan.dataB, iter=2e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999))
	save(fitB, file=file.path(outdir, "200428f_hivzeros_gp_stanfit.rda"))
	
	
	load( file.path(outdir, "200428f_hivzeros_gp_stanfit.rda") )
	load( file.path(outdir, "200428f_hivmeans_gp_stanfit_10.rda") )
	fit10 <- fit
	load( file.path(outdir, "200428f_hivmeans_gp_stanfit_11.rda") )
	fit11 <- fit
	load( file.path(outdir, "200428f_hivmeans_gp_stanfit_01.rda") )
	fit01 <- fit
	load( file.path(outdir, "200428f_hivmeans_gp_stanfit.rda") )
	
	#
	#	what we expect roughly
	ggplot(vla2, aes(x=AGE_LABEL)) +
			geom_point(aes(y= VLNZ_N/N * exp(VLNZ_MEAN+VLNZ_S2/2)))
	
	tmp <- summary(fit)$summary
	tmp[grepl('exp_m_predict',rownames(tmp)),]
	tmp[grepl('m_predict',rownames(tmp)),]
	tmp[grepl('s_predict',rownames(tmp)),]
	
	#
	#	get posterior samples
	reB <- rstan::extract(fitB)
	re01 <- rstan::extract(fit01)
	re10 <- rstan::extract(fit10)
	re11 <- rstan::extract(fit11)
	re <- rstan::extract(fit)
	re$exp_m_predict_10 <- re10$exp_m_predict
	re$exp_m_predict_01 <- re01$exp_m_predict
	re$exp_m_predict_11 <- re11$exp_m_predict
	re$m_rho_10 <- re10$m_rho
	re$m_rho_01 <- re01$m_rho	
	re$m_rho_11 <- re11$m_rho
	re$s_rho_10 <- re10$s_rho
	re$s_rho_01 <- re01$s_rho
	re$s_rho_11 <- re11$s_rho	
	
	
	#
	#	extract hyperparams rho
	ps <- c(0.025,0.25,0.5,0.75,0.975)
	tmp <- cbind( quantile(re$m_rho_00, probs=ps),
			quantile(re$m_rho_10, probs=ps),
			quantile(re$m_rho_01, probs=ps),
			quantile(re$m_rho_11, probs=ps),
			quantile(re$m_alpha_00, probs=ps),
			quantile(re$m_alpha_10, probs=ps),
			quantile(re$m_alpha_01, probs=ps),
			quantile(re$m_alpha_11, probs=ps), 	
			quantile(re$s_rho_00, probs=ps),
			quantile(re$s_rho_10, probs=ps),
			quantile(re$s_rho_01, probs=ps),
			quantile(re$s_rho_11, probs=ps),
			quantile(re$s_alpha_00, probs=ps),
			quantile(re$s_alpha_10, probs=ps),
			quantile(re$s_alpha_01, probs=ps),
			quantile(re$s_alpha_11, probs=ps)			
			)			
	colnames(tmp) <- c(	'm_rho_00','m_rho_10','m_rho_01','m_rho_11','m_alpha_00','m_alpha_10','m_alpha_01','m_alpha_11',
						's_rho_00','s_rho_10','s_rho_01','s_rho_11','s_alpha_00','s_alpha_10','s_alpha_01','s_alpha_11')
	rownames(tmp) <- c('CL','IL','M','IU','CU')
	tmp <- as.data.table(reshape2::melt(tmp))
	setnames(tmp, 'Var2', 'GP_hyper_par')
	tmp[, SEX:= as.integer(gsub('^([a-z]+_[a-z]+)_([0-9])([0-9])','\\2',GP_hyper_par))]
	tmp[, LOC:= as.integer(gsub('^([a-z]+_[a-z]+)_([0-9])([0-9])','\\3',GP_hyper_par))]
	tmp[, GP_hyper_par:= gsub('^([a-z]+_[a-z]+)_([0-9])([0-9])','\\1',GP_hyper_par)]
	tmp <- dcast.data.table(tmp, LOC+SEX+GP_hyper_par~Var1, value.var='value')
	mvl.gp.pars <- merge(unique(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), tmp, by=c('SEX','LOC'))
	ggplot(mvl.gp.pars, aes(x= paste0(GP_hyper_par, ' ', LOC_LABEL, ' ', SEX_LABEL))) +
			geom_point(aes(y=M)) +
			geom_errorbar(aes(ymin=CL, ymax=CU)) +
			coord_flip() +
			theme_bw() +
			labs(x='GP hyperparameter\n', y='')
	ggsave(file=file.path(prjdir,'results_200220','200428f_notsuppAmongInfected_gppars.pdf'), w=6, h=3)
	
	#
	#	extract MVL by gender and location
	ps <- c(0.025,0.5,0.975)
	tmp2 <- sample(seq_len(nrow(reB$p_predict_00)),nrow(re$exp_m_predict_00),replace=TRUE)	
	tmp <- cbind( (1-reB$p_predict_00[tmp2, ]) * re$exp_m_predict_00,
			(1-reB$p_predict_10[tmp2, ]) * re$exp_m_predict_10,
			(1-reB$p_predict_01[tmp2, ]) * re$exp_m_predict_01,
			(1-reB$p_predict_11[tmp2, ]) * re$exp_m_predict_11	)		
	rp <- as.data.table(reshape2::melt(tmp))
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
	rp[, LABEL:= paste0(sprintf('%2.0f',M),' (',sprintf('%2.0f',CL),'-',sprintf('%2.0f',CU),')') ]
	mvl.by.sex.loc <- copy(rp)
	#	   LOC_LABEL SEX_LABEL        CL         CU         M             LABEL
	#1:   fishing         M 6459.2055 11298.7826 8318.0350 	8318 (6459-11299)
	#2:   fishing         F 2191.5886  4269.4039 2968.2880  2968 (2192-4269)
	#3:    inland         M 1054.5387  1950.9344 1405.3776  1405 (1055-1951)
	#4:    inland         F  559.6118   887.5624  698.1173     698 (560-888)	


	#	extract MVL ratio female:male and male:female
	tmp2 <- sample(seq_len(nrow(reB$p_predict_00)),nrow(re$exp_m_predict_00),replace=TRUE)	
	tmp <- cbind( (1-reB$p_predict_00[tmp2, ]) * re$exp_m_predict_00,
			(1-reB$p_predict_10[tmp2, ]) * re$exp_m_predict_10,
			(1-reB$p_predict_01[tmp2, ]) * re$exp_m_predict_01,
			(1-reB$p_predict_11[tmp2, ]) * re$exp_m_predict_11	)		
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
	rp[, LABEL:= paste0(sprintf('%2.2f',M),' (',sprintf('%2.2f',CL),' - ',sprintf('%2.2f',CU),')') ]	
	mvl.ratio.by.loc <- copy(rp)	
	#	   LOC_LABEL variable        CL        CU         M              LABEL
	#1:   fishing    PR_FM 0.2307099 0.5531012 0.3560961 0.36 (0.23 - 0.55)
	#2:   fishing    PR_MF 1.8079873 4.3344482 2.8082310 2.81 (1.81 - 4.33)
	#3:    inland    PR_FM 0.3349882 0.7226552 0.4958447 0.50 (0.33 - 0.72)
	#4:    inland    PR_MF 1.3837858 2.9851794 2.0167604 2.02 (1.38 - 2.99)


	#
	#	extract MVL by age for gender and location
	ps <- c(0.025,0.5,0.975)
	tmp2 <- sample(seq_len(nrow(reB$p_predict_00)),nrow(re$exp_m_predict_00),replace=TRUE)	
	tmp <- cbind( (1-reB$p_predict_00[tmp2, ]) * re$exp_m_predict_00,
			(1-reB$p_predict_10[tmp2, ]) * re$exp_m_predict_10,
			(1-reB$p_predict_01[tmp2, ]) * re$exp_m_predict_01,
			(1-reB$p_predict_11[tmp2, ]) * re$exp_m_predict_11	)	
	#tmp <- cbind( re$exp_m_predict_00, re$exp_m_predict_10, re$exp_m_predict_01, re$exp_m_predict_11	)	
	tmp <- apply(tmp, 2, quantile, probs=ps)
	rownames(tmp) <- c('CL','M','CU')
	tmp <- as.data.table(reshape2::melt(tmp))	
	tmp <- dcast.data.table(tmp, Var2~Var1, value.var='value')
	mvl.by.age <- cbind( as.data.table(expand.grid(AGE_LABEL= stan.dataB$x_predict, SEX=c(0,1), LOC=c(0,1))), tmp )
	mvl.by.age <- merge(unique(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), mvl.by.age, by=c('LOC','SEX'))
	ggplot(mvl.by.age) + 		
			geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2) +
			geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL)) +
			scale_x_continuous( expand=c(0,0) ) + 
			#scale_y_log10() +
			scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
			facet_wrap(~LOC_LABEL, ncol=2, scales='free_y') +
			theme_bw() +
			labs(x='\nage at visit (years)', 
					y='mean viral load\n(95% credibility interval)\n', 
					colour='gender', 
					linetype='location')
	ggsave(file=file.path(prjdir,'results_200220','200428f_mvl_vs_age_by_gender_fishinland_stan.pdf'), w=7, h=5)

	
	#
	#	extract MVL ratio by age for gender and location
	tmp2 <- sample(seq_len(nrow(reB$p_predict_00)),nrow(re$exp_m_predict_00),replace=TRUE)	
	tmp <- cbind( (1-reB$p_predict_00[tmp2, ]) * re$exp_m_predict_00,
			(1-reB$p_predict_10[tmp2, ]) * re$exp_m_predict_10,
			(1-reB$p_predict_01[tmp2, ]) * re$exp_m_predict_01,
			(1-reB$p_predict_11[tmp2, ]) * re$exp_m_predict_11	)	
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
	rp[, LABEL:= paste0(sprintf('%2.2f',M),' (',sprintf('%2.2f',CL),' - ',sprintf('%2.2f',CU),')') ]	
	mvl.ratio.by.loc.age <- copy(rp)
	
	
		
	save(vla, reB, re, re11, re01, re10, mvl.gp.pars, mvl.by.sex.loc, mvl.ratio.by.loc, mvl.by.age, mvl.ratio.by.loc.age, file=file.path(outdir, "200428f_mvl.rda"))
	
	
	#	make table version suppressed
	mvl.by.age[, LABEL:= paste0(sprintf('%2.0f',M),' (',sprintf('%2.0f',CL),'-',sprintf('%2.0f',CU),')') ]
	set(mvl.by.age, NULL, 'SEX_LABEL', mvl.by.age[, factor(as.character(SEX_LABEL), levels=c('F','M'))])	
	setnames(mvl.ratio.by.loc.age,'LABEL', 'LABEL2')
	dt <- subset(mvl.by.age, AGE_LABEL%in%c(20.5,25.5,30.5,35.5,40.5,45.5))	
	dt <- dcast.data.table(dt, LOC_LABEL+AGE_LABEL~SEX_LABEL, value.var='LABEL')
	tmp <- subset(mvl.ratio.by.loc.age, variable=="PR_FM" & AGE_LABEL%in%c(20.5,25.5,30.5,35.5,40.5,45.5), c(LOC_LABEL, AGE_LABEL, LABEL2))
	dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
	tmp <- subset(mvl.ratio.by.loc.age, variable=="PR_MF" & AGE_LABEL%in%c(20.5,25.5,30.5,35.5,40.5,45.5), c(LOC_LABEL, AGE_LABEL, LABEL2))
	setnames(tmp, 'LABEL2', 'PR_MF')
	dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
	write.csv(dt, row.names=FALSE, file=file.path(outdir, "200428f_mvl.csv"))
	
}
