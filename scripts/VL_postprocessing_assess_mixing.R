self_relative_path <- "scripts/VL_postprocessing_assess_mixing.R"

########################
cat("\nStart of:", self_relative_path, "\n")
########################

#################
#      AIMS     #
#################

# Load a STAN model and assess diagnostics + mixing of chains

# --- Merge fit to diagnostics ---
# Error in dim(x) <- length(x) : 
#   invalid first argument, must be vector (list or atomic) 
# Calls: plot.group.mcmc.parcoord ... mcmc_parcoord_data -> prepare_mcmc_array -> as.array -> as.array.default
# 

#################
#   Libraries   #
#################

# load libraries
{
    library(data.table)
    library(ggplot2)
    library(bayesplot)
    library(posterior)
    if(system.file(package='rstan') != ""){
        library(rstan)
    } else {
        library(cmdstanr)
    }
}

# paths
self_relative_path <- "scripts/VL_postprocessing_assess_mixing.R"
if (interactive()) {
    here::i_am(self_relative_path)
    gitdir <- here::here()
} else {
    cmd <- commandArgs()
    cmd <- cmd[cmd %like% "file"]
    gitdir <- gsub(paste0("--file=(.*)/", self_relative_path, "$"), "\\1", cmd)
}

# helpers
source(file.path(gitdir, "R/paths.R"))

# options:
args <- args[names(args) %like% "jobname|run|out.dir|round"]
if(interactive()){
    args$jobname <- 'cmdstan_alpha100sharedhyper_vl_1000'
    args$run.gp.prevl <- TRUE
    args$round <- 19
}
print(args)

##################
#      Main      #
##################

catn("Locate model fit and output directories")
# _____________________________________________

outdir <- args$out.dir.exact
if(is.na(outdir)){
    outdir <- with(args, {
        model <- which( c( 
            `run-gp-prevl` = run.gp.prevl,
            `run-gp-supp-hiv` = run.gp.supp.hiv,
            `run-gp-supp-pop` = run.gp.supp.pop
        )) |> names()
        stopifnot("One of the models needs be specified" = length(model) > 0)

        if(length(model) > 1){
            warning("Only assessing mixing for the model:", model[1])
        }

        file.path(out.dir.prefix, jobname, model[1])
    })
    stopifnot("Output directory does not exist" = dir.exists(outdir))
}

path.stan.output <- list.files( 
    outdir, 
    paste0('round', args$round, '.rds'), 
    full.names =TRUE)
# outfile.prefix <- gsub( '.rds$', '-', path.stan.output)
outfile.diagnostics <- file.path(outdir, paste0('diagnostics/assess_round', args$round)) 
outfile.figures <- file.path(outdir, paste0('figures/assess_round', args$round))
sapply(c(outfile.figures, outfile.diagnostics), function(x){
    d <- dirname(x)
    if(! dir.exists(d)) dir.create(d)
})


catn("Loading model fit and samples")
# __________________________________

fit <- readRDS(path.stan.output)
fit.type <- get.fit.type(fit)
samples <- posterior::as_draws_df(fit)
n_iter <- nrow(samples)

catn("Extracting samples")
# ________________________

# note samples are only really needed for WAIC and LOO here...

diagns <- make_convergence_diagnostics_stats(
    fit=fit,
    re=samples,
    exclude_rgx = 'L_cov|rho_hyper_par|_rep',
    outfile.diagnostics)

catn("Merge fit to diagnostics")
# ________________________________

p <- plot.group.mcmc.parcoord(fit=fit, re=samples,'logit_p_predict')
ggsave(p, file=paste0(outfile.figures, '-mcmc-parcoord-p_predict.png'), w=8, h=8)

catn("Traceplots")
# ________________

regex_gppars <- c('rho_[0-1]', 'alpha_[0-1]')
p <- bayesplot::mcmc_trace(samples, regex_pars = regex_gppars, facet_args = list(ncol=4)) + 
    theme_default() 
ggsave(p, file = paste0(outfile.figures, '-mcmc-trace_plots.png'), w  = 8, h = 8)


catn("Interval plot")
# ___________________

# baseline parameters
inverse_logit <- function(x){
    return(1/(1+exp(-x)))
}

p <- bayesplot::mcmc_intervals(
    samples, 
    transformations = inverse_logit,
    regex_pars = c('sex[0-1]_loc[0-1]')
) + 
    theme_default() + 
    scale_y_discrete(labels = dict_stan_params_t) +
    labs(title="Baseline parameters (under inverse logit transform)") 
ggsave(p, file = paste0(outfile.figures, '-mcmc-intervals_baseline.png'), w  = 8, h = 8)

# hyperparameters
p <- bayesplot::mcmc_intervals(
    samples, 
    regex_pars = c('rho_[0-1]', 'alpha_[0-1]')
) + 
    scale_y_discrete(labels = dict_stan_params2) + 
    theme_default()
ggsave(p, file = paste0(outfile.figures, '-mcmc-intervals_hyper.png'), w  = 8, h = 8)


catn("Pairs plot") 
# _________________

# hyperparameters
p <- lapply(c('00', '01', '10', '11'), function(group){
    p <- bayesplot::mcmc_pairs(
        samples,
        regex_pars = paste0(c('rho_', 'alpha_'), group), 
        facet_args = list(ncol=4),
        np=bayesplot::nuts_params(fit)
    ) }
)  |> ggpubr::ggarrange(plotlist=_, ncol=2, nrow=2)
ggsave(p, file = paste0(outfile.figures, '-mcmc-pairs_plot_gphyperparams.png'), w  = 10, h = 10)

#####################
catn("End of script")
#####################
