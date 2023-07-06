#################
#      AIMS     #
#################

# Load a STAN model and assess diagnostics + mixing of chains

#################
#   Libraries   #
#################

library(data.table)
library(ggplot2)
library(rstan)
library(rstanarm)
library(bayesplot)

# paths
self_relative_path <- "src/postprocessing_assess_mixing.R"
if (interactive()) {
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

with(args, 
    if( run.gp.supp.hiv ){
        1
    }
)

##################
#      Main      #
##################

catn("Locate model fit and output directories")
# _____________________________________________

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
path.stan.output <- list.files( 
    outdir, 
    paste0('round', args$round, '.rds'), 
    full.names =TRUE)
outfile.prefix <- gsub( '.rds$', '-', path.stan.output)
outfile.figures <- file.path(outdir, 'figures')
if(!dir.exists(outfile.figures)) dir.create(outfile.figures)


catn("Loading model fit and samples")
# __________________________________

fit <- readRDS(path.stan.output)

if(! "stanfit" %in% class(fit)){
    fit <- rstan::read_stan_csv(fit$output_files())
}

samples <- rstan::extract(fit, permuted = FALSE, inc_warmup = FALSE)

catn("Extracting samples")
# ________________________

make_convergence_diagnostics_stats(
    fit=fit,
    re=samples,
    exclude_rgx = 'L_cov|rho_hyper_par',
    outfile.prefix)

catn("Merge fit to diagnostics")
# ________________________________

p <- plot.group.mcmc.parcoord('p_predict')
ggsave(p, file=file.path(outfile.figures, '-mcmc-parcoord-p_predict.png'), w=8, h=8)

catn("Traceplots")
# ________________

p <- bayesplot::mcmc_trace(fit, regex_pars = c('rho_', 'alpha_')) + 
    theme_default() +
ggsave(p, file = file.path(outfile.figures, '-mcmc-trace_plots.png'), w  = 8, h = 8)


catn("Interval plot")
# ___________________

# baseline parameters
inverse_logit <- function(x){
    return(1/(1+exp(-x)))
}

p <- bayesplot::mcmc_intervals(
    fit, 
    transformations = inverse_logit,
    regex_pars = c('sex[0-1]_loc[0-1]')
) + 
    theme_default() + 
    scale_y_discrete(labels = dict_stan_params_t) +
    labs(title="Baseline parameters (under inverse logit transform)") 
ggsave(p, file = paste0(outfile.figures, '-mcmc-intervals_baseline'), w  = 8, h = 8)

# hyperparameters
p <- bayesplot::mcmc_intervals(
    fit, 
    regex_pars = c('rho_[0-1]', 'alpha_[0-1]')
) + 
    scale_y_discrete(labels = dict_stan_params2) + 
    theme_default()
ggsave(p, file = paste0(outfile.figures, '-mcmc-intervals_hyper.png'), w  = 8, h = 8)



#
# Pairs plot
#

# hyperparameters
p <- lapply(c('00', '01', '10', '11'), function(group){
    bayesplot::mcmc_pairs(fit, regex_pars = paste0(c('rho_', 'alpha_'), group)) + theme_bw()
}) |> ggpubr::ggarrange(plotlist=_, ncol=2, nrow=2)
ggsave(p, file = paste0(outfile.figures, '-mcmc-pairs_plot_gphyperparams.png'), w  = 10, h = 10)
