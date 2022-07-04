# Prepare some fake Data
library(data.table)
library(cmdstanr)

stan.file <- 'simple_multinomial_logit_260522.stan'
mod <- cmdstan_model(stan.file)

# helpers
softmax <- function(x)
{
        x <- x - max(x)
        exp(x)/sum(exp(x))
}

# Basics
N <- 10000
K <- 2

stan_data <- list(
        N = N,
        K = K
)

# Simulate scenario:
dsim <- data.table(
                sex = sample(c('M', 'F'), N, replace=TRUE),
                comm = sample(c('fishing', 'inland', 'other'), N, replace=TRUE, prob=c(.4, .3, .1))
)

dsim[, p_base := dplyr::if_else(sex=='M', .6, .5)]
dsim[, DUMMY := 1:.N]
dsim[, deltaComm := 0] 
dsim[sex=='M', deltaComm := switch(comm, 
                                fishing=.15,
                                inland=0,
                                other=.1
        ) , by=DUMMY] 
dsim[sex=='F', deltaComm := switch(comm, 
                                fishing=.1,
                                inland=0,
                                other=-.1
        ) , by=DUMMY] 
dsim[, DUMMY := NULL]
dsim[, p := p_base + deltaComm]
dsim[, u := runif(.N)]
dsim[, y := as.integer(u < p) ]
dsim

ddata <- subset(dsim, select=c('sex', 'comm', 'y'))
cols <- c('sex', 'comm')
ddata[, (cols) := lapply(.SD, as.factor), .SDcols=cols]

mf <- model.frame(y ~ . , data=ddata)
stan_data$y <- mf$y + 1
stan_data$x <- model.matrix(y~., data=mf) 
stan_data$D <- dim(stan_data$x)[2]


# Fit stan model
# __

# rstan_options(auto_write=TRUE)
options(mc.cores = parallel::detectCores())

fit <- mod$sample(data=stan_data, seed=1)
fit$summary()
fit$sampler_diagnostics()
params <- fit$summary()$variable

# Ok posterior medians make sense!
# __
nms <- grep('^beta\\[', params, value=TRUE)
idx <- params %in% nms

pmeds <- fit$summary()[idx, c('variable', 'median')]
pmeds <- matrix(pmeds$median, nrow=stan_data$D)

tmp <- unique(stan_data$x)
tmp <- as.data.table(tmp)
setkeyv(tmp, names(tmp))
tmp <- as.matrix(tmp)

pmeds <- tmp %*% pmeds 
pmeds <- t(apply(pmeds, 1, softmax))
colnames(pmeds) <- paste0('P(K=', 1:stan_data$K, ')')
pmeds <- cbind(tmp, pmeds)
pmeds <- as.data.table(pmeds)

if( any(grepl('sex', names(pmeds))) )
{
        pmeds[, sex := dplyr::if_else(sexM == 1, 'M', 'F')]
}

if( any(grepl('comm', names(pmeds))) )
{
        cols <- grep('comm', names(pmeds), value=TRUE)
        types <- c('inland', 'fishing', 'other')
        .f <- function(x) gsub('comm', '', x)
        idx <- types %in% .f(cols)
        baseline <- types[which(!idx)]

        .check <- function(x, N) (x %% N == 0 & x > 0)
        pmeds[, comm_DUMMY := .SD[[1]]*2 + .SD[[2]]*3, .SDcols=cols]
        pmeds[, comm := dplyr::if_else(comm_DUMMY == 0, baseline, '')]
        pmeds[.check(comm_DUMMY, 2), comm :=  .f(cols[2-1])]
        pmeds[.check(comm_DUMMY, 3), comm :=  .f(cols[3-1])]
        pmeds
}
# It may be appropriate not to call pmeds as such: 
# the cool features are not the medians (which will be thrown away)
# Instead, they are translation between factors and contrast representations

dsim
merge(
      pmeds[, .SD, .SDcols=c('sex','comm', 'P(K=1)')],
      dsim[, .(p=sum(y == 0)/.N, .N), by=c('sex', 'comm')]
)
# we seem to recover nicely but there seems to be some shrinkage towards 0.5


# Extract posterior draws
# __
pdraws <- fit$draws(format='array', inc_warmup=FALSE)
params <- dimnames(pdraws)[["variable"]]

select.chains <- seq_along(dimnames(pdraws)[['chain']])
iters <- 1:(fit$metadata()[['iter_sampling']])
tot.iters <- prod(dim(pdraws)[c(1,2)])

pos <- list()
{
	# bring samples into long array format (change for beta, cause it s a matrix)
        tmp <- pdraws[, select.chains, grepl('beta\\[', params)]
        tmp <- apply(tmp, 3, rbind)
        idx <- c(dim(tmp)[1], stan_data$D, stan_data$K)
        tmp <- array(tmp, dim=idx)
        pos$betas <- tmp
}

#

mat_group = unique(stan_data$x)
mat_group
pa <- lapply(1:tot.iters, function(idx)
       {
               tmp <- mat_group %*% pos$betas[idx,,] 
               tmp <- t(apply(tmp, MARGIN=1, softmax))
               tmp
       }
)
pa <- simplify2array(pa) # dim(pa)

apply(pa, 1:2, mean)
apply(pa, 1:2, sd)
apply(pa, 1:2, median)
apply(pa, 1:2, quantile, .75)




# plot
# https://mc-stan.org/bayesplot/articles/visual-mcmc-diagnostics.html
library(bayesplot)
mcmc_hist(fit$draws())
mcmc_parcoord(fit$draws())
mcmc_pairs(fit$draws())
color_scheme_set("mix-brightblue-gray")
mcmc_trace(fit$draws())

stan_data$x
