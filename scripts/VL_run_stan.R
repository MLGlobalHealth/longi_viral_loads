# AIMS:
# - adapt older code from Oli to multi-round, longitudinal settings
# TODO: discuss: we are removing individuals with missing VLs: they are very little

library(data.table) |> suppressPackageStartupMessages()
library(ggplot2)    |> suppressPackageStartupMessages()
library(Hmisc)      |> suppressPackageStartupMessages()
library(rstan)      |> suppressPackageStartupMessages()
library(optparse)   |> suppressPackageStartupMessages()
library(here)       |> suppressPackageStartupMessages()

################
#    PATHS     #
################

# automatically finding the github directory may be complicated 
# if script is called outside from it.
self_relative_path <- 'scripts/VL_run_stan.R'
parallelise <- FALSE

if( interactive() )
{
    here::i_am(self_relative_path)
    gitdir <- here::here()
} else {
    cmd <- commandArgs()
    cmd <- cmd[cmd %like% 'file']
    gitdir <- gsub(paste0("--file=(.*)/", self_relative_path, '$'), "\\1", cmd)
}

# helpers
source(file.path(gitdir, "R/paths.R"))
source(file.path(gitdir.functions, "phsc_vl_helpers.R"))

# options (automatically sourced in R/options.R)
args <- args[names(args) %like% '^run|viral.load|jobname|indir|out.dir|refit|round|^only.firstparticipants']
if(interactive()){ # testing
    args$only.firstparticipants <- TRUE 
    args$out.dir.prefix <- '/rds/ephemeral/asdllsd.pbs/vl_1000_firstpart/run-gp-supp-pop' 
    args$run.gp.prevl <- TRUE
} 
print(args)

# parallel backend (use multiple cores if local)
if (parallelise) 
    source(file.path(gitdir.R, "local_cores_parallelisation.R"))

file.exists(path.hivstatusvl.r1520) |>
    all() |>
    stopifnot()


################
#     MAIN     #
################

# check study round exists
stopifnot(all(args$round %in% 16:19))

# set viral load thresholds
VL_DETECTABLE <- args$vl.detectable
VIREMIC_VIRAL_LOAD <- args$viremic.viral.load

# specify and create output directories as func(vl, jobname)
if(! is.na(args$out.dir.exact) ){
    vl.out.dir <- args$out.dir.exact
    dir.create(vl.out.dir, recursive = TRUE)
}else{
    stopifnot(dir.exists(args$out.dir.prefix))
    out.dir <- file.path(args$out.dir.prefix)
    suffix <- make.suffix(args)
    vl.out.dir <- file.path( out.dir, suffix)
}
cat('vl.out.dir specified as:\n ', vl.out.dir, '\n')
dir.create(vl.out.dir, showWarnings = FALSE)
stopifnot(dir.exists(vl.out.dir))

# get data
dall <- get.dall(path = path.hivstatusvl.r1520, only_firstpart = args$only.firstparticipants)

if (0) { # Info for introduction to results

    .mean2 <- function(x) paste0(round(100 * mean(x), 2), "%")

    # proportions of viraemic measurements
    dall[VL_COPIES > VIREMIC_VIRAL_LOAD, .mean2(SEX == "F"), ]
    dall[VL_COPIES > VIREMIC_VIRAL_LOAD, .mean2(SEX == "F"), by = "ROUND"]
    dall[VL_COPIES > VIREMIC_VIRAL_LOAD, .mean2(FC == "inland"), by = "ROUND"]

    # ARVs
    dall[HIV_AND_VL == 1, .mean2(is.na(ARVMED)), by = "ROUND"]
    cat("mean log10 VL across people not reporting ARVMED\n")
    dall[HIV_AND_VL == 1 & is.na(ARVMED), .(VLmean = mean(log(VL_COPIES + 1, 10))), by = "ROUND"]
    cat("mean log10 VL across people reporting ARVMED\n")
    dall[HIV_AND_VL == 1 & !is.na(ARVMED), .(VLmean = mean(log(VL_COPIES + 1, 10))), by = "ROUND"]
    #
}

if (args$run.comm.analysis) 
{

    library(patchwork)
    library(rstanarm)
    library(tidybayes)
    library(bayesplot)
    library(bayestestR)
    library(lme4)
    library(modelsummary)

    source(file.path(gitdir.functions, 'community_level.R'))

    # create output directory
    glm.out.dir <- file.path(vl.out.dir, "glm")

    # get VL data by comm
    vlc <- vl.vlprops.by.comm.gender.loc(dall, write.csv = FALSE)$DT

    # Plot VIRAEMIA by community-round + binom hypothesis testing
    # ___________________________________________________________

    p_full <- vlc |> plot.prev.viraemia.amonghiv.by.comm.round()
    filename <- file.path("220729_longitudinal_PVLNSofHIV_by_comm_sex.pdf")
    ggsave2(p_full, file = filename, LALA = glm.out.dir, w = 9, h = 12)

    # perform test on simple model assuming P of suppression for every ind.
    tmp0 <- simple.binomial.tests.by.community(vlc)
    means <- tmp0[, list(
        FC=c('inland', 'fishing'),
        P_ROUND_SEX=rep(unique(P_ROUND_SEX),times=2)
        ), by=c('SEX', 'ROUND')]
    signif <- tmp0[P_COMB_FISHER <= .05, .(SEX, COMM_NUM)] 

    p_onlysig <- vlc |> 
        merge(signif, by=c('SEX', 'COMM_NUM')) |> 
        plot.prev.viraemia.amonghiv.by.comm.round(MEANS=means) 
    p_onlysig <- vlc |> plot.prev.viraemia.amonghiv.by.comm.round()
    filename <- file.path("220729_longitudinal_PVLNSofHIV_by_comm_sex_onlysignificant.pdf")
    ggsave2(p_onlysig, file = filename, LALA = glm.out.dir, w = 9, h = 12)


    if (0) {
        p_map <- make.map.220810(copy(dall), copy(vlc))
        filename <- file.path("220729_map_PVLNSofHIV_by_comm_sex.pdf")
        ggsave2(p_map, file = filename, w = 9, h = 12)
    }

    # What kind of GLM to use?
    # Well we can start with a binomial right?

    dglm <- get.glm.data(dall)
    dglm[, ROUND := as.factor(ROUND)]
    dglm[, ROUNDi := as.integer(ROUND) - 1]

    # tmp0 <- agresti.coull.exploration(dglm)

    # What's the best way to set priors and hyperprios? Maybe look at 8 schools examples.
    # Also: http://mc-stan.org/rstanarm/articles/glmer.html
    options(mc.cores = parallel::detectCores())
    prior.pars <- list(intercept.mean = dglm[, logit(sum(VLNS_N) / sum(HIV_N))])

    # FIT AND SAVE MODELS
    # ___________________
    
    fit_allinteract <- fit.glm.model(
        formula = VLNS_N / HIV_N ~ (AGEGROUP*SEX_LABEL*FC*ROUND) + (1|COMM_NUM),
        suffix="suppofpart_allinteractions.rds"
    )

    compare.random.effects(fit_allinteract)


    #
    # ANALYSE
    #

    # loo for leave one out cross validation!

    p_slopes <- analyse_glm_model(glm_slopes, "slopes")

    loo(glm_slopes) |> plot()

    glm_model_choice <- copy(glm_5dec)

    get.effective.sample.size.glm(glm_model_choice)

    # help('prior_summary.stanreg')
    prior_summary(glm_model_choice)

    # get a feel
    launch_shinystan(glm_model_choice)
    mcmc_intervals(glm_model_choice)

    summary(glm_model_choice)

    # lala
    glm_model_choice$ses

    .f <- function(vec) data.table(NAME = names(vec), VALUE = vec)
    .f(se(glm_reffs)) |> ggplot(aes(x = NAME, y = VALUE)) +
        geom_point() +
        coord_flip()

    # COMM_NUM random effects sorted by median. Can I color by comm type
    if (0) {
        tmp <- glm_model_choice %>%
            spread_draws(b[term, group], `(Intercept)`) %>%
            tidyr::separate(group, c("RAND_EFF", "COMM_NUM"), ":")
        setDT(tmp)
        tmp <- tmp[COMM_NUM != "FC",
            {
                z <- quantile(b, probs = c(.025, .5, .975))
                list(CL = z[1], M = z[2], CU = z[3])
            },
            by = "COMM_NUM"
        ]

        setorder(tmp, "M")
        performers_bad <- tmp[CU < 0]
        performers_good <- tmp[CL > 0]
        nms <- tmp[, COMM_NUM]
        nms <- paste0("b[(Intercept) COMM_NUM:", nms, "]")

        p_tmp <- mcmc_intervals(x = glm_model_choice, pars = nms)
        .gs <- function(x) as.integer(gsub("[A-z]|\\)|\\(|:", "", x))
        lvls <- .gs(p_tmp$data$parameter)
        p_tmp$data$COMM_NUM <- lvls
        tmp <- unique(dglm[, .(COMM_NUM, FC2), ])
        p_tmp$data <- merge(p_tmp$data, tmp)

        p_ranges <- ggplot(
            p_tmp$data,
            aes(x = m, y = ordered(COMM_NUM, levels = lvls), xmin = ll, xmax = hh, col = FC2)
        ) +
            geom_vline(aes(xintercept = 0), linetype = "dotted", color = "red") +
            geom_point(size = 1.5) +
            geom_errorbarh(height = 0, size = 0.5) +
            geom_errorbarh(aes(xmin = l, xmax = h), height = 0, size = 1) +
            scale_color_manual(values = palettes$comm2) +
            labs(x = "Posterior", y = "Community level random effects", color = "Community type")
        p_ranges

        filename <- file.path("comm_stanglm_diagnostics_comm_random_effects.pdf")
        ggsave2(p_ranges, file = filename, w = 9, h = 12)
    }
}

# Estimate HIV prevalence
# ________________________

if (args$run.gp.prevl) {
    vl.prevalence.by.gender.loc.age.gp(dall, refit = args$refit)
}

# Estimate mean viral load
# ________________________

if (args$run.icar.mean.vl) {
    vl.meanviralload.by.gender.loc.age.icar(dall, refit = args$refit)
}

# Estimate suppressed pop
# _______________________

if (args$run.gp.supp.hiv) { # Among HIV positive
    vl.suppofinfected.by.gender.loc.age.gp(dall, refit = args$refit)
}


if (args$run.gp.supp.pop) { # Among Entire population
    vl.suppofpop.by.gender.loc.age.gp(dall, refit = args$refit)
}


if (0) {
    # GET POSTERIORS ON SUPP AMONG POP
    # ________________________________

    .f <- function(file) {
        tmp <- new.env()
        round <- as.integer(gsub("^.*?round([0-9][0-9]).*?$", "\\1", file))
        load(file, envir = tmp)
        ls(tmp)
        tmp <- tmp$nspop.by.age
        tmp[, ROUND := round]
        tmp
    }
    dsupp <- list.files(vl.out.dir, "220729f_suppAmongPop.*rda", full.names = T)
    dsupp <- lapply(dsupp, .f)
    dsupp <- rbindlist(dsupp)



    # CAN WE GET INCIDENCE NOW? (why?)
    # ________________________

    dinc <- file.path(git.repository, "data", "RCCS_1518_incidence.csv")
    dinc <- fread(dinc)

    dinc[, .(INCIDENCE * PY, NEWINF)]

    tmp <- grep("ROUND|COMM|AGEYRS|SEX|INCIDEN|PREVALEN", names(dinc), value = TRUE)
    dinc <- dinc[, ..tmp]

    p <- ggplot(dinc, aes(x = AGEYRS, colour = ROUND)) +
        geom_line(aes(y = INCIDENCE, group = ROUND)) +
        facet_grid(COMM ~ SEX) +
        viridis::scale_color_viridis() +
        theme_bw()


    # Want to
    # 1. get Mean Viral Load by location
    # 2. run an analysis GLM like
    # 3. Using what as a predictor?
}
