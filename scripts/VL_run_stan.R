# AIMS:
# - adapt older code from Oli to multi-round, longitudinal settings
# TODO: check why we do not have any ARVMED == 2
# TODO: discuss: we are removing individuals with missing VLs: they are very little

################
# DEPENDENCIES #
################

library(data.table) 
library(ggplot2)
library(Hmisc)
library(rstan)
# For parallelisation across cores:
# library(foreach)
# library(doParallel)

################
#    PATHS     #
################

parallelise <- FALSE

usr <- Sys.info()[['user']]
if(usr == 'andrea')
{
    git.repository <-'~/git/longi_viral_loads'
    indir.deepsequence.data <- '~/Documents/Box/ratmann_pangea_deepsequencedata'
    indir.deepanalyses.xiaoyue <- '/home/andrea/Documents/Box/ratmann_xiaoyue_jrssc2022_analyses/live'
    # out.dir.prefix <- '/media/andrea/SSD/2022/longvl'
    # out.dir.prefix <- file.path(git.repository, 'results')
    # parallelise <- TRUE
    # 
}

if(usr == 'ab1820')
{
    git.repository <-'/rds/general/user/ab1820/home/git/longi_viral_loads'
    indir.deepsequence.data <- '/rds/general/project/ratmann_pangea_deepsequencedata/live'
    indir.deepanalyses.xiaoyue <- '/rds/general/project/ratmann_xiaoyue_jrssc2022_analyses/live'
    # out.dir.prefix <- '/home/andrea/Documents/Box/ratmann_xiaoyue_jrssc2022_analyses/live'
}

path.stan <- file.path(git.repository, 'stan')
path.tests <- file.path(indir.deepsequence.data, 
                        'RCCS_R15_R20',
                        "all_participants_hivstatus_vl_220729.csv")
file.exists(
            # out.dir,
            path.stan,
            path.tests
) |> all() |> stopifnot()

################
#   OPTIONS    #
################

option_list <- list(
    optparse::make_option(
        "--refit",
        type = "logical",
        default = FALSE,
        help = "Flag on whether to re-run the stan models even if already exist [Defaults to FALSE]", 
        dest = 'refit'
    ),
    optparse::make_option(
        "--run-gp-prevl",
        type = "logical",
        default = FALSE,
        help = "Flag on whether to run the Stan model for prevalence [Defaults to FALSE] ",
        dest = 'run.gp.prevl'
    ),
    optparse::make_option(
        "--run-icar-mean-vl",
        type = "logical",
        default = FALSE,
        help = "Flag on whether to run the Stan model for population level viral loads [Defaults to FALSE] ",
        dest = 'run.icar.mean.vl'
    ),
    optparse::make_option(
        "--run-gp-supp-hiv",
        type = "logical",
        default = FALSE,
        help = "Flag on whether to run the Stan model for suppression among HIV positive [Defaults to FALSE] ",
        dest = 'run.gp.supp.hiv'
    ),
    optparse::make_option(
        "--run-gp-supp-pop",
        type = "logical",
        default = FALSE,
        help = "Flag on whether to run the Stan model for suppression among population [Defaults to FALSE] ",
        dest = 'run.gp.supp.pop'
    ),
    optparse::make_option(
        "--run-comm-analysis",
        type = "logical",
        default = FALSE,
        help = "Flag on whether to run the Stan GLM model for community-level suppression [Defaults to FALSE] ",
        dest = 'run.comm.analysis'
    ),
    optparse::make_option(
        "--viremic-viral-load",
        type = "numeric",
        default = 1000,
        help = "Duration of acute phase, in years [Defaults to 2 months]", 
        dest = 'viremic.viral.load'
    ),
    optparse::make_option(
        "--vl-detectable",
        type = "numeric",
        default = 150,
        help = "Duration of acute phase, in years [Defaults to 2 months]", 
        dest = 'viremic.viral.load'
    ),
    optparse::make_option(
        "--outdir",
        type = "character",
        default = '~/OneDrive/2022/longvl/220729_oli',
        help = "Path to output directory where to store results [Defaults to my directory]", 
        dest = 'out.dir.prefix'
    ),
    optparse::make_option(
        "--round",
        type = "numeric",
        default = c(16,17,18,19),
        help = "Path to output directory where to store results [Defaults to my directory]", 
        dest = 'round'
    )
)

args <-  optparse::parse_args(optparse::OptionParser(option_list = option_list))
print(args)

################
#    HELPERS   #
################

source( file.path(git.repository,'functions/base_utilities.R') )
source( file.path(git.repository,'scripts/phsc_vl_helpers.R'))

if(parallelise)
{   # set up parallel backend
    n.cores <- min(4, parallel::detectCores() - 1 )

    logname <- Sys.time() |> format('%m%d_%H%m')
    logname <- paste0('.parallel_log_', logname )
    my.cluster <- parallel::makeCluster(
            n.cores,
            type='FORK',
            outfile=logname)
    doParallel::registerDoParallel(cl = my.cluster)
    print(my.cluster)
}


################
#     MAIN     #
################

# check study round exists
stopifnot(all(args$round %in% 16:19))

# set viral load thresholds
VL_DETECTABLE = args$vl.detectable
VIREMIC_VIRAL_LOAD = args$viremic.viral.load

# specify and create output directories
stopifnot( dir.exists(args$out.dir.prefix))
out.dir <- file.path(args$out.dir.prefix)
vl.out.dir <- out.dir
if(usr=='andrea')
    vl.out.dir <- file.path(out.dir, paste0('vl_', VIREMIC_VIRAL_LOAD) )
dir.create(vl.out.dir, showWarnings = FALSE)

# get data
dall <- get.dall(path.tests)

if(0) 
{   # Info for introduction to results

    .mean2 <- function(x) paste0(round(100*mean(x), 2), '%')

    # proportions of viraemic measurements
    dall[VL_COPIES > VIREMIC_VIRAL_LOAD, .mean2(SEX=='F'), ]
    dall[VL_COPIES > VIREMIC_VIRAL_LOAD, .mean2(SEX=='F'), by='ROUND']
    dall[VL_COPIES > VIREMIC_VIRAL_LOAD, .mean2(FC=='inland'), by='ROUND']

    # ARVs 
    dall[HIV_AND_VL==1, .mean2(is.na(ARVMED)), by='ROUND']
    cat('mean log10 VL across people not reporting ARVMED\n')
    dall[HIV_AND_VL==1 & is.na(ARVMED), .(VLmean=mean(log(VL_COPIES+1, 10))) , by='ROUND']
    cat('mean log10 VL across people reporting ARVMED\n')
    dall[HIV_AND_VL==1 & !is.na(ARVMED), .(VLmean=mean(log(VL_COPIES+1, 10))) , by='ROUND']
    # 
}

if(0) # Study ARVMED
{
    darv <- dall[HIV_STATUS == 1]
    
    cat('Assume NA ARV means no ARV usage...\n')
    tmp <- darv[ HIV_AND_VL == 1]
    tmp[is.na(ARVMED) | ARVMED != 1, ARVMED := 0]
    by_cols <- c('ROUND', 'FC', 'SEX', 'ARVMED')
    cols <- c('M', 'CL', 'CU')
    tmp[, Y := as.integer(VL_COPIES <= VIREMIC_VIRAL_LOAD)]
    tmp <- tmp[, binconf( sum(Y) , .N, return.df=T), by=by_cols]
    names(tmp) <- c(by_cols, cols)
    setkeyv(tmp, by_cols)
    supp.prop.by.arv <- copy(tmp)
    
    ggplot(supp.prop.by.arv, aes(x=FC, colour=SEX)) + 
            geom_point(aes(y=M, pch=as.factor(ARVMED)), position=position_dodge(width=.5) ) +
            geom_linerange(aes(ymin=CL, ymax=CU, linetype=as.factor(ARVMED)), position=position_dodge(width=.5) ) +
            facet_grid( ~ ROUND) + 
            scale_y_continuous(labels=scales:::percent, limits=c(0,1), expand=c(0,0)) +
            scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
            theme(legend.position='bottom') + 
            labs(x='Community type', y='Proportion of suppressed measurements',
                 linetype='Ever reported ARV', pch='Ever reported ARV',
                 title='Suppression by ARV reporting...')


    cat('Suppression among participants reporting ARV \n')
    darv[ARVMED == 1 & is.na(VL_COPIES), STUDY_ID ] -> idx
    tmp <- darv[STUDY_ID %in% idx, any(VL_COPIES==0, na.rm=T) , by='STUDY_ID']
    tmp[, cat('Out of', .N, 'HIV + participants with NA VL measurements', sum(V1), 
              'were measured suppressed at least once\n')]

    tmp <- darv[HIV_AND_VL == 1 & ARVMED == 1]
    by_cols <- c('ROUND', 'FC', 'SEX')
    cols <- c('M', 'CL', 'CU')
    tmp[, Y := as.integer(VL_COPIES <= VIREMIC_VIRAL_LOAD)]
    tmp <- tmp[, binconf( sum(Y) , .N, return.df=T), by=by_cols]
    names(tmp) <- c(by_cols, cols)
    setkeyv(tmp, by_cols)
    supp.prop.among.report <- copy(tmp)

    # Sex plays a bigger role than community.
    ggplot(supp.prop.among.report, aes(x=FC, colour=SEX)) + 
            geom_point(aes(y=M), position=position_dodge(width=.5) ) +
            geom_linerange(aes(ymin=CL, ymax=CU), position=position_dodge(width=.5) ) +
            facet_grid( ~ ROUND) + 
            scale_y_continuous(labels=scales:::percent, limits=c(.5,1), expand=c(0,0)) +
            scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
            theme(legend.position='bottom') + 
            labs(x='Community type', y='Proportion of suppressed measurements', title='Among people reporting ever ARV...')

}

if(args$run.comm.analysis)
{
    library(patchwork)
    library(rstanarm)
    library(tidybayes)
    library(bayesplot)
    library(bayestestR)
    library(lme4)
    # tmp$p for plot and tmp$DT for 'vlc' datatable
    tmp <- vl.vlprops.by.comm.gender.loc(dall, write.csv=FALSE)
    tmp$DT
    tmp$p_among_pop + tmp$p_among_inf

    vlc <- copy(tmp$DT)
    p <- ggplot(vlc, aes(x=ROUND)) +
        # scale_x_continuous(labels=scales:::percent) +
        scale_y_continuous(labels=scales:::percent) +
        # geom_linerange(aes(ymin=PVLNSofHIV_CL, ymax=PVLNSofHIV_CU), alpha=0.2) +
        # geom_errorbarh(aes(y=PVLNS_MEAN, xmin=PHIV_CL, xmax=PHIV_CU), alpha=0.2) +
        geom_line(aes( y=PVLNSofHIV_MEAN, colour=FC, group=COMM_NUM)) +
        geom_text(aes( y=PVLNSofHIV_MEAN, label=COMM_NUM), size=4) +
        facet_grid(~SEX) +
        scale_colour_manual(values=palettes$comm) + 
        theme_bw() +
        theme(legend.position='bottom') + 
        labs(x='Survey Round', 
             y='proportion unsuppressed HIV among infected\n', 
             colour='location')

    filename <- file.path('220729_longitudinal_PVLNSofHIV_by_comm_sex.pdf')
    ggsave2(p, file=filename, w=9, h=12)

    if(0)
    {
        p_map <- make.map.220810(copy(dall), copy(vlc))
        filename <- file.path('220729_map_PVLNSofHIV_by_comm_sex.pdf')
        ggsave2(p_map, file=filename, w=9, h=12)
    }

    # What kind of GLM to use? 
    # Well we can start with a binomial right? 

    inverse.logit <- function(x){ exp(x)/(exp(x) + 1)}
    dglm <- get.glm.data(dall)
    dglm[, ROUND := as.factor(ROUND)]

    #cols <- c("COMM_NUM","ROUND","SEX", "FC","N","PHIV_MEAN","PHIV_CL","PHIV_CU","PVLNS_MEAN",
    #          "PVLNS_CL","PVLNS_CU","PVLNSofHIV_MEAN","PVLNSofHIV_CL","PVLNSofHIV_CU","FC2")
    #dglm <- subset(tmp$DT,select=cols)
    #.f <- as.integer
    #dglm[, `:=` (N_HIV=.f(N*PHIV_MEAN), N_HIV_NOTSUP=.f(N*PHIV_MEAN*PVLNSofHIV_MEAN))]

    str(dglm)
    comm_lvls <- dglm[, sort(unique(COMM_NUM)), by='FC2'][, V1, ]

    # glm_no_random_effects(PVLNSofHIV_MEAN ~ FC:ROUND)
    # glm_no_random_effects(PVLNSofHIV_MEAN ~ FC:ROUND:SEX)
    # glm_no_random_effects(PVLNSofHIV_MEAN ~ FC:ROUND + SEX) 
    # glm_random_effects(PVLNSofHIV_MEAN ~ FC:ROUND + SEX + (1|COMM_NUM) ) |> AIC()
    # glm_random_effects(PVLNSofHIV_MEAN ~ FC:SEX + ROUND + (1|COMM_NUM) ) |> AIC()

    prior.pars <- list( intercept.mean=logit(dglm[, sum(N_HIV_NOTSUP)/sum(N_HIV)]) )
    
    dglm[, ROUND_IDX := ROUND - min(ROUND + 1)]
    # What's the best way to set priors and hyperprios? Maybe look at 8 schools examples.
    # Also: http://mc-stan.org/rstanarm/articles/glmer.html
    names(dglm)
    options(mc.cores=parallel::detectCores())
    glm_reffs <- stan_glmer(data=dglm,
                            formula=VLNS_N/HIV_N ~ FC:ROUND + SEX + (1|COMM_NUM) + AGEYRS,
                            weights=HIV_N,
                            # prior_intercept = student_t(df = 7, prior.pars$intercept.mean, 5),
                            prior_intercept = normal( prior.pars$intercept.mean, 5),
                            family=binomial(link='logit'))
    # save fit:
    filename <- file.path(vl.out.dir, 'supppofpart_glmer_fit.rds')
    saveRDS(glm_reffs, filename)

    # get sample size
    tmp <- as.data.frame(summary(glm_reffs)) |> as.data.table(keep.rownames = TRUE)
    tmp <- tmp[! rn %in% 'log-posterior', {z1 <- which.min(n_eff); z2 <- which.max(Rhat); list(n=n_eff[z1], pn=rn[z1], R=Rhat[z2], pR=rn[z2])}]
    cat('Minimum effective sample size:',tmp$n, 'for parameter', tmp$pn, '\n')
    cat('Maximum Rhat:',tmp$R, 'for parameter', tmp$pR, '\n')

    help('prior_summary.stanreg')
    prior_summary(glm_reffs)

    # get a feel
    # launch_shinystan(glm_reffs)
    mcmc_intervals(glm_reffs)

    # COMM_NUM random effects sorted by median. Can I color by comm type 
    tmp <- glm_reffs %>%
        spread_draws(b[term,group], `(Intercept)`) %>%
        tidyr::separate(group, c('RAND_EFF', 'COMM_NUM'), ':')
    setDT(tmp)
    tmp <- tmp[, {z <- quantile(b, probs=c(.025, .5, .975)); list(CL=z[1], M=z[2], CU=z[3]) }, by='COMM_NUM' ]
    setorder(tmp, 'M')
    performers_bad <- tmp[CU < 0]
    performers_good <- tmp[CL > 0]
    nms <- tmp[, COMM_NUM]
    nms <- paste0('b[(Intercept) COMM_NUM:',nms,']')

    p_tmp <- mcmc_intervals(x=glm_reffs, pars=nms) 
    .gs <- function(x) as.integer(gsub('[A-z]|\\)|\\(|:','',x))
    lvls <- .gs(p_tmp$data$parameter)
    p_tmp$data$COMM_NUM <- lvls
    tmp <- unique(dglm[, .(COMM_NUM, FC2),])
    p_tmp$data <- merge(p_tmp$data, tmp)

    p_ranges <- ggplot(p_tmp$data,
                       aes(x=m, y=ordered(COMM_NUM, levels=lvls), xmin=ll, xmax=hh, col=FC2)) +
            geom_vline(aes(xintercept=0), linetype='dotted', color='red')+
            geom_point(size=1.5) +
            geom_errorbarh(height=0, size=0.5) +
            geom_errorbarh(aes(xmin=l, xmax=h), height=0, size=1) + 
            scale_color_manual(values=palettes$comm2) +
            labs(x='Posterior', y='Community level random effects', color='Community type') 
    p_ranges

    filename <- file.path( 'comm_stanglm_diagnostics_comm_random_effects.png')
    ggsave2(p_ranges, file=filename, w=9, h=12)


    # POSTERIOR PREDICTIVE CHECKS
    # https://mc-stan.org/rstanarm/reference/plot.stanreg.html
    bayesplot::available_ppc()
    p_check <- pp_check(glm_reffs)
    filename <- file.path('comm_stanglm_diagnostics_ppcheck_density.png')
    ggsave2(p_check, file=filename, w=9, h=8)

    pp_check(glm_reffs, plotfun = "boxplot", nreps = 10, notch = FALSE)
    pp_check(glm_reffs, plotfun = "scatter_avg") # y vs. average yrep
    pp_check(glm_reffs, plotfun = "hist", nreps=3) # y vs. average yrep

    # Posterior vs prior
    p <- posterior_vs_prior( glm_reffs) +
        guides(color=FALSE) + theme(axis.text.x = element_text(angle = 90))
    filename <- file.path( 'comm_stanglm_diagnostics_postvsprior.png')
    ggsave2(p, file=filename, w=9, h=6)

    # Statistical 'significance': probability of direction
    p_direction(glm_reffs)

}

# Estimate HIV prevalence
#________________________

if(args$run.gp.prevl) 
    vl.prevalence.by.gender.loc.age.gp(dall, refit=args$refit)

# Estimate mean viral load
# ________________________

if(args$run.icar.mean.vl) 
    vl.meanviralload.by.gender.loc.age.icar(dall, refit=args$refit)

# Estimate suppressed pop
# _______________________

if(args$run.gp.supp.hiv)   # Among HIV positive
    vl.suppofinfected.by.gender.loc.age.gp(dall, refit=args$refit)


if(args$run.gp.supp.pop)   # Among Entire population
    vl.suppofpop.by.gender.loc.age.gp(dall, refit=args$refit)


if(0)
{
# GET POSTERIORS ON SUPP AMONG POP
# ________________________________

.f <- function(file)
{
        tmp <- new.env()
        round <- as.integer(gsub('^.*?round([0-9][0-9]).*?$', '\\1',file))
        load(file, envir=tmp)
        ls(tmp)
        tmp <- tmp$nspop.by.age
        tmp[, ROUND:=round]
        tmp
}
dsupp <- list.files(vl.out.dir, '220729f_suppAmongPop.*rda', full.names=T)
dsupp <- lapply(dsupp, .f)
dsupp <- rbindlist(dsupp)



# CAN WE GET INCIDENCE NOW? (why?)
# ________________________

dinc <- file.path(git.repository, 'data', 'RCCS_1518_incidence.csv')
dinc <- fread(dinc)

dinc[, .(INCIDENCE*PY, NEWINF)]

tmp <- grep( 'ROUND|COMM|AGEYRS|SEX|INCIDEN|PREVALEN',names(dinc), value=TRUE)
dinc <- dinc[, ..tmp]

p <- ggplot(dinc, aes(x=AGEYRS, colour=ROUND)) + 
        geom_line(aes(y=INCIDENCE, group=ROUND)) + 
        facet_grid(COMM ~ SEX) + 
        viridis::scale_color_viridis() +
        theme_bw() 


# Want to
# 1. get Mean Viral Load by location
# 2. run an analysis GLM like
# 3. Using what as a predictor? 
}
