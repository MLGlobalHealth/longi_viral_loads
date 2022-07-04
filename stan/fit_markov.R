# AIMS:

# Load visit pairs and fit the model


################
# DEPENDENCIES #
################

library(cmdstanr)
library(data.table)
library(ggplot2)
library(loo)
library(bayesplot)

################
#    PATHS     #
################

usr <- Sys.info()[['user']]
if(usr == 'andrea')
{
        indir.repository <-'~/git/longi_viral_loads'
        indir.deepsequence.data <- '~/Documents/Box/ratmann_pangea_deepsequencedata'
        indir.deepsequence.analyses   <- '~/Documents/Box/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI'
}else{
        indir.repository <-'~/git/longi_viral_loads'
        indir.deepsequence.data <- '/rds/general/project//ratmann_pangea_deepsequencedata/live'
        indir.deepsequence.analyses   <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI'
}

out.dir <- file.path(indir.repository,'results')
file.viral.loads.processed <- file.path(indir.deepsequence.data,
                                        'RCCS_R15_R20',
                                        "viral_loads_r15r20_processed_220614.csv")



################
#    HELPERS   #
################

source(file.path(indir.repository, 'functions', 'base_utilities.R'))
source(file.path(indir.repository, 'functions', 'preprocessing_helpers.R'))

make_standata <- function(dvl, THRESHOLD = 50)
{
        cat('Chosen threshold of', THRESHOLD, 'copies/mL\n\n')

        tmp <- dvl[, .(study_id, round, sex, hivdate, hiv_vl)]
        setkey(tmp, study_id, hivdate)

        tmp[, cat( round(100 * mean(is.na(hiv_vl)), 2), '% of measurements are NA, remove them...\n')]
        tmp <- tmp[!is.na(hiv_vl)]

        tmp[, v_start := as.integer(hiv_vl > THRESHOLD)]
        tmp[, v_final := shift(v_start, -1), by=study_id]
        tmp[, round_final := shift(round, -1), by=study_id]
        tmp[, cat(sum(!is.na(v_final)), 'visit-pairs obtained out of ', .N, 'visits.\n\n')]
        tmp <- tmp[!is.na(v_final)]

        with(
             tmp,
             list(
                  NP = length(v_start),
                  v_start = v_start,
                  v_final = v_final,
                  zF = as.integer(sex == 'F')
             )
        ) -> stan_data

        list(stan=stan_data, R=tmp)
}

################
#     MAIN     #
################

stan.file <- 'markov_transitions_220614.stan'
tmp <- gsub('.stan', '', stan.file)
out.dir <- file.path(out.dir, tmp)
dir.create(out.dir)

# Load data
dvl <- fread(file.viral.loads.processed)
tmp <- make_standata(dvl)
stan_data <- tmp$stan
dvl_analysis <- tmp$R

# Load file and run
mod <- cmdstan_model(stan.file)

fit <- mod$sample(stan_data)
print(fit)

pdraws <- fit$draws(format='array', inc_warmup=FALSE)
params <- dimnames(pdraws)[["variable"]]

select.chains <- seq_along(dimnames(pdraws)[['chain']])
tot.iters <- prod(dim(pdraws)[c(1,2)])

pos <- list()
{
	# bring samples into long array format (change for beta, cause it s a matrix)
        tmp <- pdraws[, select.chains, grepl('pi\\[[0-9],[0-9]]', params)]
        tmp <- apply(tmp, 3, rbind)

        pos$piM0 <- tmp[, dimnames(tmp)$variable == 'pi[1,1]']
        pos$piF0 <- tmp[, dimnames(tmp)$variable == 'pi[2,1]']
        pos$piM1 <- tmp[, dimnames(tmp)$variable == 'pi[1,2]']
        pos$piF1 <- tmp[, dimnames(tmp)$variable == 'pi[2,2]']
}
pos


# Plots by round of first visit 
#_______________________________

.plot.probabilities.next.viremic.by.sex <- function( which.round = '' )
{

        by.var <- switch(which.round,
                start = list(col='round',       xlab='Round of first visit in the visit-pair'),
                final = list(col='round_final', xlab='Round of second visit in the visit-pair')
        )
        if(is.null(by.var)) stop('which.round should take either of "start" or "final" as values.\n')

        tmp <- dvl_analysis[, 
                     list(
                          N = length(v_final),
                          S = sum(v_final)
                     )
                     , by = c('sex', 'v_start', by.var[[1]])]

        tmp1 <- tmp[, Hmisc::binconf(x=S, n=N)]
        tmp <- cbind(tmp, tmp1)

        # Transform round to factor
        tmp[, (by.var[[1]]) := lapply(.SD, round2factor) , .SDcols=by.var[[1]]]
        tmp[, first_viremic := (v_start==1) ]

        .f <- function(x){quantile(x, probs = c(0.025, 0.5, .975))}
        tmp1 <- lapply(pos, .f)
        tmp1 <- do.call(rbind, tmp1)
        tmp1 <- cbind(rownames(tmp1),as.data.table(tmp1))
        names(tmp1)[-1] <- c('L', 'M', 'U')
        tmp1[, sex:=gsub('[a-z]|[0-1]', '', V1)]
        tmp1[, first_viremic:=as.integer(gsub('[A-z]', '', V1))]
        tmp1[, first_viremic:=as.logical(first_viremic)]


        p <- ggplot(tmp, aes_string(x = by.var[[1]], shape='first_viremic' , linetype='first_viremic', color='sex')) +
                geom_point(aes(y = PointEst), position=position_dodge(width=.5), size=2) + 
                geom_linerange(aes(ymin=Lower, ymax=Upper),
                               position=position_dodge(width=.5)) + 
                geom_hline(data=tmp1, aes(yintercept=M, color=sex)) +
#        geom_rect(data=tmp1, aes(ymin=L, ymax=U, group=first_viremic),
#                  xmin=as.numeric(tmp$round[[2]]) - .3,
#                  xmax=as.numeric(tmp$round[[1]]) + 0.3) + 
                scale_y_continuous(labels = scales::percent, 
                                   expand = c(0, 0), limits = c(0, 1)) +
                scale_color_manual(values=palettes$sex) + 
                theme_bw() +
                theme(legend.position='bottom') + 
                labs(x = by.var[[2]],
                     y = 'Proportion of viremic follow-ups', 
                     subtitle = 'Posterior means (horizontal lines) and empirical data (point + lineranges)',
                     title = 'How many follow-up measurements are viremic? ') 
                p
}

p <- .plot.probabilities.next.viremic.by.sex(which.round='start')
filename = file.path(out.dir, 'proportion_viremic_by_roundfirst_sex.png')
saveplot(filename, p, w=8, h=6)

p <- .plot.probabilities.next.viremic.by.sex(which.round='final')
filename = file.path(out.dir, 'proportion_viremic_by_roundsecond_sex.png')
saveplot(filename, p, w=8, h=6)

# Plots by round difference
#__________________________

dvl_analysis[round_final < round, study_id] -> idx

dvl[study_id %in% idx]


tmp <- copy(dvl_analysis)

tmp[, round_diff := round_final - round ]

        tmp <- tmp[, list(
                          N = length(v_final),
                          S = sum(v_final)
                     )
                     , by = c('sex', 'v_start', 'round_diff')]
        tmp1 <- tmp[, Hmisc::binconf(x=S, n=N)]
        tmp <- cbind(tmp, tmp1)

        # Transform round to factor
        # tmp[, (by.var[[1]]) := lapply(.SD, round2factor) , .SDcols=by.var[[1]]]
        tmp[, first_viremic := (v_start==1) ]

        .f <- function(x){quantile(x, probs = c(0.025, 0.5, .975))}
        tmp1 <- lapply(pos, .f)
        tmp1 <- do.call(rbind, tmp1)
        tmp1 <- cbind(rownames(tmp1),as.data.table(tmp1))
        names(tmp1)[-1] <- c('L', 'M', 'U')
        tmp1[, sex:=gsub('[a-z]|[0-1]', '', V1)]
        tmp1[, first_viremic:=as.integer(gsub('[A-z]', '', V1))]
        tmp1[, first_viremic:=as.logical(first_viremic)]
        tmp[, round_diff := as.factor(round_diff)]


        p <- ggplot(tmp, aes_string(x = 'round_diff', shape='first_viremic' , linetype='first_viremic', color='sex')) +
                geom_point(aes(y = PointEst), position=position_dodge(width=.5), size=2) + 
                geom_linerange(aes(ymin=Lower, ymax=Upper),
                               position=position_dodge(width=.5)) + 
                geom_hline(data=tmp1, aes(yintercept=M, color=sex)) +
#        geom_rect(data=tmp1, aes(ymin=L, ymax=U, group=first_viremic),
#                  xmin=as.numeric(tmp$round[[2]]) - .3,
#                  xmax=as.numeric(tmp$round[[1]]) + 0.3) + 
                scale_y_continuous(labels = scales::percent, 
                                   expand = c(0, 0), limits = c(0, 1)) +
                scale_color_manual(values=palettes$sex) + 
                theme_bw() +
                theme(legend.position='bottom') + 
                labs(x = 'Round difference',
                     y = 'Proportion of viremic follow-ups', 
                     subtitle = 'Posterior means (horizontal lines) and empirical data (point + lineranges)',
                     title = 'How many follow-up measurements are viremic? ') 
                p



