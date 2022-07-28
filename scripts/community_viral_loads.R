# AIMS:

# generalise Oli's phyloscan.viral.load.project.R analyses to a longitudinal setting.


################
# DEPENDENCIES #
################

library(data.table)
library(ggplot2)
library(Hmisc)

################
#    PATHS     #
################

usr <- Sys.info()[['user']]
if(usr == 'andrea')
{
        init()
        indir.repository <-'~/git/longi_viral_loads'
        indir.deepsequence.data <- '~/Documents/Box/ratmann_pangea_deepsequencedata'
        
}else{
        indir.repository <-'~/git/longi_viral_loads'
}

out.dir <- file.path(indir.repository,'results')
path.viral.load <- file.path(indir.deepsequence.data, 'RCCS_R15_R20',
                      "viral_loads_r15r20_processed_220727.csv")
path.r15r18_hiv_tests <- file.path(indir.deepsequence.data, 'RCCS_R15_R18',
                                   "HIV_R15_R18_VOIs_220129.csv")
path.prevalence <- file.path(indir.deepsequence.data, 'RCCS_R15_R18', 'RCCS_census_eligible_count_220719.csv')

file.exists(
        indir.repository,
        indir.deepsequence.data,
        path.viral.load,
        path.r15r18_hiv_tests,
        path.prevalence
) |> all() |> stopifnot()


################
#    HELPERS   #
################

source(file.path(indir.repository, 'functions', 'base_utilities.R'))

.age.aggregate <- function(x, width)
{
        bounds <- seq(15, 50, by=width)
        labs <- paste0(bounds, '-', shift(bounds, -1))
        labs <- head(labs, -1)
        cut(x, bounds, include.lowest=TRUE, labels=labs)
}

pvl.measures.by.region <- function(DT, by_cols=c('comm_num', 'round'))
{
        # MVL = mean viral load (Geometric Mean)
        # PDV = prevalence of detectable viremia
        # CTI = community transmission index (?). [Wilson et. al]

        .mvl <- function(x)
        {
                x <- fifelse(x==0, 1, x)
                exp(mean(log(x)))
        }

        .pdv <- function(x)
                mean(x >= VIREMIC_VIRAL_LOAD)

        .cti <- function(x)
        {
                v0 <- 150 
                beta0 <- 3e-3 
                x_ratio <- x/v0
                beta1 <- 2.45^(log( x_ratio , base=10))*beta0
                beta1[beta1 < 0] <- 0
                cti = 100*(1 - (1-beta1)^100)
                mean(cti)
        }

        .all <- function(x) 
                list(MVL=.mvl(x), PDV=.pdv(x), CTI=.cti(x), N=length(x))

        setkey(DT, region, round)

        tmp <- DT[ hiv_vl > 0, .all(hiv_vl), by= by_cols]
        cols <- setdiff(names(tmp), by_cols)
        setnames(tmp, cols, paste0(cols, '+'))

        tmp2 <- DT[, .all(hiv_vl), by= by_cols ]
        setnames(tmp2, cols, paste0(cols, 'pop'))

        merge(tmp, tmp2, by=by_cols)
}

get.empirical.pdv.by.round.region.sex.age <- function(DT, age_width=5)
{
        # specify years age group
        DT[ ageyrs %between% c(15,50), AGEGROUP:=.age.aggregate(ageyrs, width=age_width), ]

        DT[, DET := as.integer(hiv_vl > VIREMIC_VIRAL_LOAD) ]
        by_cols <- c('round', 'comm', 'sex', 'AGEGROUP')
        cols <- c('PDV_m', 'PDV_l', 'PDV_u')

        tmp <- DT[, .(DET=sum(DET), TOTAL_DET=.N), by=by_cols]
        tmp <- tmp[!is.na(AGEGROUP)]
        tmp1 <- tmp[, list(AGEGROUP='overall', 
                           DET=sum(DET, na.rm=T),
                           TOTAL_DET=sum(TOTAL_DET, na.rm=T)) ,
            by=c('comm', 'round', 'sex')]
        tmp <- rbind(tmp, tmp1)

        by_cols <- c(by_cols, 'DET', 'TOTAL_DET')
        tmp <- tmp[, binconf(DET, TOTAL_DET, return.df=T), by=by_cols]
        setnames(tmp, setdiff(names(tmp), by_cols), cols)


        # empirical Prevalence of detectable viremia
        p <- ggplot(tmp[!is.na(AGEGROUP)],
               aes(x=AGEGROUP, y=PDV_m, colour=sex) ) + 
                geom_point(position=position_dodge(width=0.2)) +
                geom_linerange(aes(ymin=PDV_l, ymax=PDV_u), position=position_dodge(width=0.2)) + 
                facet_grid(round ~ comm) + 
                scale_y_continuous(labels = scales::percent) + 
                scale_colour_manual(values=palettes$sex) + 
                theme_bw() +
                theme( legend.position='bottom') +
                labs(x='age groups', y='prevalence of detectable viremia',
                     title='Prevalence of Detectable Viremia among HIV+, by community, age and sex.')

        filename <- file.path(out.dir.day, 'pdvamongpositive_by_round_location_agegroup.png')
        if(! is.null(filename))    saveplot(filename, p, w=10,h=8)

        # 
       p <-  ggplot(tmp[!is.na(AGEGROUP)]) +
                geom_tile(aes(x=AGEGROUP, y=round, fill=PDV_m)) +
                scale_fill_gradient(low='white', high='red') + 
                scale_x_discrete(expand=c(0,0)) + 
                scale_y_continuous(expand=c(0,0)) +
                facet_grid(sex ~ comm) + 
                theme_bw() +
                theme( legend.position='bottom') +
                labs(x='', y='visit round' , fill='PDV', 
                     title='Prevalence of Detectable Viremia among HIV+, by community, age and sex.')
                
        filename <- file.path(out.dir.day, 'pdvamongpositive_heatmap_by_round_location_agegroup.png')
        if(! is.null(filename))    saveplot(filename, p, w=10,h=8)

        tmp[, round := as.character(round)]
        tmp
}

################
#     MAIN     #
################

# Create output directory with the results of the day
#____________________________________________________
date <- format(Sys.time(), '%y%m%d')
out.dir.day <- file.path(out.dir, paste0(date, '_cvl'))
dir.create(out.dir.day)

# Args
#_____
set.na.vl.to.viremic.threshold <- TRUE
VIREMIC_VIRAL_LOAD = 500

# I: study prevalence
#____________________
cat("\n\n--- I: Study prevalence ----\n\n")
dprev <- fread(path.prevalence)

if(0)
{
        # Compute empirical prevalence with CI
        cols <- c('PREVALENCE_EMPIRICAL','PREVALENCE_EMP_L', 'PREVALENCE_EMP_U' )
        dprev[ !is.na(COUNT) , (cols) := binconf(COUNT, TOTAL_COUNT, return.df=T)]


        # are the GPs maybe overconfident??
        ggplot(dprev, aes(x=AGEYRS, y=PREVALENCE_PROPORTION, col=COMM)) + 
                geom_line() + 
                geom_linerange(aes(ymin=PREVALENCE_EMP_L, ymax=PREVALENCE_EMP_U), position=position_dodge(width=0.2)) + 
                geom_ribbon(aes(ymin=PREVALENCE_PROPORTION_CL, ymax=PREVALENCE_PROPORTION_CU, fill=COMM),
                            alpha=0.2, colour=NA) + 
                facet_grid( ROUND ~ SEX) + 
                scale_colour_manual(values=palettes$comm) + 
                scale_fill_manual(values=palettes$comm) + 
                theme_bw()

        ggplot(dprev[ROUND==17, ], aes(x=AGEYRS, y=PREVALENCE_PROPORTION, col=SEX)) + 
                geom_line() + 
                geom_point(aes(y= PREVALENCE_EMPIRICAL), linetype='dotted') + 
                geom_linerange(aes(ymin=PREVALENCE_EMP_L, ymax=PREVALENCE_EMP_U), position=position_dodge(width=0.2)) + 
                geom_ribbon(aes(ymin=PREVALENCE_PROPORTION_CL, ymax=PREVALENCE_PROPORTION_CU, fill=SEX),
                           alpha=0.2, colour=NA) + 
                facet_grid( ROUND ~ COMM) + 
                scale_colour_manual(values=palettes$sex) + 
                scale_fill_manual(values=palettes$sex) + 
                theme_bw()
}

# also obtain a version of prevalence stratified by 5 y age groups

by_cols <- c('COMM', 'ROUND', 'SEX', 'AGEYRS', 'COUNT', 'TOTAL_COUNT' )
dprev2 <- dprev[, ..by_cols]
dprev2 <- dprev2[, AGEGROUP := .age.aggregate(AGEYRS, width=5)]


by_cols <- c('COMM', 'ROUND', 'SEX', 'AGEGROUP')
var_cols <- c('COUNT', 'TOTAL_COUNT')
dprev2 <- dprev2[, lapply(.SD,sum, na.rm=T) , .SDcols=var_cols, by=by_cols]
tmp <- dprev2[, list(AGEGROUP='overall', 
                     COUNT=sum(COUNT, na.rm=T),
                     TOTAL_COUNT=sum(TOTAL_COUNT, na.rm=T)) ,
              by=c('COMM', 'ROUND', 'SEX')]
dprev2 <- rbind(dprev2, tmp)

by_cols <- c(by_cols, var_cols)
var_cols <- c('PREVALENCE_EMPIRICAL','PREVALENCE_EMP_L', 'PREVALENCE_EMP_U' )
dprev2[ , (var_cols) := binconf(sum(COUNT), sum(TOTAL_COUNT), return.df=T), by=by_cols]

if(0)
{
        names(dprev2)

        # are the GPs maybe overconfident??
        p <- ggplot(dprev2, aes(x=AGEGROUP, col=SEX)) + 
                geom_point(aes(y=PREVALENCE_EMPIRICAL), position=position_dodge(width=0.2)) +
                geom_linerange(aes(ymin=PREVALENCE_EMP_L, ymax=PREVALENCE_EMP_U), position=position_dodge(width=0.2)) + 
                facet_grid( ROUND ~ COMM) + 
                scale_y_continuous(labels = scales::percent) + 
                scale_colour_manual(values=palettes$sex) + 
                scale_fill_manual(values=palettes$sex) + 
                theme_bw() +
                theme( legend.position='bottom') +
                labs(x='', y='HIV prevalence', colour='sex', 
                     title='Prevalence by sex, age, round and community type')

        filename <- file.path(out.dir.day, 'prevalence_by_round_location_agegroup.png')
        if(! is.null(filename))    saveplot(filename, p, w=10,h=8)

        p <- ggplot(dprev2) +
                geom_tile(aes(x=AGEGROUP, y=ROUND, fill=PREVALENCE_EMPIRICAL)) +
                scale_fill_gradient(low='white', high='red') + 
                scale_x_discrete(expand=c(0,0)) + 
                scale_y_discrete(expand=c(0,0)) +
                facet_grid(SEX ~ COMM) + 
                theme_bw() +
                theme( legend.position='bottom') +
                labs(y='Visit round', x='Age group', fill='HIV Prevalence', title='Empirical prevalence by sex, age, round and community type')

        filename <- file.path(out.dir.day, 'prevalence_heatmap_by_round_location_agegroup.png')
        if(! is.null(filename))    saveplot(filename, p, w=10,h=8)
}

# II: study suppression among HIV-positive
# want to do this by age and round

dvl <- fread(path.viral.load)
dvl <- dvl[round %between% c(16, 20)]

if(set.na.vl.to.viremic.threshold)
{
        cat('Setting Unknown VLs equal to viremic threshold...\n')
        dvl[is.na(hiv_vl), hiv_vl := VIREMIC_VIRAL_LOAD]
}


dpdv <- get.empirical.pdv.by.round.region.sex.age(dvl)
names(dpdv) <- toupper(names(dpdv))


# III combine I, II to obtain proportion of population with unsuppressed viral load.

# TODO: check.
# It is not necessarily the case that Numerators of the prevalence calculations
# are larger than the denominators of the DV pop...
# WHAT TO DO?

cols <- c('COMM', 'ROUND', 'SEX', 'AGEGROUP')
setkey(dpdv, COMM, ROUND, SEX, AGEGROUP)
setkey(dprev2, COMM, ROUND, SEX, AGEGROUP)
setcolorder(dpdv, cols)
setcolorder(dprev2, cols)

cols1 <- c(cols, 'DET', 'TOTAL_DET')
cols2 <- c(cols, 'COUNT', 'TOTAL_COUNT')

merge(
        dpdv[, ..cols1],
        dprev2[, ..cols2],
        by=cols
) -> tmp

tmp

# How to combine the uncertainty in the two estimates?
# well, we can assume that the 2 experiments are independent
# estimating the probabilities p1 and p2.
# Now, p1 ~ Beta(x1, n1-x1) and p2 ~ Beta(x2, n2-x2)
# so we can estimate p1*p2 by MC:

tmp[, {
        p1 <- rbeta(10000, shape1=DET, shape2=TOTAL_DET-DET)
        p2 <- rbeta(10000, shape1=COUNT, shape2=TOTAL_COUNT-COUNT)

        p <- c(2.5, 50, 97.5)
        tmp <- as.data.table(quantile(p1 * p2, probs=p/100))
        tmp <- transpose(tmp)
        names(tmp) <- c('CL', 'M', 'CU')
        tmp
},  by = cols] -> tmp

p <- ggplot(tmp, aes(x=AGEGROUP, y=M,  col=SEX)) + 
        geom_point(position=position_dodge(width=.2)) + 
        geom_linerange(aes(ymin=CL, ,ymax=CU), position=position_dodge(width=.2)) + 
        facet_grid( ROUND ~ COMM) + 
        scale_y_continuous(labels = scales::percent) + 
        scale_colour_manual(values=palettes$sex) + 
        theme_bw() +
        theme( legend.position='bottom') +
        labs(x='', y='UVL prevalence', colour='sex',
             title='Proportion of population with Unsuppressed Viral Load')

filename <- file.path(out.dir.day, 'unsuppressedvl_by_round_location_agegroup.png')
if(! is.null(filename))    saveplot(filename, p, w=10,h=8)

.f <- function(DT, m, l, u, nm='V1')
{
        .r <- function(x) as.character(round(x*100, 2))
        .g <- function(x) paste0( .r(x[[1]]), ' (', .r(x[[2]]), '-', .r(x[[3]]), ')')  

        cols <- c('COMM', 'ROUND', 'SEX', 'AGEGROUP')
        sdcols <- c(m, l, u)
        out <- DT[ AGEGROUP == 'overall', .g(.SD), by=cols, .SDcols=sdcols ]
        setnames(out, 'V1', nm)
        out
}

list(
        .f(tmp, 'M', 'CL','CU', nm='U-VL'),
        .f(dpdv, 'PDV_M', 'PDV_L','PDV_U', nm='PDV'),
        .f(dprev2, 'PREVALENCE_EMPIRICAL', 'PREVALENCE_EMP_L','PREVALENCE_EMP_U', nm='HIV-prev')
) -> tbl
tbl <- Reduce(merge,tbl)
tbl[, AGEGROUP := NULL]
tbl
setkey(tbl, ROUND, COMM, SEX, `HIV-prev`, PDV, `U-VL`)
setcolorder(tbl, c('ROUND', 'COMM', 'SEX', 'HIV-prev', 'PDV', 'U-VL'))
setnames(tbl,
         c('ROUND', 'COMM', 'SEX', 'HIV-prev', 'PDV', 'U-VL'),
         c('Round', 'Location', 'Gender', 'HIV prevalence (%)', 'Unsuppressed among infected (%)', 'Unsuppressed in population (%)'))
tbl[Round == 17]
tbl1 <- knitr::kable(tbl)

filename <- file.path(out.dir.day, 'summary_table.tex')
if(! is.null(filename))  print(xtable::xtable(tbl), file=filename)

?binconf
