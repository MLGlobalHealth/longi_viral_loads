################
# DEPENDENCIES #
################

library(data.table)
library(ggplot2)
library(ggtext)
library(ggpubr)
library(knitr)
library(Hmisc)
library(xtable)

################
#    PATHS     #
################

usr <- Sys.info()[['user']]
if(usr == 'andrea')
{
    git.repository <-'~/git/longi_viral_loads'
    indir.deepsequence.data <- '~/Documents/Box/ratmann_pangea_deepsequencedata'
    indir.deepanalyses.xiaoyue <- '/home/andrea/HPC/project/ratmann_xiaoyue_jrssc2022_analyses/live'
    out.dir.prefix <- '/home/andrea/HPC/ab1820/home/projects/2022/longvl/'
}else if(usr == 'ab1820')
{
    git.repository <-'~/git/longi_viral_loads' 
    indir.deepsequence.data <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
    indir.deepanalyses.xiaoyue <- '/rds/general/project/ratmann_xiaoyue_jrssc2022_analyses/live/'
    out.dir.prefix <- '/rds/general/user/ab1820/home/projects/2022/longvl'
}

path.stan <- file.path(git.repository, 'stan')
path.tests <- file.path(indir.deepsequence.data, 
                        'RCCS_R15_R20',
                        "all_participants_hivstatus_vl_220729.csv")

path.census.eligible <- file.path(indir.deepsequence.data, 'RCCS_R15_R18/RCCS_census_eligible_individuals_220807.csv')
path.census.eligible.count <- file.path(indir.deepsequence.data, 'RCCS_R15_R18/RCCS_census_eligible_count_220719.csv')

file.exists(
        path.stan,
        path.tests,
        path.census.eligible
) |> all() |> stopifnot()

make_paper_numbers <- TRUE
if(make_paper_numbers)
    ppr_numbers <- list()

################
#   OPTIONS    #
################

option_list <- list(
    optparse::make_option(
        "--viremic-viral-load",
        type = "numeric",
        default = 200,
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
        "--indir",
        type = "character",
        default = '/home/andrea/HPC/ab1820/home/projects/2022/longvl/',
        help = "Path to output directory where model fits are stored [Defaults to my directory]", 
        dest = 'indir'
    ),
    optparse::make_option(
        "--outdir",
        type = "character",
        default = '/home/andrea/HPC/ab1820/home/projects/2022/longvl/',
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
source( file.path(git.repository,'functions/plotting_main_figures.R') )
source( file.path(git.repository,'scripts/phsc_vl_helpers.R'))

################
#     MAIN     #
################

VL_DETECTABLE = args$vl.detectable
VIREMIC_VIRAL_LOAD = args$viremic.viral.load

# plot requisite from nature med: 
naturemed_reqs() # stores them in nm_reqs

# output directory with rda files
out.dir <- args$out.dir.prefix
vl.out.dir <- file.path(out.dir, paste0('vl_', VIREMIC_VIRAL_LOAD) )

.f <- function(x) dir.create(file.path(vl.out.dir, x))
sapply(c('figures', 'tables'), .f)
stopifnot(dir.exists(vl.out.dir))

# get rda paths for fitted models
rda_files <- list.files(args$indir, pattern='.rda', full.names=T, recursive=TRUE)
rda_files <- grep(paste0('vl_', args$viremic.viral.load), rda_files, value=T)

# get data 
dall <- get.dall(path.tests)

# Get census eligible and compare with participants 
# __________________________________________________

# we need to chose whether to use the smoothed version, or the actual counts (`.count`)
dcens <- get.census.eligible()
dcens[, .(N_ELIGIBLE=sum(ELIGIBLE)), by=c('ROUND', 'SEX_LABEL')] |> kable()
last.round <- dcens[, max(ROUND)] 

cat('--- Make UNAIDS objectives plot ---\n')
tmp <- make.unaids.plots(DT=dcens)

# Summarised analyes
# __________________

cat('--- Plot Posteriors for fishing analyses ---\n')
p <- plot.all.gps(loc='fishing')
filename <- 'main_allanalyses_fishing.pdf'
ggsave2(p, file=filename, LALA=file.path(vl.out.dir, 'figures') , w=18, h=24, u='cm')

cat('--- Plot Posteriors for inland analyses ---\n')
p <- plot.all.gps(loc='inland')
filename <- 'main_allanalyses_inland.pdf'
ggsave2(p, file=filename, LALA=file.path(vl.out.dir, 'figures'), w=18, h=24, u='cm')

# Unaids table
# ____________

cat('--- Make UNAIDS goals table ---\n')
tmp <- make.table.unaids.goals()
print(xtable::xtable(tmp), 
      include.rownames=FALSE, 
      hline.after=c(-1, seq(0, nrow(tmp), nrow(tmp)/4)),
      file=file.path(vl.out.dir, 'tables', 'main_unaids_table.tex'))
xtable::xtable(tmp)

# Compare Prevalence of suppression among HIV+
# ____________________________________________

cat('-- Comparing prevalence of suppression among HIV+ ---\n')
tmp <- vl.vlprops.by.comm.gender.loc(dall)
vlc <- tmp$DT
tmp <- make.map.220810(dall, vlc, 'PVLNSofHIV_MEAN')


# Compare Rate of Change
# ______________________


compare.suppression <- function()
{
    .load.dre <- function(lab)
    {

        lab2 <- fcase(
                      lab %like% 'prevalence', "prev.hiv.by.age",
                      lab %like% 'suppAmongInfected', "nsinf.by.age",
                      lab %like% 'suppAmongPop', "nspop.by.age",
                      default=NA
        )
        stopifnot(! is.na(lab2))

        # dre contains posterior samples
        dre <- list()

        for (file in dfiles)
        {
            tmp_env <- new.env()
            load(file, envir=tmp_env)

            idx <- tmp_env$DT[, seq(min(AGE_LABEL), max(AGE_LABEL + 1), .5) ]
            idx <- data.table(Var2=seq_along(idx), AGE_LABEL=idx)

            cols <- grep('^p_predict', names(tmp_env$re), value=TRUE)
            tmp <- tmp_env$re[cols]
            .f <- function(x) as.data.table(reshape2::melt(x))
            tmp <- lapply(tmp, .f)

            .f <- function(x)
            {
                y <- merge(x, idx, by='Var2')
                y[, Var2 := NULL]
            }
            tmp <- lapply(tmp, .f)

            dre[[file]] <- tmp
        }
        rm(tmp_env)

        .f <- function(x) as.integer(gsub('^.*?round([0-9]{2}).*?$', '\\1', x))
        names(dre) <- .f(names(dre)) 

        dre
    }

    .label.sex.loc <- function(lst)
    {
        idx <- data.table(
                          .id=c('00', '01', '10', '11'),
                          SEX_LABEL=c('F', 'F', 'M', 'M'),
                          LOC_LABEL=c('inland', 'fishing', 'inland','fishing')
        )

        tmp <- rbindlist(lst, idcol=T)
        tmp[, .id := gsub('[a-z]|_','',.id)]
        tmp <- merge(tmp, idx, by='.id')
        tmp[, .id := NULL]
        tmp
    }

    rda_files <- list.files(args$indir, pattern='.rda', full.names=T)
    rda_files <- grep('16|19', rda_files, value=TRUE)

    # get posterior samples ratios of supp19/supp16
    dsupp_pop <- .load.dre('suppAmongPop')
    dsupp_pop <- lapply( dsupp_pop, .label.sex.loc)
    dsupp_pop <- rbindlist(dsupp_pop, idcol="ROUND")
    by_cols <- setdiff(names(dsupp_pop), c('value', 'ROUND'))
    setkey(dsupp_pop, ROUND)
    dsupp_pop <- dsupp_pop[, .(value=value[2]/value[1]) , by=by_cols]

    # now get quantiles
    .f <- function(x)
    {
        ps <- c(0.025,0.25,0.5,0.75,0.975)
        cols <- c('CL','IL','M','IU','CU')
        out <- as.list(quantile(x, probs=ps))
        names(out) <- cols
        out
    }

    by_cols <- setdiff(by_cols, "iterations")
    dsupp_pop <- dsupp_pop[, .f(value), by=by_cols ]

    dsupp_pop[, SEX:=fifelse(SEX_LABEL=='F', 'Women', 'Men')]

    p <- ggplot(data=dsupp_pop, aes(x=AGE_LABEL, color=LOC_LABEL, fill=LOC_LABEL)) +
        geom_line(aes(y=M)) + 
        geom_ribbon(aes(ymin=CL, ymax=CU), alpha=.5, color=NA) +
        geom_hline(aes(yintercept=1), linetype='dashed', color='darkred') + 
        facet_grid(~SEX) +
        scale_y_log10(expand=c(0,0)) +
        scale_x_continuous(expand=c(0,0)) +
        scale_fill_manual(values=palettes$comm) + 
        scale_color_manual(values=palettes$comm) + 
        theme(legend.position='bottom') +
        labs(color='Community type', fill='Community type',
             x='Age', 
             y='Ratio between prevalence of viremia\n at the end and start of study period\n(95% credible interval)')

    filename <- 'main_prevratio_r19r16_by_sex_loc.pdf'
    ggsave2(p, file=filename, w=18, h=12, u='cm')
}

if(0)
{
    # IT seems like the 'outliers' are legit...
    DT <- .preprocess.ds.oli(dall)
    vla1 <- .preprocess.make.vla(DT)
    vla1[ROUND == 19 & AGE_LABEL %between% c(30, 35) & LOC_LABEL == 'fishing' & SEX==1,
         .(AGE_LABEL, N, HIV_N, VLNS_N)]

    if(1)
    {
            tmp <- DT[, sort(unique(AGEYRS))]
            tmp1 <- DT[, sort(unique(ROUND))]
            tmp2 <- DT[, sort(unique(COMM_NUM))]
        vla <- as.data.table(expand.grid(ROUND=tmp1,
                                              FC=c('fishing','inland'),
                                              COMM_NUM=tmp2,
                                              SEX=c('M','F'),
                                              AGEYRS=tmp))
            vla <- vla[, {		
             z <- which(DT$ROUND==ROUND & DT$FC==FC & DT$SEX==SEX & DT$AGEYRS==AGEYRS & DT$COMM_NUM==COMM_NUM)	
             list(N          = length(z),
                  HIV_N      = sum(DT$HIV_STATUS[z]==1),
                  VLNS_N     = sum(DT$VLNS[z]==1),
                  ARV_N      = sum(DT$ARVMED[z]==0 & DT$HIV_STATUS[z]==1 & !is.na(DT$ARVMED[z])),
                  HIV_N      = sum(DT$HIV_STATUS[z]==1)
             )
            }, by=names(vla)]
     
        setnames(vla, c('FC','SEX','AGEYRS'), c('LOC_LABEL','SEX_LABEL','AGE_LABEL'))
        vla[, LOC:= as.integer(LOC_LABEL=='fishing')]
        vla[, SEX:= as.integer(SEX_LABEL=='M')]
        vla[, AGE:= AGE_LABEL-14L]
        vla[, ROW_ID:= seq_len(nrow(vla))]
     
             # Extract selected fields
             cols <- c('ROUND', 'LOC_LABEL', 'COMM_NUM', 'SEX_LABEL', 'AGE_LABEL', 'LOC', 'SEX', 'AGE', 'ROW_ID')
             tmp <- c('N', 'HIV_N', 'VLNS_N', 'ARV_N')
             cols <- unique(c(cols, tmp))
             vla <- vla[, ..cols]
    }
     
    tmp <- vla[ ROUND==19 & SEX_LABEL=="M" & LOC_LABEL == loc, ]
    cols <- c('M', 'CU', 'CL')
    tmp[, (cols) := binconf(VLNS_N, N, return.df=T), ]

    ggplot(data=tmp, aes(x=AGE_LABEL, y=M, color=as.factor(COMM_NUM))) +
            geom_point(position = position_dodge(width=.4)) + 
            geom_linerange(aes(ymin=CL, ymax=CU), position = position_dodge(width=.4)) + 
            facet_wrap(~ROUND)
    tmp
}
