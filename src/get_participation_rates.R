# TODO: change name of this script. 
# This should besome preliminary plots

library(data.table)
library(ggplot2)
kable <- knitr::kable

gitdir <- here::here()
source(file.path(gitdir, 'R/paths.R'))
source(file.path(gitdir.functions, 'phsc_vl_helpers.R'))
source(file.path(gitdir.functions, 'plotting_functions.R'))
source(file.path(gitdir.functions, 'paper_statements.R'))

# set output directory
outdir <- "/home/andrea/HPC/ab1820/home/projects/2022/longvl"
outdir.figures <- file.path(outdir, 'figures')
outdir.tables <- file.path(outdir, 'tables')

overwrite <- !interactive()
make_plots <- !interactive()
make_tables <- !interactive()

VL_DETECTABLE <- 150
VIREMIC_VIRAL_LOAD <- 1000 


# load stuff
cols <- c( "ROUND", "TYPE", "SEX", "AGEYRS", "ELIGIBLE")
ncen <- fread(path.census.eligible, select = cols)
setnames(ncen, c('TYPE'), c('FC'))

cols <- c("STUDY_ID", "ROUND", "SEX", "AGEYRS", "FC", 'HIV_STATUS', 'VL_COPIES', 'FIRST_PARTICIPATION')
dall <- get.dall(path.hivstatusvl.r1520) |> 
    subset(select=cols) |> unique()

# subset to rounds of interest
rounds_subsetted <- 16:19
ncen <- ncen[ ROUND %in% rounds_subsetted ]
dall <- dall[ ROUND %in% rounds_subsetted ]
ncen[, ROUND := as.integer(ROUND)]
dall[, ROUND := as.integer(ROUND)]

#################################################
catn("Make table with study pop characteristics")
#################################################

npar <- summarize.aggregates.dall()
dprop <- merge(npar, ncen, all.x=TRUE, all.y=TRUE)
stopifnot( any(is.na(dprop)) )

check_more_elig_than_part <- dprop[, all(N_PART < ELIGIBLE) ] 
stopifnot(check_more_elig_than_part)
tab_el <- paper_statements_average_participation(dprop)
tab <- make.table.eligible.participants(tab_el) 

# aggregate over age group
npar_agegroup <- subset(dprop, 
    AGEYRS %between% c(15, 49),
    select=c('ROUND', 'FC', 'SEX', 'AGEYRS', 'ELIGIBLE', 'N_PART', 'N_FIRST', 'N_HIV', 'N_HASVL')
)
npar_agegroup[, AGEGROUP := split.agegroup(AGEYRS, breaks=c(15, 25, 35, 50))]
npar_agegroup <- cube(npar_agegroup, 
    lapply(.SD, sum),
    .SDcols=names(npar_agegroup) %which.like% '^N|ELIGIBLE',
    by=c('ROUND', 'FC', 'SEX', 'AGEGROUP') ) |> 
    subset( ! is.na(ROUND) & !is.na(FC) )

if(make_tables){

    cols_lab <- c('ROUND_LAB', 'FC_LAB', 'SEX_LAB', 'AGEGROUP')
    tab <- copy(npar_agegroup)
    prettify_labels(tab)
    setcolorder(tab, cols_lab)
    setkeyv(tab, cols_lab)
    tab[, `:=` (ROUND=NULL, SEX=NULL, FC=NULL, AGEGROUP = as.character(AGEGROUP)) ]
    tab <- subset(tab, ! (is.na(SEX_LAB) & ! is.na(AGEGROUP))) 
    tab[, `:=` (SEX_LAB = fcoalesce(SEX_LAB, "Both"), AGEGROUP = fcoalesce(AGEGROUP, 'All')) ]
    tab <- delete.repeated.table.values(tab, cols = setdiff(cols_lab, 'AGEGROUP'))

    .cell <- function(NUM, DEN){
        prettify_cell( fmt_skeleton = "%s (%.1f%%)", format(NUM, big.mark=','), 100*NUM/DEN )
    }
    tab[, `:=` ( 
        N_PART = .cell(N_PART, ELIGIBLE),
        N_FIRST = .cell(N_FIRST, N_PART),
        N_HIV = .cell(N_HIV, N_PART),
        N_HASVL = .cell(N_HASVL, N_HIV)
    )]
    
    filename_tex <- file.path(outdir.tables, 'table_characteristics_participants.tex')
    write.to.tex(tab, file=filename_tex)
    filename <- 'table_characteristics_participants.pdf'
    p <- table.to.plot(tab)
    ggsave2(p=p, file=filename, LALA=outdir.tables, w=11, h=21)
}

#############
catn("Plots")
#############

catn("get participant proportion")
# ________________________________

# and aggregated participations rates rate 
key_cols <- c('ROUND', 'SEX', 'FC', 'AGEYRS' )
.ex <- expr(list(PART_RATE=round(100*sum(N_PART)/sum(ELIGIBLE),2)))
dprop[, eval(.ex), by=setdiff(key_cols, c('AGEYRS'))] |> kable()
dprop[, eval(.ex), by=setdiff(key_cols, c('AGEYRS', 'ROUND'))] |> kable()
dprop[, eval(.ex), by=setdiff(key_cols, c('AGEYRS', 'ROUND', 'FC'))] |> kable()

.ex <- expr(list(FIRST_RATE = round(100 * sum(N_FIRST)/sum(N_PART), 2)))
dprop[, eval(.ex), by=setdiff(key_cols, c('AGEYRS'))] |> kable()
dprop[, eval(.ex), by=setdiff(key_cols, c('AGEYRS', 'ROUND'))] |> kable()
dprop[, eval(.ex), by=setdiff(key_cols, c('AGEYRS', 'ROUND', 'FC'))] |> kable()

if(make_plots){
    # different age-pyramid plots
    p_pyramid_eligibleparticipants  <- plot.pyramid.bysexround( dprop, 
        .ylab = 'Number of participants among census-eligible individuals',
        NUM="N_PART",
        DEN='ELIGIBLE')
    filename <- 'pyramid_Neligible_Nparticipants.pdf'
    cmd <- ggsave2(p_pyramid_eligibleparticipants, file=filename, LALA=outdir.figures, w=10, h=11)

    p_pyramid_firstpart  <- plot.pyramid.bysexround( dprop, 
        .ylab = "Number of first-time participants among participants",
        NUM="N_FIRST",
        DEN='N_PART')
    filename <- 'pyramid_Nparticipants_Nfirst.pdf'
    cmd <- ggsave2(p_pyramid_firstpart , file=filename, LALA=outdir.figures, w=10, h=11)

    p_pyramid_hivp  <- plot.pyramid.bysexround( dprop, 
        .ylab = "Number of unsuppressed among HIV positive participants", 
        NUM="N_HIV",
        DEN='N_PART')
    filename <- 'pyramid_Nparticipants_Nhivp.pdf'
    cmd <- ggsave2(p_pyramid_hivp , file=filename, LALA=outdir.figures, w=10, h=11)


    p_pyramid_unsup  <- plot.pyramid.bysexround( dprop, 
        .ylab = "Number of unsuppressed among HIV positives", 
        NUM="N_VLNS",
        DEN='N_HIV')
    filename <- 'pyramid_Nhivp_Nunsup.pdf'
    cmd <- ggsave2(p_pyramid_unsup , file=filename, LALA=outdir.figures, w=10, h=11)

    p_pyramid_unsup_part <- plot.pyramid.bysexround( dprop, 
        .ylab = "Number of unsuppressed among participants", 
        NUM="N_VLNS",
        DEN="N_PART")
    filename <- 'pyramid_Nparticipants_Nunsup.pdf'
    cmd <- ggsave2(p_pyramid_unsup_part , file=filename, LALA=outdir.figures, w=10, h=11)
}

# get loess proportions

loess_ratepart <- dprop[, {
    ageyrs_topredict <- unique(AGEYRS) |> sort()
    .f <- function(x)
        loess(N_PART/ELIGIBLE ~ AGEYRS, span = x) |> predict(ageyrs_topredict)

    list(
        AGEYRS = ageyrs_topredict,
        N_PART = N_PART,
        ELIGIBLE = ELIGIBLE, 
        PARTRATE_RAW = N_PART / ELIGIBLE, 
        PARTRATE_SMOOTH.25 = .f(0.25),
        PARTRATE_SMOOTH.50 = .f(0.50),
        PARTRATE_SMOOTH.75 = .f(0.75)
    )
    
}, by=c('ROUND', 'FC', 'SEX')]

# .25 is wiggly, but so is raw data...
if(make_plots){
    p_partrate <- plot.smoothed.participation.rates(loess_ratepart)
    filename <- "smoothed_participationrates_bycommroundgenderage.pdf"
    cmd <- ggsave2(p=p_partrate, file=filename, LALA=outdir.figures, w=10, h=11 )
}

filename <- file.path( gitdir.data, "participation_rates_230517.rds" )
if(! file.exists(filename)){
    cat("Saving file", filename, "...\n")
    saveRDS(object=loess_ratepart, file=filename)
}else{
    cat("File", filename, "already exists...")
}


catn("what about the age composition/contribution of different pops")
# ___________________________________________________________________


# contribution to HIV
plot.agecontribution.fromN.stratby(dprop, 'N_HIV', agegroup=FALSE,
    .ylab="Contribution to PLHIV among participants") 

plot.agecontribution.fromN.stratby(dprop, 'N_HIV', agegroup = TRUE,
    .ylab="Contribution to PLHIV among participants");

# contribution to viraemia.
p_contrib_viraemia_parts <- plot.agecontribution.fromN.stratby(dprop, 
    'N_VLNS', 
    agegroup = TRUE,
    .ylab="Contribution to viraemia among participants") 
filename <- 'bars_contrib_viraemia_participants.pdf'
cmd <- ggsave2(p_contrib_viraemia_parts, file=filename, LALA=outdir.figures, w=10, h=11)


if (0) # Study ARVMED
{
    darv <- dall[HIV_STATUS == 1]

    cat("Assume NA ARV means no ARV usage...\n")
    tmp <- darv[HIV_AND_VL == 1]
    tmp[is.na(ARVMED) | ARVMED != 1, ARVMED := 0]
    by_cols <- c("ROUND", "FC", "SEX", "ARVMED")
    cols <- c("M", "CL", "CU")
    tmp[, Y := as.integer(VL_COPIES <= VIREMIC_VIRAL_LOAD)]
    tmp <- tmp[, binconf(sum(Y), .N, return.df = T), by = by_cols]
    names(tmp) <- c(by_cols, cols)
    setkeyv(tmp, by_cols)
    supp.prop.by.arv <- copy(tmp)

    ggplot(supp.prop.by.arv, aes(x = FC, colour = SEX)) +
        geom_point(aes(y = M, pch = as.factor(ARVMED)), position = position_dodge(width = .5)) +
        geom_linerange(aes(ymin = CL, ymax = CU, linetype = as.factor(ARVMED)), position = position_dodge(width = .5)) +
        facet_grid(~ROUND) +
        scale_y_continuous(labels = scales:::percent, limits = c(0, 1), expand = c(0, 0)) +
        scale_colour_manual(values = c("M" = "royalblue3", "F" = "deeppink2")) +
        theme(legend.position = "bottom") +
        labs(
            x = "Community type", y = "Proportion of suppressed measurements",
            linetype = "Ever reported ARV", pch = "Ever reported ARV",
            title = "Suppression by ARV reporting..."
        )


    cat("Suppression among participants reporting ARV \n")
    darv[ARVMED == 1 & is.na(VL_COPIES), STUDY_ID] -> idx
    tmp <- darv[STUDY_ID %in% idx, any(VL_COPIES == 0, na.rm = T), by = "STUDY_ID"]
    tmp[, cat(
        "Out of", .N, "HIV + participants with NA VL measurements", sum(V1),
        "were measured suppressed at least once\n"
    )]

    tmp <- darv[HIV_AND_VL == 1 & ARVMED == 1]
    by_cols <- c("ROUND", "FC", "SEX")
    cols <- c("M", "CL", "CU")
    tmp[, Y := as.integer(VL_COPIES <= VIREMIC_VIRAL_LOAD)]
    tmp <- tmp[, binconf(sum(Y), .N, return.df = T), by = by_cols]
    names(tmp) <- c(by_cols, cols)
    setkeyv(tmp, by_cols)
    supp.prop.among.report <- copy(tmp)

    # Sex plays a bigger role than community.
    ggplot(supp.prop.among.report, aes(x = FC, colour = SEX)) +
        geom_point(aes(y = M), position = position_dodge(width = .5)) +
        geom_linerange(aes(ymin = CL, ymax = CU), position = position_dodge(width = .5)) +
        facet_grid(~ROUND) +
        scale_y_continuous(labels = scales:::percent, limits = c(.5, 1), expand = c(0, 0)) +
        scale_colour_manual(values = c("M" = "royalblue3", "F" = "deeppink2")) +
        theme(legend.position = "bottom") +
        labs(x = "Community type", y = "Proportion of suppressed measurements", title = "Among people reporting ever ARV...")
}
