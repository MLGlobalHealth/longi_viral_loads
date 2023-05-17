# TODO: change name of this script. 
# This should besome preliminary plots

library(data.table)
library(ggplot2)

gitdir <- here::here()
source(file.path(gitdir.R, 'paths.R'))
source(file.path(gitdir.functions, 'phsc_vl_helpers.R'))

# set output directory
outdir <- "/home/andrea/HPC/ab1820/home/projects/2022/longvl"
outdir.figures <- file.path(outdir, 'figures')
outdir.tables <- file.path(outdir, 'tables')

VL_DETECTABLE <- 150
VIREMIC_VIRAL_LOAD <- 1000 

# load stuff
cols <- c( "ROUND", "TYPE", "SEX", "AGEYRS", "ELIGIBLE")
ncen <- fread(path.censusN, select = cols)
cols <- c("STUDY_ID", "ROUND", "SEX", "AGEYRS", "FC", 'HIV_STATUS', 'VL_COPIES', 'FIRST_PARTICIPATION')
dall <- get.dall(path.hivstatusvl.r1520) |> 
    subset(select=cols) |> unique()

# attribute some NA first parts
idx <- dall[is.na(FIRST_PARTICIPATION), STUDY_ID]
if(length(idx) > 0)
{
    sprintf('%s participants have no first participation status. Imputing...', 
        length(idx)) |> warning()
    # impute non first if other participations, else assign as first
    idx_notfirst <- dall[STUDY_ID %in% idx, any(FIRST_PARTICIPATION == 1), by='STUDY_ID'][
        V1==TRUE, STUDY_ID]
    dall[STUDY_ID %in% idx_notfirst & is.na(FIRST_PARTICIPATION), FIRST_PARTICIPATION := 0]
    dall[ is.na(FIRST_PARTICIPATION), FIRST_PARTICIPATION := 0]
}

# make compatible
setnames(ncen, c('TYPE'), c('FC'))

# subset to rounds of interest
rounds_subsetted <- 16:19
ncen <- ncen[ ROUND %in% rounds_subsetted ]
dall <- dall[ ROUND %in% rounds_subsetted ]
ncen[, ROUND := as.integer(ROUND)]
dall[, ROUND := as.integer(ROUND)]

# get participant proportion
# __________________________

key_cols <- c('ROUND', 'SEX', 'FC', 'AGEYRS' )
npar <- dall[, .(
    N_PART=.N,
    N_FIRST=sum(FIRST_PARTICIPATION),
    N_HASVL=sum(!is.na(VL_COPIES) & HIV_STATUS == 1),
    N_HIV=sum(HIV_STATUS), 
    N_VLNS=sum( VL_COPIES >= VIREMIC_VIRAL_LOAD ,na.rm=TRUE),
    N_LOWVL=sum( VL_COPIES >= 200)
    ), by=key_cols]
dprop <- merge(npar, ncen)
check_more_elig_than_part <- dprop[, all(N_PART < ELIGIBLE) ] 
stopifnot(check_more_elig_than_part)

# and aggregated participations rates rate 
.ex <- expr(list(PART_RATE=round(100*sum(N_PART)/sum(ELIGIBLE),2)))
.ex <- expr(list(FIRST_RATE = round(100 * sum(N_FIRST)/sum(N_PART), 2)))
dprop[, eval(.ex), by=]
dprop[, eval(.ex), by=setdiff(key_cols, c('AGEYRS'))]
dprop[, eval(.ex), by=setdiff(key_cols, c('AGEYRS', 'ROUND'))]
dprop[, eval(.ex), by=setdiff(key_cols, c('AGEYRS', 'ROUND', 'FC'))]



# different age-pyramid plots
p_pyramid_eligibleparticipants  <- plot.pyramid.bysexround( dprop, 
    .ylab = 'Number of participants among census eligible individuals',
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


# plot.pyramid.bysexround( dprop, 
#     .ylab = "Number of unsuppressed among participants", 
#     NUM="N_HASVL",
#     DEN="N_HIV")


# what about the age composition/contribution of different pops
# _____________________________________________________________


# contribution to HIV
plot.agecontribution.fromN.stratby(dprop, 'N_HIV', agegroup=FALSE,
    .ylab="Contribution to PLHIV among participants") 

plot.agecontribution.fromN.stratby(dprop, 'N_HIV', agegroup = TRUE,
    .ylab="Contribution to PLHIV among participants") 

# contribution to viraemia.
p_contrib_viraemia_parts <- plot.agecontribution.fromN.stratby(dprop, 
    'N_VLNS', 
    agegroup = TRUE,
    .ylab="Contribution to viraemia among participants") 
filename <- 'bars_contrib_viraemia_participants.pdf'
cmd <- ggsave2(p_contrib_viraemia_parts, file=filename, LALA=outdir.figures, w=10, h=11)
cmd


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
