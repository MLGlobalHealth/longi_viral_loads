library(data.table)
library(ggplot2)
library(scales)
library(lubridate)
library(rstan)
library("haven")
library(here)

gitdir <- here::here()
source(file.path(gitdir, "R/paths.R"))

naturemed_reqs()

path.tests <- file.path(indir.deepdata.r1520, "all_participants_hivstatus_vl_220729.csv")

make.plots <- FALSE

# load files
community.keys <- fread(path.community.types)
quest <- fread(path.quest_r1519_221207) # NOT USE: path.quest_r1520
hiv <- fread(path.hivstatusvl.r1520)

###############################
# FIND SELF-REPORTED ART USE  #
###############################

# keep variable of interest
rin.cols <- c("ageyrs", "round", "study_id", "sex", "comm_num", "intdate", "arvmed", "cuarvmed")
rin <- subset(quest, select = rin.cols, subset = !(round %like% "R015"))
rin[, round := round2numeric(round)]


# find  community
community.keys[, comm := ifelse(strsplit(as.character(COMM_NUM_A), "")[[1]][1] == "f", "fishing", "inland"), by = "COMM_NUM_A"]
rinc <- merge(rin, community.keys, by.x = "comm_num", by.y = "COMM_NUM_RAW")

# to upper
colnames(rinc) <- toupper(colnames(rinc))

# find index of round
rinc <- rinc[order(STUDY_ID, ROUND)]
rinc[, INDEX_ROUND := 1:length(ROUND), by = "STUDY_ID"]

# restric age
rinc <- rinc[AGEYRS > 14 & AGEYRS < 50]

# get hiv status
rhiv <- hiv[ROUND %between% c(16, 19), .(STUDY_ID, ROUND, HIV = HIV_STATUS)]
hivs <- merge(rhiv, rinc, by = c("STUDY_ID", "ROUND"), all.x = TRUE, all.y = TRUE)

# get vl count

# study missing information
cube(hivs,
    by = "ROUND",
    j = list(
        NO_RINC = sum(is.na(COMM)),
        NO_RHIV = sum(is.na(HIV)),
        TOT = uniqueN(STUDY_ID)
    )
) |>
    knitr::kable(caption = "Number of participants in testing who do not appear in questionnaire")

# keep HIV positive
rprev <- hivs[HIV == 1]

# get ART status

# I would do it as follows
rprev[, `:=`(
    EVER_ART = (ARVMED == 1 | CUARVMED == 1),
    CURRENT_ART = CUARVMED == 1
)]

# check EVER_ART is monotonically increasing.
# ie: no one reports never ART use after reporting it (2% do)
setkey(rprev, STUDY_ID, ROUND)
idx <- rprev[,
    {
        stopifnot(!is.unsorted(ROUND))
        all(EVER_ART == cummax(EVER_ART))
    },
    by = STUDY_ID
][V1 == FALSE, STUDY_ID]
rprev[, table(STUDY_ID %in% idx)] |> proportions()


########################
# ADD VIRAL LOAD DATA  #
########################

# for round with suppressed set art to true if indiv is suppressed

# tuning
VL_DETECTABLE <- 400
VIREMIC_VIRAL_LOAD <- 1000 # WHO standards

source(file.path(gitdir.functions, "phsc_vl_helpers.R"))

dall <- get.dall(path.hivstatusvl.r1520, make_flowchart = FALSE) |>
    subset(AGEYRS %between% c(15, 49))
# dall[, uniqueN(STUDY_ID)]


# keep infected
dall <- dall[HIV_STATUS == 1]

# merge to self-reported data
tmp <- dall[, .(STUDY_ID, ROUND, SEX, FC, VLNS = VL_COPIES >= VIREMIC_VIRAL_LOAD, AGEYRS)]

# tmp[, ROUND := paste0("R0", ROUND)]
setnames(tmp, "AGEYRS", "AGEYRS2")
rprev <- merge(rprev, tmp, by.x = c("STUDY_ID", "ROUND", "SEX", "COMM"), by.y = c("STUDY_ID", "ROUND", "SEX", "FC"), all.x = T, all.y = T)

# set ageyrs to the viral load data if available
rprev[!is.na(AGEYRS2), AGEYRS := AGEYRS2]
set(rprev, NULL, "AGEYRS2", NULL)

# set art to true if suppressed viral load
rprev[VLNS == 0, EVER_ART := TRUE]
rprev[, sprintf("Removing %d individuals with NA EVER_ART", sum(is.na(EVER_ART)))]
rprev <- rprev[!is.na(EVER_ART)]

# remove na vlns for round 16 onwards otherwise it leads to % art < % suppressed
rprev[, .(`NO Viral load test` = sum(is.na(VLNS)) / sum(!is.na(VLNS))), by = "ROUND"] |>
    knitr::kable(caption = "Prop. of NA VL test among HIV participants")

rprev <- rprev[!is.na(VLNS)]

#################################

# KEEP INDIVIDUALS SEEN FOR THE FIRST TIME
# THAT ARE THE CLOSEST TO NON-PARTICIPANTS

#################################

# get participants in R09 -> R014
study_ids_in_rounds_0914 <- fread(path.hivres.r0914, select = "study_id") |>
    unlist() |>
    unique()
# get participants in R15 (btw 60% new in R15?)
study_ids_in_rounds_15 <- quest[round %like% 15, unique(study_id)]
# Merge rounds 09-14 to 15
study_ids_in_rounds_0915 <- c(
    study_ids_in_rounds_0914,
    study_ids_in_rounds_15
) |> unique()

get.first.participants.round.16.19 <- function(DT) {
    # count study_id which appear first
    dfirst <- DT |> subset(ROUND %between% c(16, 19))

    # define first participants
    dfirst[STUDY_ID %in% study_ids_in_rounds_0915, FIRST_VISIT := FALSE]
    dfirst[is.na(FIRST_VISIT), FIRST_VISIT := {
        stopifnot(!is.unsorted(ROUND))
        c(TRUE, rep(FALSE, .N - 1))
    }, by = "STUDY_ID"]

    cols <- c("FIRST_VISIT", "AGEYRS", "STUDY_ID", "SEX", "ROUND")
    dfirst |> subset(select = cols)
}

first_participants <- get.first.participants.round.16.19(hiv)
# first_participants[ AGEYRS %between% c(15, 49), comma(.N), by='FIRST_VISIT']
# first_participants[ AGEYRS %between% c(15, 49) & ROUND %in% c(16, 19),
#     any(FIRST_VISIT), by="STUDY_ID"][, comma(table(V1))]
#
# # dall[, any(FIRST_PARTICIPATION == 1), by="STUDY_ID"][, comma(table(V1))]

rprev <- merge(x = rprev, y = first_participants[, .(STUDY_ID, ROUND, FIRST_VISIT)], all.x = TRUE, by = c("STUDY_ID", "ROUND"))
stopifnot(!(rprev$FIRST_VISIT |> is.na() |> any()))
# first participants seem to be ~ 21% of the participants.
cube(rprev, mean(FIRST_VISIT == TRUE), by = "ROUND") |>
    knitr::kable(caption = "Proportion of new participants by round")


# plot
if (make.plots) {
    # NOTE: not accounting for Community types, easily fixable

    plot.first.participant.histogram <- function(DT, y_lab = "TODO") {
        by_cols <- c("FIRST_VISIT", "AGEYRS", "SEX", "ROUND")
        # get counts
        DT[!is.na(AGEYRS) & AGEYRS %between% c(15, 49),
            any(FIRST_VISIT),
            by = STUDY_ID
        ][, .(comma(.N), comma(sum(V1)))]
        DT[, FIRST_VISIT, ]
        browser()
        dcounts <- DT[, .N, by = by_cols]
        dcounts <- merge(dcounts,
            dcounts[, CJ(
                AGEYRS = unique(AGEYRS),
                SEX = unique(SEX),
                ROUND = unique(ROUND)
            )],
            by = c("AGEYRS", "SEX", "ROUND"), all.y = TRUE
        )
        dcounts[is.na(N), N := 0]

        dcounts |>
            prettify_labels() |>
            subset(AGEYRS %between% c(15, 50)) |>
            ggplot(aes(x = AGEYRS, fill = SEX_LAB, alpha = as.integer(FIRST_VISIT), y = N)) +
            geom_col(position = position_stack(reverse = FALSE), aes(group = interaction(SEX_LAB, FIRST_VISIT)), color = "black") +
            facet_grid(ROUND_LAB ~ SEX_LAB) +
            theme_default() +
            scale_fill_manual(values = palettes$sex, labels = sex_dictionary2) +
            scale_y_expand_lower0 +
            scale_alpha_continuous(breaks = c(0, 1)) +
            guides(alpha = guide_legend(order = 1)) +
            my_labs(y = y_lab, alpha = "First participants") +
            theme_default()
    }

    p1 <- first_participants |>
        plot.first.participant.histogram(y_lab = "Number of first-time participants among all participants")
    filename <- "hist_propfirstparticipants_bysexround.pdf"
    ggsave2(p = p1, file = filename, LALA = file.path(OUTDIR, "figures"), w = 18, h = 18, u = "cm")

    # now do analogous same among positives, taking FIRST_VISIT classification from first_participants
    p2 <- rprev |> plot.first.participant.histogram(y_lab = "Number of first-time participants among HIV positive")
    filename <- "hist_propfirstparticipants_amongpositives_bysexround.pdf"
    ggsave2(p = p2, file = filename, LALA = file.path(OUTDIR, "figures"), w = 18, h = 18, u = "cm")

    # what about non suppressed
    p_viraemic <- rprev |>
        subset(VLNS == TRUE) |>
        plot.first.participant.histogram(y_lab = "Number of first-time participants among viraemic")
    filename <- "hist_propfirstparticipants_amongviraemic_bysexround.pdf"
    ggsave2(p = p_viraemic, file = filename, LALA = file.path(OUTDIR, "figures"), w = 18, h = 18, u = "cm")
    # can I visualise both on the same plot?
    # NOTE: assumning those not appearing in dplot are negative
    # ggpubr::ggarrange(p1, p2, ncol=2, common.legend = TRUE, legend='bottom')
    dplot_all <- merge(
        first_participants,
        dplot[, .(STUDY_ID, ROUND, HIV)],
        by = c("STUDY_ID", "ROUND"), all.x = TRUE
    )[is.na(HIV), HIV := 0]

    # find counts
    by_cols <- c("ROUND", "SEX", "HIV", "FIRST_VISIT", "AGEYRS")
    dcounts <- dplot_all[, .N, by = by_cols]
    dcounts[, CJ(
        ROUND = unique(ROUND),
        SEX = unique(SEX),
        AGEYRS = unique(AGEYRS),
        HIV = unique(HIV),
        FIRST_VISIT = unique(FIRST_VISIT)
    )] |>
        merge(x = dcounts, y = _, by = by_cols, all.y = TRUE)
    dcounts[is.na(N), N := 0]

    dcounts[, .N, by = c("ROUND", "SEX", "HIV", "FIRST_VISIT", "AGEYRS")][, all(N == 1)]
    dcounts[, HIV := as.logical(HIV)]
    dcounts[, P := N / sum(N), by = c("ROUND", "SEX", "FIRST_VISIT", "AGEYRS")]
    dcounts <- dcounts |>
        subset(AGEYRS %between% c(15, 50)) |>
        prettify_labels()

    p3 <- dcounts |>
        ggplot(aes(group = interaction(HIV_LAB, FIRST_VISIT), fill = HIV_LAB, alpha = FIRST_VISIT, x = AGEYRS, y = N)) +
        geom_col(position = "stack", color = "black") +
        facet_grid(ROUND_LAB ~ SEX_LAB) +
        theme_default() +
        scale_alpha_manual(values = c(.4, 1)) +
        scale_fill_manual(values = palettes$hivstatus) +
        scale_y_expand_lower0 +
        guides(alpha = guide_legend(order = 1)) +
        my_labs(y = "Number of RCCS participants", alpha = "First participants") +
        theme_default()
    filename <- "hist_Nfirstparticipantsandpositives_bysexround.pdf"
    ggsave2(p = p3, file = filename, LALA = file.path(OUTDIR, "figures"), w = 18, h = 18, u = "cm")

    # NOTE: this may be better with points and linerange (AC binomial intervals)
    p4 <- dcounts |>
        ggplot(aes(group = interaction(HIV_LAB, FIRST_VISIT), fill = HIV_LAB, alpha = FIRST_VISIT, x = AGEYRS, y = P)) +
        geom_col(data = subset(dcounts, FIRST_VISIT == 1), aes(x = AGEYRS - .25), position = "stack", color = "black", width = .5) +
        geom_col(data = subset(dcounts, FIRST_VISIT == 0), aes(x = AGEYRS + .25), position = "stack", color = "black", width = .5) +
        facet_grid(ROUND_LAB ~ SEX_LAB) +
        theme_default() +
        scale_y_percentage +
        scale_alpha_manual(values = c(.4, 1)) +
        scale_fill_manual(values = palettes$hivstatus) +
        guides(alpha = guide_legend(order = 1)) +
        my_labs(y = "Number of RCCS participants", alpha = "First participants") +
        theme_default()
    filename <- "hist_Pfirstparticipantsandpositives_bysexround.pdf"
    ggsave2(p = p4, file = filename, LALA = file.path(OUTDIR, "figures"), w = 18, h = 18, u = "cm")


    # compare first participants vs other
    plot.comparison.firstparticipants.andnon <- function(DT) {
        DT[, AGEGROUP := split.agegroup(AGEYRS)]

        DT[!is.na(HIV),
            {
                z <- Hmisc::binconf(x = sum(VLNS), n = .N, return.df = TRUE)
                names(z) <- c("M", "CL", "CU")
                z
            },
            by = c("ROUND", "COMM", "SEX", "AGEGROUP", "FIRST_VISIT")
        ] |>
            prettify_labels() |>
            ggplot(aes(x = AGEGROUP, y = M, ymin = CL, ymax = CU, color = SEX_LAB, pch = FIRST_VISIT, linetype = FIRST_VISIT)) +
            geom_point(position = position_dodge(width = .4)) +
            # geom_line(position=position_dodge(width=.4)) +
            geom_linerange(position = position_dodge(width = .4)) +
            facet_grid(ROUND_LAB ~ SEX_LAB * COMM) +
            scale_y_percentage +
            scale_color_manual(values = palettes$sex, labels = sex_dictionary2) +
            theme_default() +
            theme(axis.text.x = element_text(angle = 45, vjust = .5, hjust = .5)) +
            my_labs(y = "Prevalence of suppression among HIV positive") +
            p_reqs
    }
    p_comp_suppamonghiv <- plot.comparison.firstparticipants.andnon(rprev)
    filename <- "points_Psuppamonghiv_byfirstpart_sexround.pdf"
    ggsave2(p = p_comp_suppamonghiv, file = filename, LALA = file.path(OUTDIR, "figures"), w = 22, h = 18, u = "cm")
}



##########################################
# aggregate by round, sex, comm and age  #
##########################################

# KEEP INDIVIDUALS SEEN FOR THE FIRST TIME
# find self reported under art for participant
rart <- rprev[FIRST_VISIT == TRUE,
    list(
        COUNT = sum(EVER_ART == T),
        TOTAL_COUNT = length(EVER_ART)
    ),
    by = c("ROUND", "SEX", "COMM", "AGEYRS")
]

############################
# save de-identified data  #
############################

# is ART coverage what we want? Or may we want something else?
file.name <- file.path(gitdir.data, "aggregated_newlyregistered_count_art_coverage.csv")
if (!file.exists(file.name)) {
    cat("\n Saving", file.name, "...\n")
    fwrite(rart, file = file.name, row.names = F)
} else {
    cat("\n Output file", file.name, "already exists.\n")
}



###############################
# get suppression proportions #
###############################

rprev[, table(HIV)]

cols <- c("STUDY_ID", "ROUND", "HIV_STATUS", "HIV_VL")

fread(path.hivstatusvl.r1520)
