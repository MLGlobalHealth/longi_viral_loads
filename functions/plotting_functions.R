plot_visitpairs <- function(DT) {
    ggplot(DT) +
        geom_segment(aes(x = hivdate, xend = hivdate_2, y = hiv_vl, yend = hiv_vl_2, colour = type)) +
        geom_hline(yintercept = log10(THRESHOLD + 0.01)) +
        labs(x = "Visit date", y = "log10(viral load)") +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            legend.position = "bottom"
        )
}

plot_rounds_vl_collection <- function(DT, rounds = 16:19, atleast = NULL, community = NULL, filename = NULL) 
{
    cols <- c("study_id", "comm", "round")
    dcounts <- DT[round %in% rounds & !is.na(hiv_vl), .SD, .SDcols = cols]

    subtitle_comm <- ""
    if (!is.null(community)) {
        subtitle_comm <- paste(" in", community, "communitites")
        dcounts <- participant_comm2mode(dcounts)
        dcounts <- dcounts[comm %in% community]
    }

    subtitle_atleast <- ""
    if (!is.null(atleast)) { # extract individuals with sufficient VL measurements
        subtitle_atleast <- paste(" with at least", atleast, "VL measurements")
        tmp <- dcounts[, .N, by = "study_id"]
        idx <- tmp[N >= atleast, unique(study_id)]
        dcounts <- dcounts[study_id %in% idx]
    }

    nms <- dcounts[, sort(unique(round))]
    tmp1 <- dcounts[, list(round = nms), by = c("study_id", "comm")]
    tmp1 <- fsetdiff(tmp1, dcounts)
    dcounts[, measurement := TRUE]
    tmp1[, measurement := FALSE]
    dcounts <- rbind(dcounts, tmp1)
    setkey(dcounts, study_id, round)
    dcounts[, round := round2factor(round)]

    p1 <- ggplot(data = dcounts, aes(x = as.factor(round), fill = measurement)) +
        geom_bar(stat = "count") +
        theme_bw() +
        theme(
            axis.title = element_text(),
            legend.position = "bottom"
        ) +
        labs(
            x = "Rounds", y = "Number of participants", fill = "VL collected?",
            title = "In which rounds were participants VL collected?",
            subtitle = paste0(
                "(",
                dcounts[, uniqueN(study_id), ], " participants",
                subtitle_comm,
                subtitle_atleast, ")"
            )
        )

    if (!is.null(filename)) ggsave2(filename, p1, 7, 5)
    p1
}


plot_classification <- function(DT, DCLASS, filename_gif = NULL, filename = NULL) 
{
    tmp <- merge(DT, DCLASS, by = "study_id")
    tmp[, hiv_vl2 := log10(hiv_vl + 0.02 + runif(.N, -0.015, 0.015))]

    require(ggthemes)

    tmp[, class := .dict.class(class)]
    p1 <- ggplot(tmp, aes(x = hivdate, y = hiv_vl2, group = study_id, colour = class)) +
        geom_line(size = .3) +
        geom_hline(
            aes(yintercept = log10(THRESHOLD + .02)),
            linetype = "dotted",
            size = 2
        ) +
        theme_fivethirtyeight() +
        theme(
            axis.title = element_text(),
            text = element_text(family = "Rubik"),
            legend.text = element_text(size = 10)
        ) +
        labs(
            title = "VL Trajectory classification",
            x = "Visit date",
            y = "Viral Load (log base 10)",
            color = "Class"
        ) +
        guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
        scale_colour_tableau()


    if (is.null(filename_gif)) {
        filename_gif <- file.path(out.dir, "viral_trajectories_classes.gif")
    }

    if (usr == "andrea" & !file.exists(filename_gif) & dir.exists(dirname(filename_gif))) {
        require(gganimate)
        p1.animated <- p1 +
            transition_states(class) +
            shadow_mark(alpha = .1) +
            labs(subtitle = "Class: {previous_state}")

        p1.animated <- animate(p1.animated, height = 800, width = 1200)
        anim_save(filename_gif, p1.animated)
    }

    p1
}

plot_classes_by_sex_age <- function(DCLASS, DVL, filename = NULL) {
    # subtitle <- paste0('(Measurements from round, ', DVL[, min(round)], ' onwards)')
    setkey(DVL, study_id, round)
    tmp <- DVL[, hivdate[.N], by = "study_id"]
    tmp <- merge(tmp, dbirth, by = "study_id", all.x = TRUE)
    tmp[, age_last_measurement := round((V1 - birthdate) / 365.25)]

    tmp1 <- DCLASS[, .(study_id, class)]
    tmp1 <- merge(tmp1, tmp, all.x = TRUE)
    tmp2 <- tmp1[, .N, by = c("class", "sex")]
    tmp2[, `:=`(
        xpos = dplyr::if_else(sex == "F", Inf, Inf),
        ypos = dplyr::if_else(sex == "F", -Inf, Inf),
        hjustvar = dplyr::if_else(sex == "F", 1.5, 1.5),
        vjustvar = dplyr::if_else(sex == "F", 1.5, 1.5),
        lab = paste0(sex, ": N=", N)
    )]

    p <- ggplot(
        data = tmp1,
        aes(x = age_last_measurement, fill = sex, y = ifelse(sex == "M", yes = 1, no = -1))
    ) +
        geom_bar(stat = "identity") +
        geom_hline(yintercept = 0) +
        geom_text(
            data = tmp2[sex == "M"],
            aes(x = xpos, y = ypos, hjust = 1.5, vjust = 1.5, label = lab)
        ) +
        geom_text(
            data = tmp2[sex == "F"],
            aes(x = xpos, y = ypos, hjust = -.5, vjust = 1.5, label = lab)
        ) +
        scale_y_continuous() +
        coord_flip() +
        facet_wrap(~ .dict.class(class), scales = "free_x") +
        theme_bw() +
        scale_fill_manual(values = palettes$sex) +
        theme(legend.position = "bottom") +
        labs(
            y = "count", x = "age at last measurement",
            title = "VL Trajectories classifications: age and sex composition"
        )

    if (!is.null(filename)) ggplot2(filename, p, w = 10, h = 8)
    p
}

plot_classes_by_comm_age <- function(DCLASS, DVL, filename = NULL) {
    # subtitle <- paste0('(Measurements from round, ', DVL[, min(round)], ' onwards)')
    # DVL <- copy(dvl_15)
    # DCLASS <- copy(dclass_15)
    tmp <- participant_comm2mode(DVL, quietly = TRUE)
    stopifnot(tmp[, uniqueN(comm), by = "study_id"][, all(V1)])

    setkey(tmp, study_id, round)
    tmp <- tmp[, hivdate[.N], by = c("study_id", "comm")]
    tmp <- merge(tmp, dbirth, by = "study_id", all.x = TRUE)
    tmp[, age_last_measurement := round((V1 - birthdate) / 365.25)]

    tmp1 <- DCLASS[, .(study_id, class)]
    tmp1 <- merge(tmp1, tmp, all.x = TRUE)
    tmp2 <- tmp1[, .N, by = c("class", "comm")]
    tmp2[, `:=`(
        xpos = dplyr::if_else(comm == "fishing", Inf, Inf),
        ypos = dplyr::if_else(comm == "fishing", -Inf, Inf),
        hjustvar = dplyr::if_else(comm == "fishing", 1.5, 1.5),
        vjustvar = dplyr::if_else(comm == "fishing", 1.5, 1.5),
        lab = paste0(comm, ": N=", N)
    )]

    p <- ggplot(
        data = tmp1,
        aes(
            x = age_last_measurement, fill = comm,
            y = ifelse(comm == "inland", yes = 1, no = -1)
        )
    ) +
        geom_bar(stat = "identity") +
        geom_hline(yintercept = 0) +
        geom_text(
            data = tmp2[comm == "inland"],
            aes(x = xpos, y = ypos, hjust = 1.5, vjust = 1.5, label = lab)
        ) +
        geom_text(
            data = tmp2[comm == "fishing"],
            aes(x = xpos, y = ypos, hjust = -.5, vjust = 1.5, label = lab)
        ) +
        scale_y_continuous() +
        coord_flip() +
        facet_wrap(~ .dict.class(class), scales = "free_x") +
        theme_bw() +
        scale_fill_manual(values = palettes$comm) +
        theme(legend.position = "bottom") +
        labs(
            y = "count", x = "age at last measurement", fill = "community type",
            title = "VL Trajectories classifications: age and community composition"
        )

    if (!is.null(filename)) ggplot2(filename, p, w = 10, h = 8)
    p
}

plot_scatter_cd4vl_bygroup <- function(DCD4 = dcd4, DVL = dvl, group, filename = NULL) {
    # group <- 'arvmed'
    # DCD4 <- copy(dcd4)
    # DVL <- copy(tmp)
    dmerged <- merge(DCD4, DVL, by = c("study_id", "round"))
    dmerged <- dmerged[!is.na(cd4) & !is.na(hiv_vl), ]
    # dmerged[, plot(hist((hivdate - cd4date)/365))]

    .f <- function(n) {
        dmerged[, table(abs(hivdate - cd4date) > n)]
    }
    lapply(c(1, 10, 100, 365, 365 * 2, 365 * 3), .f)

    setnames(dmerged, group, "group")
    txt <- dmerged[hiv_vl != 0, .N, by = "group"][, lab := paste0(group, ": ", N, "\n")]
    txt <- txt[, list(
        lab = paste0(txt$lab, collapse = ""),
        xpos = Inf, ypos = Inf
    )]
    txt

    # text <- dmerged[hiv_vl != 0, .N, by='comm'][, text:=paste0(comm,': ', N)]

    rt3 <- function(x) {
        x^(1 / 3)
    }
    p <- ggplot(data = dmerged, aes(x = rt3(cd4), y = log10(hiv_vl), color = group, fill = group)) +
        geom_point() +
        geom_smooth(data = dmerged[hiv_vl != 0], method = lm) +
        geom_text(
            data = txt,
            aes(x = xpos, y = ypos, hjust = 1.5, vjust = 1.5, label = lab),
            inherit.aes = FALSE
        ) +
        theme_bw() +
        scale_color_manual(values = `[[`(palettes, group)) +
        scale_fill_manual(values = `[[`(palettes, group)) +
        theme(
            legend.position = "bottom",
            # axis.title.x = element_blank(),
            # axis.title.y = element_blank(),
            # axis.text.x = element_blank(),
            # axis.text.y = element_blank(),
            plot.margin = unit(c(3, -5.5, 4, 3), "mm")
        ) +
        labs(x = "third rooted cd4 count(copies/mL)", y = "log10 viral load (copies/mL)", fill = group, color = group)

    dens_x <- ggplot(dmerged, aes(x = rt3(cd4), fill = group)) +
        geom_histogram(alpha = .8, position = "identity") +
        theme_void() +
        scale_fill_manual(values = `[[`(palettes, group)) +
        theme(legend.position = "none")

    dens_y <- ggplot(dmerged, aes(x = log10(hiv_vl), fill = group)) +
        geom_histogram(alpha = .8, position = "identity") +
        theme_void() +
        coord_flip() +
        scale_fill_manual(values = `[[`(palettes, group)) +
        theme(legend.position = "none")

    require(patchwork)
    p1 <- dens_x + plot_spacer() + p + dens_y +
        plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
    p1

    if (!is.null(filename)) ggplot2(filename, p1, 8, 8)
    p1
}


plot.n.census.eligible.smooth <- function(DT) {
    tmp <- subset(DT, !ROUND %like% "15",
        select = c("TYPE", "SEX", "AGEYRS", "ROUND", "ELIGIBLE", "ELIGIBLE_SMOOTH.25", "ELIGIBLE_SMOOTH.50", "ELIGIBLE_SMOOTH.75")
    ) |>
        melt.data.table(
            id.vars = c("TYPE", "SEX", "AGEYRS", "ELIGIBLE", "ROUND"),
            variable.name = "SMOOTH_PAR"
        ) |>
        prettify_labels()

    tmp[, SMOOTH_PAR := gsub("^.*\\.([0-9]+)", "\\1", SMOOTH_PAR)]

    df_label <- tmp[, .(
        diff = round(abs(sum(ELIGIBLE) - sum(value))),
        ylevel = 900 - 3 * as.integer(SMOOTH_PAR)
    ),
    by = c("TYPE", "SEX", "SMOOTH_PAR", "ROUND")
    ]

    p <- ggplot(tmp, aes(x = AGEYRS)) +
        geom_bar(data = unique(tmp[, .(TYPE, SEX_LAB, AGEYRS, ELIGIBLE, ROUND_LAB)]), aes(y = ELIGIBLE), stat = "identity", alpha = 0.5) +
        geom_line(aes(y = value, col = SMOOTH_PAR)) +
        labs(y = "Census eligible individuals", x = "Age", color = "Loess parameter") +
        facet_grid(ROUND_LAB ~ TYPE + SEX_LAB) +
        theme_default() +
        geom_label(data = df_label, aes(x = 49, y = ylevel, label = diff, col = SMOOTH_PAR), size = 3, label.size = NA)
    p
}

plot.proportion.firstparticipants.by.vars <- function(DT, by_cols = c("ROUND", "TYPE")) {
    # DT <- copy(dfirstround); by_cols = c('ROUND','TYPE')
    tmp <- cube(dfirstround,
        j = list(
            N_FIRST = sum(ROUND == QST_1STRND),
            N_REPEATED = sum(ROUND != QST_1STRND),
            N_TOT = .N
        ), by = by_cols
    )
    tmp[, (by_cols) := lapply(.SD, function(x) {
        x[is.na(x)] <- "Total"
        x
    }), .SDcols = by_cols]
    tmp[, P_FIRST := N_FIRST / N_TOT]

    p_prop_newparts <- tmp |>
        subset(ROUND != "Total") |>
        prettify_labels() |>
        ggplot(aes(x = as.integer(ROUND), y = P_FIRST, color = TYPE)) +
        geom_line() +
        scale_y_continuous(
            labels = scales::label_percent(),
            limits = c(0, .4),
            expand = expansion(mult = 0)
        ) +
        theme_default() +
        my_labs(x = "Interview round", y = "Proportion of new participants", color = "Community type")

    p_prop_newparts
}

plot.pyramid.eligible.participants <- function(DT) {
    cols <- c("N_PART", "ELIGIBLE")
    new_cols <- paste(cols, "PYR", sep = "_")

    DT[, (new_cols) := lapply(.SD, function(x) (-1 + 2 * (SEX == "F")) * x), .SDcols = cols]

    dlabs <- DT[,
        .(P = round(100 * sum(N_PART) / sum(ELIGIBLE), 2)),
        by = c("ROUND", "FC", "SEX")
    ] |>
        prettify_labels()
    dlabs[, `:=`(
        P_LAB = paste0(P, "%"),
        XPOS = Inf,
        YPOS = Inf * (-1 + 2 * as.integer(SEX == "F"))
    )]
    dlabs[, `:=`(
        HJUST  = 1 / 2 + sign(XPOS) / 2,
        VJUST  = 1 / 2 + sign(YPOS) / 2
    )]

    DT |>
        prettify_labels() |>
        ggplot(aes(x = AGEYRS, fill = SEX_LAB)) +
        geom_col(aes(y = ELIGIBLE_PYR), fill = "white", color = "grey60") +
        geom_col(aes(y = N_PART_PYR), color = "grey60") +
        geom_text(data = dlabs[SEX == "F"], aes(x = XPOS, y = YPOS, hjust = 1.2, vjust = 2, label = P_LAB)) +
        geom_text(data = dlabs[SEX == "M"], aes(x = XPOS, y = YPOS, hjust = -0.2, vjust = 2, label = P_LAB)) +
        coord_flip() +
        facet_grid(ROUND_LAB ~ FC_LAB, scales = "free_x") +
        scale_fill_manual(values = palettes$sex) +
        scale_color_manual(values = palettes$sex) +
        scale_y_continuous(labels = abs, expand = c(0.05, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        theme_default() +
        my_labs(y = "Number of census eligible individuals") +
        NULL
}

plot.pyramid.bysexround <- function(DT, NUM, DEN, .ylab) {
    dplot <- copy(DT)
    cols <- c(NUM, DEN)
    new_cols <- paste(c("NUM", "DEN"), "PYR", sep = "_")

    dplot[, (new_cols) := lapply(.SD, function(x) (-1 + 2 * (SEX == "F")) * x), .SDcols = cols]

    dlabs <- dplot[,
        .(P = round(100 * sum(get(NUM)) / sum(get(DEN)), 2)),
        by = c("ROUND", "FC", "SEX")
    ] |>
        prettify_labels()
    dlabs[, `:=`(
        P_LAB = paste0(P, "%"),
        XPOS = Inf,
        YPOS = Inf * (-1 + 2 * as.integer(SEX == "F"))
    )]
    dlabs[, `:=`(
        HJUST  = 1 / 2 + sign(XPOS) / 2,
        VJUST  = 1 / 2 + sign(YPOS) / 2
    )]

    dplot |>
        prettify_labels() |>
        ggplot(aes(x = AGEYRS, fill = SEX_LAB)) +
        geom_col(aes(y = DEN_PYR), fill = "white", color = "grey60") +
        geom_col(aes(y = NUM_PYR), color = "grey60") +
        geom_text(data = dlabs[SEX == "F"], aes(x = XPOS, y = YPOS, hjust = 1.2, vjust = 2, label = P_LAB)) +
        geom_text(data = dlabs[SEX == "M"], aes(x = XPOS, y = YPOS, hjust = -0.2, vjust = 2, label = P_LAB)) +
        coord_flip() +
        facet_grid(ROUND_LAB ~ FC_LAB, scales = "free_x") +
        scale_fill_manual(values = palettes$sex) +
        scale_color_manual(values = palettes$sex) +
        scale_y_continuous(labels = abs, expand = c(0.05, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        theme_default() +
        my_labs(y = .ylab) +
        NULL
}

plot.agecontribution.fromN.stratby <- function(DT, var, by_cols = c('ROUND', 'FC'), .ylab=NA_character_, agegroup = FALSE) 
{

    .xlab <- 'Age'

    if(agegroup){
        .xlab <- 'Age group'
        DT[, AGEGROUP := split.agegroup(AGEYRS)]
        DT <- DT[, lapply(.SD,sum), by=c(by_cols, 'AGEGROUP', 'SEX'), .SDcols=var]
        dcontr <-  DT[, .(P = get(var)/sum(get(var)), AGE=AGEGROUP, SEX=SEX), by=by_cols]
    }else{
        dcontr <- DT[, .(P = get(var)/sum(get(var)), AGE=AGEYRS, SEX=SEX), by=by_cols]
    }

    dcontr |> 
        prettify_labels() |>
        ggplot(aes(x=AGE, y=P, fill=SEX_LAB, color=SEX_LAB)) + {
            if(agegroup){
                geom_col(position='dodge') 
            }else{
                geom_line() 
            }
        } +
        facet_grid(ROUND_LAB~FC_LAB) +
        scale_fill_manual(values=palettes$sex)+ 
        scale_color_manual(values=palettes$sex)+ 
        scale_y_continuous(labels=scales::label_percent(), expand=expansion(mult=c(0,.05))) +
        theme_default() +
        my_labs(y=.ylab, x=.xlab)
}

plot.smoothed.participation.rates <- function(DT)
{
    # DT <- copy(loess_ratepart)
    dplot <- melt.data.table(DT,
        id.vars = c('ROUND' , 'FC' , 'SEX' , 'AGEYRS' ),
        variable.name = "LOESS_SPAN", 
        value.name = "PARTRATE", 
        measure.vars = c('PARTRATE_RAW', 'PARTRATE_SMOOTH.25',  'PARTRATE_SMOOTH.50',  'PARTRATE_SMOOTH.75')) |>
        prettify_labels()
    dplot[, LOESS_SPAN := gsub("PARTRATE_|RAW|SMOOTH", "", LOESS_SPAN) , ]

    ggplot(dplot[LOESS_SPAN!=""], aes(x=AGEYRS, y=PARTRATE, linetype=LOESS_SPAN, color=SEX_LAB)) +
        geom_line() +
        geom_point(data=dplot[LOESS_SPAN==""], aes(linetype=NULL)) + 
        scale_color_manual(values=palettes$sex) + 
        facet_grid( ROUND_LAB ~ FC_LAB ) +
        scale_y_percentage +
        theme_default() +
        my_labs()
}
