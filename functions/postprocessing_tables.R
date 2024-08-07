tablify.posterior.Nunsuppressed <- function(list, model = "run-gp-supp-pop", CELLname = "UNSUPP") {
    # number unsuppressed:
    tab_list <- copy(list)
    tab_list <- lapply(tab_list, remove.ILIU)
    tab_list <- lapply(tab_list, function(DT) DT[MODEL == model])

    # report as integer, with 1,000 comma separators
    cols <- c("M", "CL", "CU")
    tab_list$round_totals[.SDcols=cols, 
        j = (cols) := lapply(.SD, function(x){ comma(as.integer(x))})]

    with(tab_list, {
        round_totals[, CELL_N := sprintf("%s (%s-%s)", M, CL, CU)]
        percent_reduction[, CELL_P := prettify_cell(M * 100, CL * 100, CU * 100, newline = FALSE, percent = TRUE)]

        # select what to show per round
        percent_reduction[, ROUND := gsub("redR([0-9]{2})", "\\1", variable) |> as.integer(), ]

        .s <- function(DT) {
            subset(DT, select = intersect(c("SEX", "LOC", "ROUND", "CELL_N", "CELL_P"), names(DT)))
        }
        merge(.s(round_totals), .s(percent_reduction), by = c("SEX", "LOC", "ROUND"), all.x = TRUE)
    }) -> tab

    prettify_labels(tab)

    keys <- c("LOC_LAB", "SEX_LAB", "ROUND_LAB")
    tab <- subset(tab, select = c(keys, c("CELL_N", "CELL_P"))) |>
        setkeyv(keys)

    setnames(tab, c("CELL_N", "CELL_P"), paste("CELL", CELLname, c("N", "P"), sep = "_"))

    # tab2 <- delete.repeated.table.values(tab, cols = keys) |>
    #     table.na.to.empty(cols=keys)
    return(tab)
}

tablify.agecontributions <- function(DT, model = "run-gp-supp-pop") {
    # DT <- copy(dcontrib_agegroup) ; model= 'run-gp-supp-pop'
    dtab <- subset(DT, MODEL == model)
    remove.ILIU(dtab) |>
        prettify_labels()
    dtab[, CELL := prettify_cell(M * 100, CL * 100, CU * 100, percent = TRUE, newline = TRUE)]

    keys <- c("LOC_LAB", "ROUND_LAB", "SEX_LAB", "AGEGROUP")
    dtab <- subset(dtab, select = c(keys, c("CELL"))) |>
        setkeyv(keys)

    dtab <- dcast.data.table(dtab,
        LOC_LAB + ROUND_LAB ~ SEX_LAB + AGEGROUP,
        value.var = "CELL"
    ) |>
        delete.repeated.table.values(cols = c("LOC_LAB", "ROUND_LAB"))
    return(dtab)
}

make.table.Nhivpositive <- function(DT = joint_ageagrr_list$round_totals, DC = dcens) {
    # aggregate number of census-eligible individuals by gender, loc, sex
    dcens_aggr <- DC[
        j = lapply(.SD, sum),
        .SDcols = c("ELIGIBLE_SMOOTH", "ELIGIBLE"),
        by = c("ROUND", "LOC", "SEX")
    ]

    tab <- merge(
        DT[MODEL == "run-gp-prevl"],
        dcens_aggr,
        by = c("ROUND", "LOC", "SEX")
    ) |> remove.ILIU()
    cols <- c("CL", "M", "CU")
    tab[, (cols) := lapply(.SD, function(x) x / ELIGIBLE), .SDcols = cols]
    tab[, CELL_NHIV := prettify_cell(M * 100, CL * 100, CU * 100, percent = TRUE)]
    tab[ROUND == 19, {
        sprintf(
            "In round 19, %s of men and %s of women were estimated to be HIV positive in fishing communities, \n compared to %s and %s among men and women in inland, respectively.\n",
            CELL_NHIV[LOC == "fishing" & SEX == "M"], CELL_NHIV[LOC == "fishing" & SEX == "F"],
            CELL_NHIV[LOC == "inland" & SEX == "M"], CELL_NHIV[LOC == "inland" & SEX == "F"]
        ) |> cat()
    }]

    # setnames(tab,
    #     c(""),
    #     c("")
    # )


    tab <- subset(tab, select = c("ROUND", "LOC", "SEX", "CELL_NHIV"))
    prettify_labels(tab)
    return(tab)
}

make_main_table_contributions <- function(
    DTPOP = ncen, 
    DJOINT = djoint_agegroup,
    DCONTRIB = dcontrib_agegroup,
    DSUPP = dsupp_agegroup,
    add_asterisks_unaids = TRUE,
    include_totals = c("SEX", "AGEGROUP"),
    round = 19
) {
    include_totals <- intersect(include_totals, c("SEX", "AGEGROUP"))
    if (all(include_totals == "SEX")) {
        warning("SEX totals can only be included if AGEGROUP is. Including both...")
        include_totals <- c("SEX", "AGEGROUP")
    }

    # get info for census-eligible pop size
    DTPOP[, AGEGROUP := split.agegroup(AGEYRS)]
    ncen_agegroup <- DTPOP[ROUND == 19, list(ELIGIBLE = sum(ELIGIBLE)), by = c("ROUND", "FC", "SEX", "AGEGROUP")] |> setnames("FC", "LOC")
    # get totals
    ncen_totals <- groupingsets(ncen_agegroup,
        j = list(ELIGIBLE = sum(ELIGIBLE)),
        sets = list("LOC", c("LOC", "SEX")),
        by = c("LOC", "SEX", "AGEGROUP"),
    )[, lapply(.SD, na2string, str = "Total")]
    ncen_totals[, ROUND := 19]
    # compute nice cells N(P%)
    expr_eligible_cell <- expr({
        sprintf("%s (%.2f%%)", 
            format(ELIGIBLE, big.mark = ","), 
            100 * proportions(ELIGIBLE))
    })
    cols <- c("ROUND", "LOC")
    ncen_agegroup[SEX != "Total" & AGEGROUP != "Total", ELIGIBLE_CELL := eval(expr_eligible_cell), by = cols]
    ncen_totals[SEX != "Total" & AGEGROUP == "Total", ELIGIBLE_CELL := eval(expr_eligible_cell), by = c(cols)]
    ncen_totals[SEX == "Total" & AGEGROUP == "Total", ELIGIBLE_CELL := eval(expr_eligible_cell), by = c(cols)]
    # merge together
    ncen_agegroup <- rbind(ncen_agegroup, ncen_totals, use.names = TRUE)

    e_sex <- expr(ROUND == round & SEX != "Total")
    if ("SEX" %in% include_totals) {
        e_sex <- expr(ROUND == round)
    }
    e_age <- expr(AGEGROUP != "Total")
    if ("AGEGROUP" %in% include_totals) {
        e_age <- expr(TRUE)
    }

    # get info from gp fits
    r19_comp_usnupp <- DCONTRIB[MODEL == "run-gp-supp-pop" & eval(e_sex) & eval(e_age)]
    r19_comp_hiv <- DCONTRIB[MODEL == "run-gp-prevl" & eval(e_sex) & eval(e_age)]
    r19_hivprev <- DJOINT[MODEL == "run-gp-prevl" & eval(e_sex) & eval(e_age)]
    r19_prop_unsupp <- DJOINT[MODEL == "run-gp-supp-pop" & eval(e_sex) & eval(e_age)]
    r19_supphiv <- DSUPP[eval(e_sex) & eval(e_age)]
    ncen_agegroup <- ncen_agegroup[eval(e_sex) & eval(e_age)]

    # go from supp to unsupp
    negate.percent.quantiles(r19_supphiv) 
    if (add_asterisks_unaids) {
        r19_supphiv[, `:=`(
            UNAIDS_achieved = fifelse(M > .95^3, yes = " * ", no = ""),
            M = NULL, CL = NULL, CU = NULL, IL = NULL, IU = NULL
        )]
    } else {
        NULL
    }

    # check proportions sum to 100% 
    if (length(include_totals) == 0) {
        check <- lapply(
            list(r19_comp_usnupp, r19_comp_hiv),
            function(DT) {
                DT[AGEGROUP != "Total" & SEX != "Total", 
                sum(M) %between% c(.95, 1.05) |> stopifnot(),
                by = c("ROUND", "LOC")]
            }
        )
    }

    # prettify cells 
    null <- lapply(
        list(r19_hivprev, r19_prop_unsupp, r19_comp_usnupp, r19_supphiv, r19_comp_hiv),
        function(DT) {
            DT[, CELL := prettify_cell(M * 100, CL * 100, CU * 100, percent = TRUE)]
            DT[, `:=`(MODEL = NULL, M = NULL, CL = NULL, CU = NULL, IU = NULL, IL = NULL)]
        }
    ); rm(null)

    totals <- expr(AGEGROUP == "Total" & SEX == "Total")

    # prettify all
    longname <- FALSE
    `%+%` <- function(x, y){ paste0(x, y)}
    endings <- c(
        N = "\n\n(n, (% of total population))",
        p = "\n\n(posterior median\nestimate, (95%CrI))"
    )

    names_comp <- if(longname) {
        c(
            hiv = "Age profile of HIV burden",
            unsupp = "Age profile of viraemia burden"
        )
    } else {
        c(
            hiv = "Age profile of people\nliving with HIV",
            unsupp = "Age profile of people\nexhibiting viraemia"
        )
    }
    setnames(ncen_agegroup, "ELIGIBLE_CELL", "Census-eligible\nindividuals in\neach age group" %+% endings["N"])
    if( longname ){
        setnames(r19_hivprev, "CELL", "HIV seroprevalence in age band" %+% endings["p"])
        setnames(r19_prop_unsupp, "CELL", "Proportion of census-eligible individuals\nwho are unsuppressed in age  band" %+% endings["p"])
        setnames(r19_supphiv, "CELL", "Proportion of people\nliving with HIV\nexhibiting viraemia in\neach age group" %+% endings["p"])
    }else{
        setnames(r19_hivprev, "CELL", "Proportion of\npeople living with\nHIV in each age\ngroup" %+% endings["p"])
        setnames(r19_prop_unsupp, "CELL", "Proportion of all\npersons with viraemia\nin each age group" %+% endings["p"])
        setnames(r19_supphiv, "CELL", "Proportion of people\nliving with HIV\nexhibiting viraemia in\neach age group" %+% endings["p"])
    }
    setnames(r19_comp_hiv, "CELL", names_comp["hiv"] %+% endings["p"])
    setnames(r19_comp_usnupp, "CELL", names_comp["unsupp"] %+% endings["p"])
    dtable <- Reduce(
        f = function(x, y) merge(x, y, all = TRUE, by = c("ROUND", "LOC", "SEX", "AGEGROUP")),
        x = list(
            ncen_agegroup,
            r19_hivprev,
            r19_supphiv,
            r19_prop_unsupp,
            r19_comp_hiv,
            r19_comp_usnupp
        )
    )
    # if ("SEX" %in% include_totals) {
    #     cols <- names(dtable) %like% paste(names_comp, collapse = "|")
    #     dtable[SEX == "Total", 
    #     (cols) := lapply(.SD, function(x) return("100.00%")), .SDcols = cols]
    # }
    prettify_labels(dtable)
    dtable[, SEX_LAB := sex_dictionary2[SEX_LAB]]
    cols <- c("LOC_LAB", "SEX_LAB", "AGEGROUP")
    setcolorder(dtable, cols)
    if (add_asterisks_unaids) {
        dtable[, (names_comp["unsupp"]) := lapply(.SD, \(x) paste0(x, UNAIDS_achieved)), .SDcols = names_comp["unsupp"]]
    }
    dtable[, `:=`(ROUND_LAB = NULL, LOC = NULL, SEX = NULL, ELIGIBLE = NULL, ROUND = NULL, UNAIDS_achieved = NULL)]

    # put first column last (so that Totals are on top)
    dtable <- dtable[, rbind(.SD[.N], .SD[1:(.N - 1)]), by = c("LOC_LAB", "SEX_LAB")]
    dtable <- unique(dtable)
    dtable <- delete.repeated.table.values(dtable, cols)
    # change col names for _LAB variables 
    setnames(dtable,
        names(my_labs_dictionary),
        unname(my_labs_dictionary),
        skip_absent = TRUE
    )

    # final changes
    setnames(dtable, "Age group", "Age\ngroup")
    dtable
}

make.supp.table.meanage <- function(DT=dmeanage, label=NA_character_, newline=TRUE, include_date=TRUE){

    if(include_date)
        newline <- TRUE

    tab <- if(! is.na(label)){
        subset(DT,  MODEL == label)
    }else{ subset(DT) } |> subset(TYPE %like% "AGE[0-9]{2}")

    # get entries
    tab <- tab[, CELL := prettify_cell(M, IL, IU)] |>
        remove_quantiles() |>
        prettify_labels() |>
        remove.nonpretty() |>
        dcast(  LOC_LAB + SEX_LAB + ROUND_LAB ~ TYPE + MODEL, value.var = "CELL") |>
        delete.repeated.table.values(cols=c("LOC_LAB", "SEX_LAB", "ROUND_LAB"))

    # change colnames 
    setnames(tab, 
        old=names(dict_table_names$mean_ages),
        new=unname(dict_table_names$mean_ages), 
        skip_absent = TRUE)

    if(newline){
        nms <- setdiff(names(tab), c("Community type", "Gender", "Round"))
        tab[, DUMMY := 1:.N]
        tab <- tab[, lapply(.SD, function(x){ 
            z <- stringr::str_split_1(x, " \\(") 
            if(length(z)==2){
                z[2] <- paste0("(", z[2]) 
            }else{z[2] <- ""}
            z
        }),  by=DUMMY]
        tab[, DUMMY := NULL]
    }

    if(include_date){
        .labs <- drounds[, paste(START, END, sep="-")]
        .labs <- gsub("20","", .labs)
        tab[Round == "", Round := rep(.labs, times=.N/length(.labs)) ]
    }

    # sort colnames
    order <- intersect(unname(dict_table_names$mean_ages), names(tab))
    setcolorder(tab, order)
    print(
        xtable(tab),
        hline.after = which(tab$`Location` == "Inland") - 1,
        add.to.row = list(
            pos = list(which(tab$Gender == "Male") - 1),
            command = sprintf("\\cline{2-%i}", ncol(tab))
        )
    )
    return(tab)
}

make.table.malefemale.diff.viraemia <- function(DT){

    dtable <- subset(DT, TYPE == "SUP-DIFF" & AGEGROUP == "Total")
    prettify_labels(dtable)
    cols <- c("M", "CL", "CU")
    dtable[, (cols) := lapply(.SD, function(x) x * 100), .SDcols=cols, by=c("LOC_LAB", "ROUND_LAB")]
    dtable[, CELL_MF_GAP := sprintf("%.1f%% (%.1f, %.1f)", M, CL, CU)]
    dtable[, SEX_LAB := "Male"]
    tmp0 <- subset(dtable, select=c("LOC_LAB","SEX_LAB","ROUND_LAB", "CELL_MF_GAP"))
    tmp1 <- copy(tmp0)[, `:=` (SEX_LAB = "Female", CELL_MF_GAP = "--")]
    return( rbind(tmp0, tmp1) )
}

