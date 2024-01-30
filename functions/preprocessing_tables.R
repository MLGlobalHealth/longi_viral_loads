make.table.firstparticipant.NPhiv.NPunsupp <- function(DT=dall){

    dt_sub <- subset(dall, FIRST_PARTICIPATION == 1)

    dtable <- dt_sub[, j=list( 
        N = .N,
        N_HIV=sum(HIV_STATUS, na.rm=TRUE),
        N_UNSUPP=sum(VL_COPIES >= VIREMIC_VIRAL_LOAD, na.rm=TRUE)
    ), by=c('FC', 'SEX', "ROUND") ]

    dtable[, `:=` (
        RATIO_HIV_R16 = round((N_HIV/N_HIV[ROUND==16]),2) ,
        RATIO_UNSUPP_R16 = round((N_UNSUPP/N_UNSUPP[ROUND==16]),2)
    ), by=.(FC, SEX)]

    out <- prettify_labels(dtable) |> 
        remove.nonpretty() |>
        setkey(FC_LAB, SEX_LAB, ROUND_LAB)|>
        setcolorder(c("FC_LAB", "SEX_LAB", "ROUND_LAB")) |>
        delete.repeated.table.values() |>
        setnamesdict(dict=my_labs_dictionary) |>
        setnamesdict(dict=dict_table_names$percent_reduction)

    tmp_dict <- c(
        "Community type" = "Community \\ type \\  \\  \\  \\  \\  \\  \\ (n)",
        "Gender" = "Gender \\  \\  \\  \\  \\  \\  \\  \\ (n)",
        "Survey round" = "Survey \\ round \\  \\  \\  \\  \\  \\  \\ (n)",
        "N" = "First-time \\ participants \\  \\  \\  \\  \\  \\  \\ (n)",
        "N_HIV" = "First-time \\ participants \\ with HIV \\  \\  \\  \\  \\  \\ (n)",
        "N_UNSUPP" = "First-time \\ participants \\ with \\ unsuppressed \\ virus \\  \\  \\  \\ (n)",
        "RATIO_HIV_R16" = "Ratio of \\ firt time \\ participants \\ with HIV \\ relative to \\ round 16 \\  \\  \\ (p)",
        "RATIO_UNSUPP_R16" = "Ratio of \\ firt time \\ participants \\ with \\ unsuppressed \\ virus \\ relative to \\ round 16 \\ (p)"
    )

    # could add \textbf here and tabular here...

    setnamesdict(out, dict=tmp_dict)

    print(xtable(out, type = "latex"))

    return(out)
}

make.table.eligible.participants <- function(DT=tab_el){

    # DT <- copy(tab_el)
    dtable <- copy(DT)
    setkey(dtable, FC, SEX, ROUND)
    cols <- c("EL", "PART")
    dtable[, (cols) := lapply(.SD, comma), .SDcols=cols]
    dtable[, PERC := paste0(PERC, "%")]
    dtable[, ROUND := round_labs3[as.character(ROUND)]]
    out <- prettify_labels(dtable) |> 
        set(j=c("ROUND_LAB", "FC", "SEX"), value=NULL) |>
        setkey(FC_LAB, SEX_LAB, ROUND)|>
        setcolorder(c("FC_LAB", "SEX_LAB", "ROUND")) |>
        delete.repeated.table.values() |>
        setnamesdict(dict=dict_table_names$eligible_participants)
    out[, Location := gsub("([g|d])$", "\\1 communities", Location) ]
    out[, Gender := gsub("Female", "Women", Gender) |> 
                gsub("Male", "Men", x=_) ]
    xtable(out) |> print()
    return(out)
}
