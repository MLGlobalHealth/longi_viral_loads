remove_quantiles <- function(DT){   
    set(DT, i=NULL, j=c("CL", "CU", "IL", "IU", "M"), value=NULL)
}

reverse_quantiles <- function(DT){
    set(DT, i=NULL, j=c("CL", "CU", "IL", "IU"), value=DT[, .(CU, CL, IU, IL)])
    cols <- c("CL", "CU", "M", "IL", "IU")
    DT[, (cols) := lapply(.SD, function(x) 1-x), .SDcols=cols]
    return(DT)
}

paper_statements_contributions_viraemia_round <- function(DT= contrib_viraemia_custom, round=19, agegroup="25-39"){
    tmp <- copy(DT)
    tmp[, CELL := prettify_cell(M*100, CL*100, CU*100, percent=TRUE)]
    .round_lab <- drounds[ROUND == as.integer(round), END2]
    tmp[ AGEGROUP == agegroup & ROUND == round, {
        sprintf("By %s, contributions to viraemia by %s aged %s were:\n - %s in Fishing\n - %s in Inland\n", 
            .round_lab,
            unique(SEX), agegroup,
            CELL[which(LOC == 'fishing')],
            CELL[which(LOC == 'inland')] 
        ) |> cat()
        NULL
    }, by=c('SEX')]
} 

paper_statements_overall_prevalence <- function(DT=djoint_agegroup, round=19){

    tmp <- subset(DT, MODEL == 'run-gp-prevl' & AGEGROUP == 'Total')
    tmp[, CELL := prettify_cell(M*100, CL*100, CU*100, percent=TRUE)]


    tmp[ROUND == round & SEX == "Total", {
        sprintf(
            "In round %s, HIV prevalence was estimated to be:\n - %s in inland\n - %s in fishing\n",
            unique(ROUND), CELL[ LOC == 'inland'], CELL[ LOC == 'fishing']
        ) |> cat()
        NULL
    }]
    tmp
}

paper_statements_female_prevalence <- function(DT=djoint_agegroup){

    tmp <- subset(DT, ROUND == 19 & AGEGROUP == 'Total' & SEX == "F" & MODEL == "run-gp-prevl")
    tmp[, CELL := prettify_cell(M*100, CL*100, CU*100, percent=TRUE)]
    tmp[, {
        sprintf("while on average %s of all women had HIV in inland communities .... . Vice versa, while on average %s of all women had HIV in fishing communities, we found that ...\n", 
        CELL[ LOC == 'inland'], CELL[ LOC == 'fishing' ]) |> cat()
    },]
    tmp
}

paper_statements_prevalence_viraemia <- function(DT=djoint_agegroup, model="run-gp-supp-pop"){

    lab <- fifelse(model == "run-gp-supp-pop", yes="census-eligible", no = "PLHIV")
    tmp <- subset(DT, ROUND %in% c(16, 19) & SEX=="Total" & AGEGROUP == 'Total' & MODEL == model)
    tmp[, `:=` ( CELL = prettify_cell(M*100, CL*100, CU*100, percent=TRUE), M=NULL, CL=NULL, CU=NULL, IL=NULL, IU=NULL)]
    tmp[, {
        sprintf("in %s communities, the contribution of unsuppressed individuals to the %s population decreased from %s in round 16 to %s in round 19.\n", 
            unique(LOC),
            lab,
            CELL[ ROUND == 16],
            CELL[ ROUND == 19 ]
        ) |> cat()
    }, by=c("LOC", "SEX")]
    tmp
}

paper_statements_prevalence_viraemia2 <- function(DT=djoint_agegroup, model="run-gp-supp-pop", round = 19){

    model <- match.arg(model, c("run-gp-supp-hiv", "run-gp-supp-pop"))
    if(model == 'run-gp-supp-hiv'){
        stop("this is not the function you want, refer to ... ")
    }
    msg <- fifelse( model == "run-gp-supp-hiv", 
        yes="Suppression among PLHIV was larger in Fishing communities",
        no="Suppression among eligible-individuals was larger in Inland communities")
    dtable <- subset(DT, 
        AGEGROUP %like% 'Total'  & SEX != "Total" &
        ROUND == round & MODEL == model)
    dtable[  , CELL := prettify_cell(M*100, CL*100, CU*100, percent=TRUE) ] |>
        remove_quantiles() |>
        prettify_labels() |>
        remove.nonpretty() |>
        setcolorder(c("LOC_LAB", "ROUND_LAB")) |>
        setkey(LOC_LAB, ROUND_LAB)

    cat(msg, "\n")
    dtable[, {
        sprintf("%s in Fishing versus %s in Inland in %s\n",
            CELL[LOC_LAB %like% 'Fish'],
            CELL[LOC_LAB %like% "Inland"],
            unique(SEX_LAB)) |> cat()
    }, by=c("SEX_LAB")]
    return(dtable)
}

paper_statements_contributions_PLHIV_custom <- function(DT=contrib_plhiv_custom){

    # DT <- copy(contrib_plhiv_custom)
    tmp <- subset(DT, SEX=="Total")
    tmp[, CELL := prettify_cell(M*100, CL*100, CU*100, percent=TRUE)]
    tmp[ AGEGROUP == '30-49', {
        sprintf("In %s communities, the contribution of individuals aged 30-49 increased from:\n%s at baseline to %s in the last round\n", unique(LOC),CELL[ROUND == 16], CELL[ROUND == 19]) |> cat()
    }, by=c('LOC')]
    tmp
}

paper_statements_suppression_above_959595 <- function(DT=djoint){

    tmp <- subset(DT, MODEL == 'run-gp-supp-hiv' & ROUND == 19)
    prettify_labels(tmp)
    sprintf("By 2019, medians for prevalence of suppression exceeded the 95-95-95 goals in:\n") |> cat()
    tmp[ M > .95^3, {
        sprintf("- %s %s aged %s or more\n", unique(LOC_LAB), unique(SEX), min(AGEYRS)) |> cat()
        NULL
    }, by=.(SEX, LOC) ]

    tmp[ LOC == 'inland' & AGEYRS == 25, {
        z <- prettify_cell(M * 100, CL*100, CU*100, percent=TRUE) 
        sprintf("%s for %s aged 25 \n", z,SEX_LAB) |> cat()
        NULL
    }, by = 'SEX_LAB']
    return(tmp)
}

paper_statements_female_contributions_prevalence <- function(DT=dcontrib_agegroup){
    tmp <- subset(DT, MODEL == 'run-gp-prevl' & ROUND == 19 & AGEGROUP == 'Total' & SEX == "F") 
    tmp[, CELL := prettify_cell(M*100, CL*100, CU*100, percent=TRUE)]
    tmp[, {
        sprintf("We found that %s of all individuals with HIV in %s communities were female\n",
        CELL, LOC) |> cat()
        NULL
    }, by=LOC]
    tmp
}

paper_statements_meanage_population <- function(DT=dmeanage, label="run-gp-supp-pop", type = "AGEMEAN"){

    type <- match.arg(type, c("AGEMEAN", "AGESTD"))

    dtable <- subset(DT, ROUND %in% c(16,19) & MODEL == label & TYPE == type)
    dtable[, CELL:=prettify_cell( M, CL, CU, precision=1)] |>
        remove_quantiles()

    fmt <- fifelse(type == "AGEMEAN", 
        yes= "In %s communities, the average age of unsuppressed %s was %s in round 16 and %s in round 19.\n",
        no =  "In %s communities, the standard deviation around the age of unsuppressed %s was %s in round 16 and %s in round 19.\n",
    )
    
    dtable[, { sprintf(
        fmt,
        LOC,
        sex_dictionary2[unique(SEX)],
        CELL[ROUND == 16],
        CELL[ROUND == 19]
    ) |> cat()
    } , by=c('LOC','SEX')] 
    return(dtable)
}

paper_statements_viraemic_among_hiv <- function(DT=djoint, age=30, negate=FALSE, rounds=c(16,19)){

    dtable <- subset(DT, MODEL == 'run-gp-supp-hiv' & AGEYRS == age & ROUND %in% rounds)
    if(negate){ dtable <- reverse_quantiles(dtable) }
    dtable[, CELL := prettify_cell(100*M, 100*CL, 100*CU, percent=TRUE) ]
    remove_quantiles(dtable); prettify_labels(dtable)

    out <- NULL
    if( uniqueN(rounds) == 2){

        word <- fifelse(negate, yes="unsuppressed decreased", no="suppressed increased")
        out <- dtable[ j= {
            sprintf('In %s communities, the proportion of HIV positive %s aged %s years old who were %s from %s in %s to %s in %s\n',
            unique(LOC_LAB), sex_dictionary2[unique(SEX_LAB)], age, word, 
            CELL[ROUND == 16], ROUND_LAB[ROUND == 16],
            CELL[ROUND == 19], ROUND_LAB[ROUND == 19]
        ) |> cat()
            list(CELL=CELL)
        }, by=c("LOC_LAB" ,"SEX_LAB")]

    }else if( uniqueN(rounds) == 1){

        word <- fifelse(negate, yes="unsuppressed was", no="suppressed was")
        out <- dtable[ j= {
            sprintf('In %s communities, the proportion of HIV positive %s aged %s years old who were %s %s\n',
                unique(LOC_LAB), sex_dictionary2[unique(SEX_LAB)], age, word, CELL
            ) |> cat()
            list(CELL=CELL)
        }, by=c("LOC_LAB" ,"SEX_LAB")]
        
    }
}

paper_statements_malefemaleratio_suppression <- function(DT=dmf_ratios,reverse=FALSE){
    .type <- fifelse(reverse==TRUE, yes="VIR", no="SUP")
    dtable <- subset(DT, ROUND %in% c(19) & AGEGROUP == "15-24" & TYPE == .type)
    dtable[, CELL := prettify_cell( M*100, CL*100, CU*100, percent=TRUE)]
    prettify_labels(dtable) |> remove_quantiles()
    dtable[ j=sprintf(
        "In %s communities, the prevalence of %s among 15-25 years old men was %s the one for HIV positive women\n",
        LOC_LAB, 
        fifelse(reverse==TRUE, yes="viraemia", no="suppression"),
        CELL
    ), by=c('ROUND_LAB', "LOC_LAB")]
}

paper_statements_malefemaleratio_suppression2 <- function(DT=dmf_ratios){
    dtable <- subset(DT, AGEGROUP %like% "Total" & TYPE == "VIR") 
    dtable[  , CELL := prettify_cell(M, CL, CU, precision=1) ] |>
        remove_quantiles() |>
        prettify_labels() |>
        remove.nonpretty() |>
        setcolorder(c("LOC_LAB", "ROUND_LAB")) |>
        setkey(LOC_LAB, ROUND_LAB)

    .round_lab <- drounds[ROUND == 19, END2]
    dtable[ ROUND_LAB %like% '19', {
        sprintf("In %s communities, by %s, an estimated %s  times more men than women remained unsuppressed\n",
            LOC_LAB,.round_lab, CELL) |> cat()
    }, by = c("LOC_LAB")]
    return(dtable)
}

paper_statements_suppression_PLHIV_aggregated <- function(DT=dsupp_agegroup, reverse=FALSE, round=19){
    word <- fifelse(reverse==TRUE, yes="non-suppression", no="suppression")
    dtable <- subset(DT, AGEGROUP == "Total" & SEX != "Total" & ROUND %in% round)
    if(reverse){ dtable <- reverse_quantiles(dtable) }
    dtable[, CELL := prettify_cell( M*100, CL * 100, CU * 100, percent=TRUE)]
    remove_quantiles(dtable) |> prettify_labels()
    if(length(round) == 2){
        dtable[ , sprintf("In %s communities, prevalence of %s in HIV positive %s was %s by 2013 and %s by 2019\n",
            unique(LOC_LAB), unique(word), unique(SEX_LAB),
            CELL[ROUND == 16], CELL[ROUND == 19]
        ) |> cat()
        , by = c("LOC_LAB", "SEX_LAB")]
    }else if(length(round) == 1){
        dtable[ , sprintf("In %s communities, prevalence of %s in HIV positive %s was %s in Round %s\n",
            unique(LOC_LAB), unique(word), unique(SEX_LAB),
            CELL[ROUND == round], round
        ) |> cat()
        , by = c("LOC_LAB", "SEX_LAB")]
    }
    return(dtable)
}

print_statements_half_plhiv <- function(DT=dcontrib_50p_PLHIV, DAGES=plhiv_contributors){
    DT[, CELL := prettify_cell( M * 100, CL * 100, CU * 100) ]
    tab <- DT |> 
        prettify_labels() |>
        remove_quantiles() |>
        set(j="INAGEGROUP", value=NULL) |>
        merge(DAGES, c('SEX', 'LOC')) 
    catn("In round 19:")
    tab[, sprintf(
        "In %s communities, %s of %s with HIV were aged between %s and %s",
        LOC_LAB, CELL, SEX_LAB, MIN, MAX), 
    by=c("SEX_LAB", "LOC_LAB")]
}

paper_statements_contributions_census_eligible <- function(DT=dcens, comm="inland", round=19, smooth=FALSE, agegroup=NA_character_, sex=NA_character_){

    sex <- match.arg(sex, c(NA_character_, "F", "M"))
    comm <- match.arg(comm, c("inland", "fishing"))
    n_var <- fifelse(smooth, yes="ELIGIBLE_SMOOTH", no="ELIGIBLE")

    stopifnot(agegroup %like% '[0-9]{2}-[0-9]{2}')
    tmp <- strsplit(agegroup, "-")[[1]]
    age_min <- as.numeric(tmp[1])
    age_max <- as.numeric(tmp[2])

    contrib <- DT[ ROUND == round & LOC == comm, {
        idx <- AGEYRS %between% c(age_min, age_max) & ( is.na(sex) | SEX == sex ) 
        sum(.SD[idx])/sum(.SD)
    }, .SDcols = n_var ]

    sprintf(
        "In round %s, %s aged %i %i accounted for %.1f%% of the CE population\n",
        round, sex, age_min, age_max, contrib*100) |> cat()
}

paper_statements_average_participation <- function(DT=dprop){
    fmt <- paste0(
        "there were on average %.0f censused inds, of who on average", 
        " %.0f  (%.1f %%) individuals participated in the RCCS\n"
        )

    tmp <- DT[, 
        .(EL = sum(ELIGIBLE), PART=sum(N_PART)),
    by=c("FC","SEX","ROUND")][, PERC := PART/EL ]
    # tmp[j=lapply(.SD, mean), .SDcols = c("EL", "ROUND")][j=sprintf(fmt, EL, PART, 100*PERC)] |> cat()
    tmp[, PERC := round(100*PERC, 1)]
    return(tmp)
}

paper_statements_prevalence_viraemia_maximum <- function(DT, round = 19){
    
    # DT <- dsupp_agegroup_custom ; round <- 19
    fmt <- "In %s living %s communities, prevalence of viraemia was highest in %s yo [ %s ]\n"
    dtable <- subset(DT, ROUND == round) |> 
        reverse_quantiles() |> 
        prettify_labels() 
    dtable[, CELL := prettify_cell(M*100, CL*100, CU*100, percent=TRUE)]
    dtable[, {
        idx <- which.max(M)
        sprintf( fmt, unique(SEX_LAB), unique(LOC_LAB), AGEGROUP[idx], CELL[idx]) |> 
        cat()
    }, by=.(LOC_LAB, SEX_LAB)]

}

make.table.eligible.participants <- function(DT=tab_el, cols=c("EL", "PART"), splitrow=TRUE){

    # DT <- copy(tab_el)
    dtable <- copy(DT)
    if('COMM' %in% names(dtable)){
        setnames(dtable, "COMM", "FC")
    }
    setkey(dtable, FC, SEX, ROUND)
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

    if(splitrow){
        out[, DUMMY := 1:.N]
        out <- out[, {
            filler <- "communities"
            filler <- fifelse(Location %like% filler, yes=filler, no="")
            rbind(.SD, data.table(Location=filler), fill=TRUE)
        }, by=DUMMY]
        out[,  `Survey\nRound` := strsplit(`Survey\nRound`, "\\n"), by="DUMMY"]
        out[, DUMMY:=NULL]
        out[is.na(out)] <- ""
        out[, `Location`:=gsub(" communities$", "", x=`Location`)]
    }
    xtable(out) |> print()
    return(out)
}

make.table.first.time.participants <- function(DT=check_r1619, splitrow=TRUE){

    by_cols <- c("COMM", "SEX", "ROUND")
    dtable_ftp <- DT[ FIRST_PARTICIPATION == TRUE, .(
        N=  uniqueN(STUDY_ID), 
        N_HIV= sum(HIV_STATUS == 1, na.rm = TRUE),
        N_UV = sum(HIV_VL > 1000, na.rm = TRUE)
    ), by=by_cols]
    setkeyv(dtable_ftp,by_cols)
    dtable_ftp[, `:=` (
        R_HIV = round(N_HIV/N_HIV[1], 2),
        R_UV = round(N_UV/N_UV[1], 2)
    ), by=c("COMM", "SEX")]
    cols <- c("N", "N_HIV", "N_UV")
    dtable_ftp[, (cols) := lapply(.SD, comma), .SDcols=cols]
    dtable_ftp[, ROUND := round_labs3[as.character(ROUND)]]

    out <- prettify_labels(dtable_ftp) |> 
        set(j=c("ROUND_LAB", "SEX"), value=NULL) |>
        setkey(COMM, SEX_LAB, ROUND)|>
        setcolorder(c("COMM", "SEX_LAB", "ROUND")) |>
        delete.repeated.table.values() |>
        setnamesdict(dict=dict_table_names$first_time_participants)
    out[, Location := gsub("([g|d])$", "\\1 communities", Location) ]
    out[, Gender := gsub("Female", "Women", Gender) |> 
                gsub("Male", "Men", x=_) ]
    out

    if(splitrow){
        out[, DUMMY := 1:.N]
        out <- out[, {
            filler <- "communities"
            filler <- fifelse(Location %like% filler, yes=filler, no="")
            rbind(.SD, data.table(Location=filler), fill=TRUE)
        }, by=DUMMY]
        out[,  `Survey\nRound` := strsplit(`Survey\nRound`, "\\n"), by="DUMMY"]
        out[, DUMMY:=NULL]
        out[is.na(out)] <- ""
        out[, `Location`:=gsub(" communities$", "", x=`Location`)]
    }
    xtable(out) |> print()
    return(out)
}
