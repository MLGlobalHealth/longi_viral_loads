remove_quantiles <- function(DT){   
    set(DT, i=NULL, j=c("CL", "CU", "IL", "IU", "M"), value=NULL)
}

reverse_quantiles <- function(DT){
    set(DT, i=NULL, j=c("CL", "CU", "IL", "IU"), value=DT[, .(CU, CL, IU, IL)])
    cols <- c("CL", "CU", "M", "IL", "IU")
    DT[, (cols) := lapply(.SD, function(x) 1-x), .SDcols=cols]
    return(DT)
}

paper_statements_contributions_viraemia_round <- function(DT= dmf_ratios, round=19, agegroup="25-39"){
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

paper_statements_prevalence_viraemia <- function(DT=djoint_agegroup){

    tmp <- subset(DT, ROUND %in% c(16, 19) & SEX=="Total" & AGEGROUP == 'Total' & MODEL == "run-gp-supp-pop")
    tmp[, `:=` ( CELL = prettify_cell(M*100, CL*100, CU*100, percent=TRUE), M=NULL, CL=NULL, CU=NULL, IL=NULL, IU=NULL)]
    tmp[, {
        sprintf("in %s communities, the contribution of viraemic individuals to the census eligible population decreased from %s in round 16 to %s in round 19.\n", 
            unique(LOC),
            CELL[ ROUND == 16],
            CELL[ ROUND == 19 ]
        ) |> cat()
    }, by=c("LOC", "SEX")]
    tmp
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
        yes= "In %s communities, the average age of %s with viraemia was %s in round 16 and %s in round 19.\n",
        no =  "In %s communities, the standard deviation around the age of %s with viraemia was %s in round 16 and %s in round 19.\n",
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

paper_statements_suppression_PLHIV_aggregated <- function(DT=dsupp_agegroup, reverse=FALSE){
    word <- fifelse(reverse==TRUE, yes="non-suppression", no="suppression")
    dtable <- subset(DT, AGEGROUP == "Total" & SEX != "Total" & ROUND %in% c(16, 19))
    if(negate){ dtable <- reverse_quantiles(dtable) }
    dtable[, CELL := prettify_cell( M*100, CL * 100, CU * 100, percent=TRUE)]
    remove_quantiles(dtable) |> prettify_labels()
    dtable[ , sprintf("In %s communities, prevalence of %s in HIV positive %s was %s by 2013 and %s by 2019",
        unique(LOC_LAB), unique(word), unique(SEX_LAB),
        CELL[ROUND == 16], CELL[ROUND == 19]
    )
    , by = c("LOC_LAB", "SEX_LAB")]
    dtable
}
