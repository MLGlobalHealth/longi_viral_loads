paper_statements_contributions_viraemia_round19 <- function(DT= contrib_viraemia_custom){
    tmp <- copy(DT)
    tmp[, CELL := prettify_cell(M*100, CL*100, CU*100, percent=TRUE)]
    tmp[ AGEGROUP == '25-39', {
        sprintf("By 2019, contributions to viraemia by %s aged 25-40 were:\n - %s in Fishing\n - %s in Inland\n", 
        unique(SEX), CELL[which(LOC == 'fishing')], CELL[which(LOC == 'inland')] ) |> cat()
        NULL
    }, by=c('SEX')]
} 

paper_statements_overall_prevalence <- function(DT=djoint_agegroup){

    tmp <- subset(DT, MODEL == 'run-gp-prevl' & AGEGROUP == 'Total')
    tmp[, CELL := prettify_cell(M*100, CL*100, CU*100, percent=TRUE)]

    tmp[ROUND == 16 & SEX == "Total", {
        sprintf(
            "At baseline, HIV prevalence was estimated to be:\n - %s in inland\n - %s in fishing\n",
            CELL[ LOC == 'fishing'], CELL[ LOC == 'inland' ]
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
    sprintf("By 2019, medians for prevalnce of suppression exceeded the 95-95-95 goals in:\n") |> cat()
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
