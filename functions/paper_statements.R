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
    tmp[, .(ROUND, LOC, SEX, CELL)]

    tmp[ROUND == 16 & SEX == "Total", {
        sprintf(
            "At baseline, HIV prevalence was estimated to be:\n - %s in inland\n - %s in fishing\n",
            CELL[ LOC == 'fishing'], CELL[ LOC == 'inland' ]
        ) |> cat()
        NULL
    }]
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
