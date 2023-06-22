paper_statements_contributions_viraemia_round19 <- function(DT= contrib_viraemia_custom){
    tmp <- copy(DT)
    tmp[, CELL := prettify_cell(M*100, CL*100, CU*100, percent=TRUE)]
    tmp[ AGEGROUP == '25-39', {
        sprintf("By 2019, contributions to viraemia by %s aged 25-40 were:\n - %s in Fishing\n - %s in Inland\n", 
        unique(SEX), CELL[which(LOC == 'fishing')], CELL[which(LOC == 'inland')] ) |> cat()
        NULL
    }, by=c('SEX')]
} 
