tablify.posterior.Nunsuppressed <- function(list, model = 'run-gp-supp-pop'){

    # number unsuppressed: 
    tab_list <- copy(list)
    tab_list <- lapply(tab_list, remove.ILIU)
    tab_list <- lapply(tab_list, function(DT) DT[MODEL==model])

    with(tab_list, {

        round_totals[ , CELL_N := prettify_cell(M, CL, CU, precision = 0, newline = FALSE)]
        percent_reduction[ , CELL_P := prettify_cell(M*100, CL*100, CU*100, newline = FALSE, percent=TRUE)]
        
        # select what to show per round
        percent_reduction[, ROUND := gsub( 'redR([0-9]{2})', '\\1', variable) |> as.integer() ,]

        .s <- function(DT) 
            subset(DT, select=intersect(c('SEX', 'LOC', 'ROUND', 'CELL_N', 'CELL_P'), names(DT)))
        tab <<- merge(.s(round_totals), .s(percent_reduction), by=c('SEX', 'LOC', 'ROUND'), all.x=TRUE )
    })

    prettify_labels(tab)
    keys <- c('LOC_LAB', 'SEX_LAB', 'ROUND_LAB')
    tab <- subset(tab, select=c(keys, c('CELL_N', 'CELL_P'))) |> 
        setkeyv(keys)

    tab2 <- delete.repeated.table.values(tab, cols = keys) |> 
        table.na.to.empty(cols=keys)
    return(tab2)
}

tablify.agecontributions <- function(DT, model= 'run-gp-supp-pop'){

    # DT <- copy(dcontrib_agegroup) ; model= 'run-gp-supp-pop'
    dtab <- subset(DT, MODEL==model)
    remove.ILIU(dtab) |> 
        prettify_labels()
    dtab[ , CELL := prettify_cell(M*100, CL*100, CU*100, percent=TRUE, newline=TRUE) ]

    keys <- c('LOC_LAB',  'ROUND_LAB', 'SEX_LAB', 'AGEGROUP')
    dtab <- subset(dtab, select=c(keys, c('CELL'))) |> 
        setkeyv(keys)

    dtab <- dcast.data.table(dtab,
        LOC_LAB + ROUND_LAB ~ SEX_LAB + AGEGROUP,
        value.var = 'CELL'
    ) |> 
        delete.repeated.table.values(cols=c('LOC_LAB','ROUND_LAB'))
    return(dtab)
}
