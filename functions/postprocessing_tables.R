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

make.main.table.contributions <- function(DTPOP=ncen, DJOINT=djoint_agegroup, DCONTRIB=dcontrib_agegroup, DSUPP=dsupp_agegroup, add_asterisks_unaids=TRUE){

    # get info for census eligible pop size
    DTPOP[, AGEGROUP:=split.agegroup(AGEYRS)]
    ncen_agegroup <- DTPOP[ ROUND == 19, list(ELIGIBLE=sum(ELIGIBLE)), by=c('ROUND', 'FC', 'SEX','AGEGROUP')] |>
        setnames("FC", "LOC")
    ncen_agegroup[, ELIGIBLE_CELL:= sprintf(
            "%s (%.2f%%)",
            ELIGIBLE,
            100*proportions(ELIGIBLE)
        ), by=c("ROUND", "LOC")]

    # get info from gp fits
    r19_hivprev <- DJOINT[ ROUND == 19 & MODEL == 'run-gp-prevl' & AGEGROUP != "Total"] 
    r19_prop_unsupp <- DJOINT[ ROUND == 19 & MODEL == "run-gp-supp-pop" & AGEGROUP != "Total"]
    r19_comp_usnupp <- DCONTRIB[ ROUND == 19 & MODEL=='run-gp-supp-pop' & AGEGROUP != "Total"]
    r19_comp_hiv <- DCONTRIB[ ROUND == 19 & MODEL=='run-gp-prevl' & AGEGROUP != "Total"]
    r19_supphiv <- DSUPP[ROUND == 19] 
    if(add_asterisks_unaids){
        r19_supphiv[, `:=` ( UNAIDS_achieved = fifelse(M >  .95^3, yes=" * ", no=""), M=NULL, CL=NULL, CU=NULL, IL=NULL, IU=NULL)] 
    }else{NULL}

    check <- lapply(
        list(r19_comp_usnupp, r19_comp_hiv),
        function(DT) {
            DT[, sum(M) %between% c(.95, 1.05) |> stopifnot(), by=c("ROUND", "LOC")]
        })
    null <- lapply(
        list(r19_hivprev, r19_prop_unsupp, r19_comp_usnupp, r19_comp_hiv),
        function(DT) {
            DT[, CELL := prettify_cell(M*100, CL*100, CU*100, percent=TRUE)]
            DT[, `:=` (MODEL=NULL,M=NULL, CL=NULL, CU=NULL, IU=NULL,IL=NULL)]
        }
    ); rm(null)
    # prettify all
    setnames(ncen_agegroup, "ELIGIBLE_CELL", "Population in age band")
    setnames(r19_hivprev, "CELL", "People with HIV in age band")
    setnames(r19_prop_unsupp, "CELL", "People with unsuppressed HIV in age  band")
    setnames(r19_comp_hiv, "CELL", "Age composition of people with HIV")
    setnames(r19_comp_usnupp, "CELL", "Age composition of people with unsuppressed HIV")
    dtable <- Reduce(
        f=function(x,y) merge(x,y, by=c("ROUND", "LOC", "SEX", "AGEGROUP")),
        x=list(ncen_agegroup, r19_hivprev, r19_prop_unsupp, r19_comp_hiv, r19_comp_usnupp, r19_supphiv)
    )
    dtable[, paste0()]
    prettify_labels(dtable)
    cols <- c("LOC_LAB","SEX_LAB", "AGEGROUP" )
    setcolorder(dtable, cols)
    if(add_asterisks_unaids){
        dtable[, `Age composition of people with unsuppressed HIV` := paste0(`Age composition of people with unsuppressed HIV`, UNAIDS_achieved)]
    }
    dtable[,`:=` (ROUND_LAB=NULL, LOC=NULL, SEX=NULL,ELIGIBLE = NULL, ROUND = NULL, UNAIDS_achieved=NULL)]  

    dtable_plot <- delete.repeated.table.values(dtable, cols)
    dtable_plot 
}
