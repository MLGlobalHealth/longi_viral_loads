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

make.main.table.contributions <- function(DTPOP=ncen, DJOINT=djoint_agegroup, DCONTRIB=dcontrib_agegroup, DSUPP=dsupp_agegroup, add_asterisks_unaids=TRUE, include_totals=c("SEX", "AGEGROUP")){

    include_totals <- intersect(include_totals, c("SEX", "AGEGROUP"))
    if(all(include_totals == "SEX")){warning("SEX totals can only be included if AGEGROUP is")}

    # get info for census eligible pop size
    DTPOP[, AGEGROUP:=split.agegroup(AGEYRS)]
    ncen_agegroup <- DTPOP[ ROUND == 19, list(ELIGIBLE=sum(ELIGIBLE)), by=c('ROUND', 'FC', 'SEX','AGEGROUP')] |> setnames("FC", "LOC")
    # get totals
    ncen_totals <- groupingsets( ncen_agegroup, 
        j = list(ELIGIBLE = sum(ELIGIBLE)),
        sets = list("LOC", c("LOC", "SEX")),
        by = c("LOC", "SEX", "AGEGROUP"),
    )[, lapply(.SD,na2string, str="Total")]
    ncen_totals[, ROUND := 19]
    # compute nice cells N(P%)
    expr_eligible_cell <- expr({ sprintf( "%s (%.2f%%)", ELIGIBLE, 100*proportions(ELIGIBLE)) })
    cols <- c("ROUND", "LOC")
    ncen_agegroup[ SEX != "Total" & AGEGROUP != "Total", ELIGIBLE_CELL := eval(expr_eligible_cell) , by=cols]
    ncen_totals[ SEX != "Total" & AGEGROUP == "Total", ELIGIBLE_CELL := eval(expr_eligible_cell), by=c(cols)]
    ncen_totals[ SEX == "Total" & AGEGROUP == "Total", ELIGIBLE_CELL := eval(expr_eligible_cell), by=c(cols)]
    # merge together
    ncen_agegroup <- rbind(ncen_agegroup, ncen_totals, use.names=TRUE)

    e_sex <- expr(ROUND == 19 & SEX != "Total")
    if("SEX" %in% include_totals){
        e_sex <- expr(ROUND == 19)
    }
    e_age <- expr(AGEGROUP != "Total")
    if("AGEGROUP" %in% include_totals){
        e_age <- expr(TRUE)
    }
    
    # get info from gp fits
    r19_comp_usnupp <- DCONTRIB[  MODEL=='run-gp-supp-pop' & eval(e_sex) & eval(e_age)]
    r19_comp_hiv <- DCONTRIB[  MODEL=='run-gp-prevl' & eval(e_sex) & eval(e_age)]
    r19_hivprev <- DJOINT[  MODEL == 'run-gp-prevl' & eval(e_sex)& eval(e_age)] 
    r19_prop_unsupp <- DJOINT[  MODEL == "run-gp-supp-pop" & eval(e_sex)& eval(e_age)]
    r19_supphiv <- DSUPP[ eval(e_sex) & eval(e_age) ] 
    ncen_agegroup <- ncen_agegroup[ eval(e_sex) & eval(e_age) ] 
    
    if(add_asterisks_unaids){
        r19_supphiv[, `:=` ( UNAIDS_achieved = fifelse(M >  .95^3, yes=" * ", no=""), M=NULL, CL=NULL, CU=NULL, IL=NULL, IU=NULL)] 
    }else{NULL}

    if(length(include_totals) == 0){
        check <- lapply(
            list(r19_comp_usnupp, r19_comp_hiv),
            function(DT) {
                DT[AGEGROUP != "Total" & SEX != "Total", sum(M) %between% c(.95, 1.05) |> stopifnot(), by=c("ROUND", "LOC")]
            })
    }
    null <- lapply(
        list(r19_hivprev, r19_prop_unsupp, r19_comp_usnupp, r19_comp_hiv),
        function(DT) {
            DT[, CELL := prettify_cell(M*100, CL*100, CU*100, percent=TRUE)]
            DT[, `:=` (MODEL=NULL,M=NULL, CL=NULL, CU=NULL, IU=NULL,IL=NULL)]
        }
    ); rm(null)
    # prettify all
    names_comp <- c(
        hiv = "Age composition of people with HIV",
        unsupp = "Age composition of people with unsuppressed HIV"
    )
    setnames(ncen_agegroup, "ELIGIBLE_CELL", "Population in age band")
    setnames(r19_hivprev, "CELL", "People with HIV in age band")
    setnames(r19_prop_unsupp, "CELL", "People with unsuppressed HIV in age  band")
    setnames(r19_comp_hiv, "CELL", names_comp['hiv'])
    setnames(r19_comp_usnupp, "CELL", names_comp['unsupp'])
    dtable <- Reduce(
        f=function(x,y) merge(x,y, all=TRUE, by=c("ROUND", "LOC", "SEX", "AGEGROUP")),
        x=list(ncen_agegroup, r19_hivprev, r19_prop_unsupp, r19_comp_hiv, r19_comp_usnupp, r19_supphiv)
    )
    if('SEX' %in% include_totals ){
        dtable[SEX == "Total", (names_comp) := lapply(.SD, function(x) "100.00%"), .SDcols=names_comp]
    }
    prettify_labels(dtable)
    cols <- c("LOC_LAB","SEX_LAB", "AGEGROUP" )
    setcolorder(dtable, cols)
    if(add_asterisks_unaids){
        dtable[, (names_comp['unsupp']) := lapply(.SD, \(x) paste0(x, UNAIDS_achieved)), .SDcols=names_comp['unsupp']]
    }
    dtable[,`:=` (ROUND_LAB=NULL, LOC=NULL, SEX=NULL,ELIGIBLE = NULL, ROUND = NULL, UNAIDS_achieved=NULL)]  

    dtable_plot <- delete.repeated.table.values(dtable, cols)
    dtable_plot 
}
