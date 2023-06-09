###
# Helpers for Latex tables 
### 

delete.repeated.table.values <- function(DT, cols){

    .repeated.to.empty <- function(x){
        if( !is.character(x))
            return(x)
        x[which(x==shift(x))] <- ""
        return(x)
    }
    table <- copy(DT)
    table[, (cols) := lapply(.SD, .repeated.to.empty), .SDcols = cols ]
    return(table)
}

remove.ILIU <- function(DT){
    DT[, `:=` (IL = NULL, IU=NULL)]
} 

table.na.to.empty <- function(DT, cols){

    .na2empty <- function(x){
        if( !is.character(x))
            return(x)
        x[is.na(x)] <- ""
        return(x)
    }
    table <- copy(DT)
    table[, (cols) := lapply(.SD, .na2empty), .SDcols = cols ]
    return(table)
}

table.to.plot <- function(DT){
    gridExtra::tableGrob(DT) |> gridExtra::grid.arrange()
}

prettify_cell <- function(..., parenthesis="(", precision=2, newline=FALSE, percent=FALSE, fmt_skeleton = NA_character_){

    # writes M and CI expression in nice format ready to be put in a table.

    n_vars <- enquos(...) |> length()

    percent_sub <- fifelse(percent, yes='%%', no='') 
    line_sub <- fifelse(newline, yes='\n', no=' ')
    parenthesis_closed <- fcase(
        parenthesis == "(", ")",
        parenthesis == "[", "]",
        parenthesis == "{", "}"
    )

    if(is.na(fmt_skeleton) ){
        fmt_skeleton <- fcase(
            n_vars == 1, "%._prec_fPERCENT",
            n_vars == 2, "(%._prec_f - %._prec_f)",
            n_vars == 3, "%._prec_fPERCENT_line_(%._prec_f - %._prec_f)"
        )
    }
    fmt_skeleton <- gsub('_prec_',precision,fmt_skeleton)
    fmt_skeleton <- gsub('PERCENT', percent_sub, fmt_skeleton)
    fmt_skeleton <- gsub("\\(", parenthesis, fmt_skeleton)
    fmt_skeleton <- gsub("\\)", parenthesis_closed, fmt_skeleton)
    fmt_skeleton <- gsub("_line_", line_sub, fmt_skeleton)

    sprintf(fmt_skeleton,...)
}


write.to.tex <- function(DT, file){

    # write to file the "body" of the table
    tab_latex <- xtable::xtable(DT, type='latex') 
    print(tab_latex, 
        file=file,
        hline.after=c(),
        include.rownames=FALSE,
        include.colnames=FALSE,
        only.contents=TRUE,
        comment=FALSE)

    # however, we want to be able to quickly copy-paste the tabular env
    # so print a skeleton
    
    skeleton <- xtable::xtable(DT[0], type='latex') |>
        print(hline.after=c(), comment=FALSE) 

    cmd_input <- sprintf("\n \\\\input{%s} \n \\\\end{tabular}", file.path("TODO", basename(file)))
    cmd <- sub('\\n \\\\end\\{tabular\\}', cmd_input, skeleton) 
    cat('\n', cmd, '\n')
    cmd
}
