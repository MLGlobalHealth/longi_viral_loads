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
    
    skeleton <- xtable::xtable(DT[0], type='latex') %>% print(hline.after=c(), comment=FALSE) 

    cmd_input <- sprintf("\n \\\\input{%s} \n \\\\end{tabular}", file.path("TODO", basename(file)))
    cmd <- sub('\\n \\\\end\\{tabular\\}', cmd_input, skeleton) 
    cat('\n', cmd, '\n')
    cmd
}
