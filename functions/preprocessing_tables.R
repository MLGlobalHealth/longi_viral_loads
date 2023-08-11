make.table.firstparticipant.NPhiv.NPunsupp <- function(DT=dall){

    tmp <- subset(dall, FIRST_PARTICIPATION == 1)
    tmp1 <- tmp[, j=list( 
        N_HIV=sum(HIV_STATUS, na.rm=TRUE),
        N_UNSUPP=sum(VL_COPIES >= VIREMIC_VIRAL_LOAD, na.rm=TRUE)
    ), by=c('FC', 'SEX', "ROUND") ]
    tmp1[, `:=` ( 
        P_HIV_REDUCTION= ( (1-(N_HIV/N_HIV[ROUND==16])) * 100) |> prettify_cell(percent=TRUE),
        P_UNSUPP_REDUCTION=( (1-(N_UNSUPP/N_UNSUPP[ROUND==16])) * 100) |> prettify_cell(percent=TRUE)
    ), by=c('FC', 'SEX') ]
    tmp1 <- prettify_labels(tmp1) |> 
        remove.nonpretty() |>
        setkey(FC_LAB, SEX_LAB, ROUND_LAB)|>
        setcolorder(c("FC_LAB", "SEX_LAB", "ROUND_LAB")) |>
        delete.repeated.table.values() |>
        setnamesdict(dict=my_labs_dictionary) |>
        setnamesdict(dict=dict_table_names$percent_reduction)
    
    return(tmp1)
}
