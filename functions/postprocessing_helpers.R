load.summarised.draws.from.rdafiles <- function(lab, files, include.raw)
{
    # lab <- 'prevalence' ; files <- rda_files; include.raw <- FALSE; 

    lab2 <- fcase(
        lab %like% 'prevalence', "prev.hiv.by.age",
        lab %like% 'suppAmongInfected', "nsinf.by.age",
        lab %like% 'suppAmongPop', "nspop.by.age",
        default=NA)
    stopifnot(! is.na(lab2))
    dfiles <- grep( lab, files, value=TRUE)

    # dfit contains the posterior
    dfit <- list()
    draw <- list()
    for (file in dfiles)
    {
        cat('Loading file:', basename(file), '...\n')
        tmp_env <- new.env()
        load(file, envir=tmp_env)
        # rlang::env_has(tmp_env, 'DT')
        dfit[[file]] <- tmp_env[[lab2]]
        ls(tmp_env)
        if(include.raw)
            draw[[file]] <- tmp_env[["DT"]]
    }
    rm(tmp_env)

    .f <- function(x) as.integer(gsub('^.*?round([0-9]{2}).*?$', '\\1', x))

    dfit <- rbindlist(dfit, idcol='ROUND')
    dfit[, ROUND := .f(ROUND)]

    # Fix age groups
    dfit <- dfit[AGE_LABEL %% 1 == .5]
    dfit[ , AGEYRS := AGE_LABEL - .5 ]
    # dfit[ , AGE_LABEL := NULL]

    if(include.raw) {
        draw <- rbindlist(draw)
        stopifnot("ROUND" %in% names(draw))
        return( list(raw=draw, fit=dfit) )
    }

    return(dfit)
}

store.rda.environments.in.list.by.round.ftpstatus <- function(DFILES)
{
    env_list <- list()
    with(DFILES, {

        rounds <- unique(ROUND)
        ftp_bools <- unique(FTP)

        for( round in rounds ){

            env_for_one_round <- list()

            for( ftp_bool in ftp_bools ){

                ftp_label <- fifelse(ftp_bool == TRUE, yes='ftp',no='allp')

                e <- new.env()
                idx <- (ROUND == round & FTP == ftp_bool) |> which()
                stopifnot(uniqueN(idx) == uniqueN(MODEL))

                for(i in idx){
                    cat("Loading", F[i], "...\n")
                    rda.path <- file.path(D[i], F[i])
                    load(rda.path, envir=e )
                }
                e$re <- NULL 
                e$re2 <- NULL 
                env_for_one_round[[ftp_label]] <- e
            }

            env_list[[as.character(round)]] <- env_for_one_round
        }

        return(env_list)
    })
}

plot.comparison.ftptype.colftp <- function(DT, ylab)
{
    dplot <- copy(DT)
    prettify_labels(dplot)

    p <- ggplot(dplot, aes(x=AGE_LABEL, y=M, ymin=CL, ymax=CU, color=FTP_LAB, fill=FTP_LAB)) +
        geom_ribbon(alpha=.2, color=NA) +
        geom_line() +
        facet_grid(LOC_LAB ~ SEX_LAB) + 
        scale_y_continuous(labels=scales::label_percent(), limits=c(0, NA),expand=c(0,0)) +
        scale_x_continuous(limits=c(15,50), expand=c(0,0)) +
        scale_color_manual(values=palettes$rakailogo) + 
        scale_fill_manual(values=palettes$rakailogo) + 
        theme_default() +
        my_labs(y=ylab) +
        NULL
    return(p)
}


plot.comparison.ftptype.colsex <- function(DT, ylab)
{
    dplot <- copy(DT)
    prettify_labels(dplot)

    p <- ggplot(dplot, aes(x=AGE_LABEL, y=M, ymin=CL, ymax=CU, color=SEX_LAB, fill=SEX_LAB)) +
        geom_ribbon(alpha=.2, color=NA) +
        geom_line() +
        facet_grid(LOC_LAB ~ FTP_LAB) + 
        scale_y_continuous(labels=scales::label_percent(), limits=c(0, NA),expand=c(0,0)) +
        scale_x_continuous(limits=c(15,50), expand=c(0,0)) +
        scale_color_manual(values=palettes$sex) + 
        scale_fill_manual(values=palettes$sex) + 
        theme_default() +
        my_labs(y=ylab) +
        NULL
    return(p)
}
