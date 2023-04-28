load.summarised.draws.from.files <- function(lab, files, include.raw)
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
