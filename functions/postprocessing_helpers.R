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

get.weighted.average.p_predict <- function(fit1, fit2, round)
{
    round <- unique(round)
    stopifnot(length(round) == 1)
    warning("Need to check which is FTP. prefertably fit2")

    params <- fit@model_pars %which.like% "^p_predict_[0-9]{2}"
    re1 <- rstan::extract(fit1, pars=params) 
    re2 <- rstan::extract(fit2, pars=params) 
    n_ages <- re1[[1]] |> ncol()
    n_iters <- with(fit2@stan_args[[1]], iter - warmup)
    n_iters <- min(n_iters, with(fit1@stan_args[[1]], iter - warmup))
    
    # check colnames are identical (probably useless)
    (names(re1) == names(re2)) |> all() |> stopifnot()

    # now, get age-sex-loc dependented weighted averages...
    get.quantiles.joint <- function(par){

        sex.idx = gsub("^.*([0-1])[0-1]$", "\\1", par) |> as.integer()
        loc.idx = gsub("^.*[0-1]([0-1])$", "\\1", par) |> as.integer()

        # get weights for age group
        dweights <- dfirst_prop |> 
            subset(
                stan_dicts$SEXtoINT[SEX] == sex.idx & 
                stan_dicts$LOCtoINT[FC]==loc.idx &
                ROUND == round,
                select=c('AGEYRS', 'P')
            ) 
        n_weights <- nrow(dweights)
        weights_ftp <- dweights$P

        # account for 'half-ages' in gp
        if( n_weights == (n_ages - 1)/2 ){
            tmp <- rep(NA, n_ages)
            tmp[2*(1:n_weights)] <- weights_ftp
            weights_ftp <- tmp
        }

        # get weighted average (to check: which is FTP) and quantiles
        re_out <- `+`(
            apply(re1[[par]], 1, `*`, 1 - weights_ftp) |> t(),
            apply(re2[[par]], 1, `*`, weights_ftp) |> t()
        )
        re_out <- re_out[,!is.na(weights_ftp)]
        colnames(re_out) <- paste0("A_", dweights$AGEYRS)
        re_out <- apply(re_out, 2, quantile2) |> 
            t() |> 
            as.data.table(keep.rownames = "AGE_LABEL")
        re_out[, `:=` ( 
            AGE_LABEL=sub('A_', '', AGE_LABEL) |> as.integer(),
            SEX=sex.idx,
            LOC=loc.idx)]
        re_out
    }

    lapply(params, get.quantiles.joint) |> rbindlist()
}

plot.fit.weighted.by.ftpstatus  <- function(DT, label)
{
    dplot <- subset(DT, MODEL == label) |> 
        prettify_labels()

    .f <- function(lab, intercept){
        if(label == lab)
            geom_hline(yintercept=intercept, linetype='dashed', color='grey80')
    }

    ggplot(dplot, aes(x=AGE_LABEL, y=M, ymin=CL, ymax=CU, fill=SEX_LAB, color=SEX_LAB)) +
        .f(lab = 'run-gp-supp-hiv', intercept=.95^3) +
        .f(lab = 'run-gp-supp-hiv', intercept=.95^2) +
        geom_ribbon(alpha=.2,color=NA) +
        geom_line() +
        facet_grid(ROUND_LAB ~ LOC_LAB) +
        scale_color_manual(values=palettes$sex) +
        scale_fill_manual(values=palettes$sex) +
        scale_y_percentage +
        scale_x_continuous(expand=c(0,0)) +
        theme_default() +
        my_labs(y=model_dict[label]) +
        NULL
}

plot.estimated.number.viraemic.among.census.eligible <- function(DT)
{
    ggplot(data=DT, aes(x=AGE_LABEL, y=M, color=SEX_LAB, fill=SEX_LAB)) + 
        geom_col(position=position_dodge()) +
        facet_grid(ROUND_LAB ~ LOC_LAB, scales='free_y')  +
        theme_default() +
        scale_color_manual(values=palettes$sex) + 
        scale_fill_manual(values=palettes$sex) + 
        scale_y_continuous(expand=expansion(mult=c(0, .1))) + 
        my_labs(y="Estimated number of viraemic individuals among census eligible") +
        NULL
}


plot.estimated.contribution.viraemic.among.census.eligible <- function(DT)
{
    ggplot(data=DT, aes(x=AGE_LABEL, y=P, color=SEX_LAB, fill=SEX_LAB)) + 
        geom_line() +
        geom_point() +
        facet_grid(ROUND_LAB ~ LOC_LAB)  +
        theme_default() +
        scale_color_manual(values=palettes$sex) + 
        scale_fill_manual(values=palettes$sex) + 
        scale_y_percentage + 
        my_labs(y="Estimated proportion of viraemic individuals among census eligible") +
        NULL

}
