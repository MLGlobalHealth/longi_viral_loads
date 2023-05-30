quantile2 <- function(x, ps = c(CL = .025, IL = .25, M = .5, IU = .75, CU = .975)) {
    out <- quantile(x = x, probs = ps, names = TRUE) |> t()
    if (!is.null(names(ps))) {
        colnames(out) <- names(ps)
    }
    as.data.table(out)
}

stanindices2vars <- function(names) {
    .gs <- function(reg) {
        out <- gsub(reg, "\\1", names)
        return(out)
    }
    index_sex <- .gs("^.*_([0-1])[0-1].*$")
    index_loc <- .gs("^.*_[0-1]([0-1]).*$")
    index_age <- .gs("^.*\\[([0-9]+)\\].*$") |> as.integer()

    out <- with(stan_dicts, list(
        SEX = unname(INTtoSEX[index_sex]),
        LOC = unname(INTtoLOC[index_loc])
    ))

    cond <- uniqueN(index_age) %between% c(71, 73) | all(is.na(index_age))
    stopifnot("stanindices2vars expects 71 or 73 age groups" = cond)

    if (!all(is.na(index_age))) {
        out$AGE <- 14 + index_age / 2
    }
    out
}


load.summarised.draws.from.rdafiles <- function(lab, files, include.raw) {
    # lab <- 'prevalence' ; files <- rda_files; include.raw <- FALSE;

    lab2 <- fcase(
        lab %like% "prevalence", "prev.hiv.by.age",
        lab %like% "suppAmongInfected", "nsinf.by.age",
        lab %like% "suppAmongPop", "nspop.by.age",
        default = NA
    )
    stopifnot(!is.na(lab2))
    dfiles <- grep(lab, files, value = TRUE)

    # dfit contains the posterior
    dfit <- list()
    draw <- list()
    for (file in dfiles)
    {
        cat("Loading file:", basename(file), "...\n")
        tmp_env <- new.env()
        load(file, envir = tmp_env)
        # rlang::env_has(tmp_env, 'DT')
        dfit[[file]] <- tmp_env[[lab2]]
        ls(tmp_env)
        if (include.raw) {
            draw[[file]] <- tmp_env[["DT"]]
        }
    }
    rm(tmp_env)

    .f <- function(x) as.integer(gsub("^.*?round([0-9]{2}).*?$", "\\1", x))

    dfit <- rbindlist(dfit, idcol = "ROUND")
    dfit[, ROUND := .f(ROUND)]

    # Fix age groups
    dfit <- dfit[AGE_LABEL %% 1 == .5]
    dfit[, AGEYRS := AGE_LABEL - .5]
    # dfit[ , AGE_LABEL := NULL]

    if (include.raw) {
        draw <- rbindlist(draw)
        stopifnot("ROUND" %in% names(draw))
        return(list(raw = draw, fit = dfit))
    }

    return(dfit)
}

store.rda.environments.in.list.by.round.ftpstatus <- function(DFILES) {
    env_list <- list()
    with(DFILES, {
        rounds <- unique(ROUND)
        ftp_bools <- unique(FTP)

        for (round in rounds) {
            env_for_one_round <- list()

            for (ftp_bool in ftp_bools) {
                ftp_label <- fifelse(ftp_bool == TRUE, yes = "ftp", no = "allp")

                e <- new.env()
                idx <- (ROUND == round & FTP == ftp_bool) |> which()
                stopifnot(uniqueN(idx) == uniqueN(MODEL))

                for (i in idx) {
                    cat("Loading", F[i], "...\n")
                    rda.path <- file.path(D[i], F[i])
                    load(rda.path, envir = e)
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

rbind.ftpstatus.datatables.by.round <- function(DTname, round, envir_list = env_list) {
    if (is.numeric(round)) {
        round <- as.character(round)
    }

    DT.allp <- get(DTname, envir = envir_list[[round]][["allp"]])
    DT.ftp  <- get(DTname, envir = envir_list[[round]][["ftp"]])

    DT.allp[, FTP_LAB := "All participants"]
    DT.ftp[, FTP_LAB := "First-time participants"]

    rbind(DT.allp, DT.ftp)
}

plot.comparison.ftptype.colftp <- function(DT, ylab) {
    dplot <- copy(DT)
    prettify_labels(dplot)

    p <- ggplot(dplot, aes(x = AGE_LABEL, y = M, ymin = CL, ymax = CU, color = FTP_LAB, fill = FTP_LAB)) +
        geom_ribbon(alpha = .2, color = NA) +
        geom_line() +
        facet_grid(LOC_LAB ~ SEX_LAB, scales = "free_y") +
        scale_y_continuous(labels = scales::label_percent(), limits = c(0, NA), expand = c(0, 0)) +
        scale_x_continuous(limits = c(15, 50), expand = c(0, 0)) +
        scale_color_manual(values = palettes$rakailogo) +
        scale_fill_manual(values = palettes$rakailogo) +
        theme_default() +
        my_labs(y = ylab) +
        NULL
    return(p)
}


plot.comparison.ftptype.colsex <- function(DT, ylab) {
    dplot <- copy(DT)
    prettify_labels(dplot)

    p <- ggplot(dplot, aes(x = AGE_LABEL, y = M, ymin = CL, ymax = CU, color = SEX_LAB, fill = SEX_LAB)) +
        geom_ribbon(alpha = .2, color = NA) +
        geom_line() +
        facet_grid(LOC_LAB ~ FTP_LAB, scales = "free_y") +
        scale_y_continuous(labels = scales::label_percent(), limits = c(0, NA), expand = c(0, 0)) +
        scale_x_continuous(limits = c(15, 50), expand = c(0, 0)) +
        scale_color_manual(values = palettes$sex) +
        scale_fill_manual(values = palettes$sex) +
        theme_default() +
        my_labs(y = ylab) +
        NULL
    return(p)
}

get.weighted.average.p_predict <- function(fit1, fit2, round, expression_prereturn={}) {
    # round <- 19
    dot.cols <- c(".chain", ".iteration", ".draw") demo.cols <- c("SEX", "LOC", "AGEYRS")

    # extract age-specific draws from both ftp and allparts
    draws_parts <- posterior::as_draws_df(fit1) |>
        posterior::subset_draws(variable = "^p_predict", regex = TRUE)
    names(draws_parts) <- gsub("p_predict", "parts", names(draws_parts))

    draws_ftp <- posterior::as_draws_df(fit2) |>
        posterior::subset_draws(variable = "^p_predict", regex = TRUE)
    names(draws_ftp) <- gsub("p_predict", "ftp", names(draws_ftp))

    # merge together
    draws_all <- posterior::bind_draws(draws_parts, draws_ftp) |> setDT()
    draws_all <- melt(draws_all,
        id.vars = dot.cols,
        variable.name = "participation",
        value.name = "value"
    )
    draws_all[, (demo.cols) := stanindices2vars(participation)]

    draws_all <- merge(
        draws_all[participation %like% "parts", list(.chain, .iteration, .draw, SEX, LOC, AGEYRS, parts = value)],
        draws_all[participation %like% "ftp", list(.chain, .iteration, .draw, SEX, LOC, AGEYRS, ftp = value)],
        by = c(dot.cols, demo.cols),
        all.x = TRUE, all.y = TRUE
    )

    # weight by participation rates to obtain estimates across full pop
    draws_all <- merge(draws_all, dpartrates[ROUND == round], by = c("LOC", "SEX", "AGEYRS"))
    draws_all[, joint := parts * PARTRATE + ftp * (1 - PARTRATE)]

    # evaluate postprocessing expression: 'enriching' is necessary before evaluation
    expression_prereturn <- enexpr(expression_prereturn)

    if(is.null(expression_prereturn))
        return(draws_all[, quantile2(joint), by = c("SEX", "LOC", "AGEYRS")])

    eval(expression_prereturn)
}

get.weighted.average.p_predict.deprecated <- function(fit1, fit2, round) {
    round <- unique(round)
    stopifnot(length(round) == 1)
    warning("Need to check which is FTP. prefertably fit2")

    params <- fit1@model_pars %which.like% "^p_predict_[0-9]{2}"
    re1 <- rstan::extract(fit1, pars = params)
    re2 <- rstan::extract(fit2, pars = params)
    n_ages <- re1[[1]] |> ncol()
    n_iters <- with(fit2@stan_args[[1]], iter - warmup)
    n_iters <- min(n_iters, with(fit1@stan_args[[1]], iter - warmup))

    # check colnames are identical (probably useless)
    (names(re1) == names(re2)) |>
        all() |>
        stopifnot()

    # now, get age-sex-loc dependented weighted averages...
    get.quantiles.joint <- function(par) {
        sex.idx <- gsub("^.*([0-1])[0-1]$", "\\1", par) |> as.integer()
        loc.idx <- gsub("^.*[0-1]([0-1])$", "\\1", par) |> as.integer()

        # get weights for age group
        dweights <- dfirst_prop |>
            subset(
                stan_dicts$SEXtoINT[SEX] == sex.idx &
                    stan_dicts$LOCtoINT[FC] == loc.idx &
                    ROUND == round,
                select = c("AGEYRS", "P")
            )
        n_weights <- nrow(dweights)
        weights_ftp <- dweights$P

        # account for 'half-ages' in gp
        if (n_weights == (n_ages - 1) / 2) {
            tmp <- rep(NA, n_ages)
            tmp[2 * (1:n_weights)] <- weights_ftp
            weights_ftp <- tmp
        }

        # get weighted average (to check: which is FTP) and quantiles
        re_out <- `+`(
            apply(re1[[par]], 1, `*`, 1 - weights_ftp) |> t(),
            apply(re2[[par]], 1, `*`, weights_ftp) |> t()
        )
        re_out <- re_out[, !is.na(weights_ftp)]
        colnames(re_out) <- paste0("A_", dweights$AGEYRS)
        re_out <- apply(re_out, 2, quantile2) |>
            t() |>
            as.data.table(keep.rownames = "AGE_LABEL")
        re_out[, `:=`(
            AGE_LABEL = sub("A_", "", AGE_LABEL) |> as.integer(),
            SEX = sex.idx,
            LOC = loc.idx
        )]
        re_out
    }

    lapply(params, get.quantiles.joint) |> rbindlist()
}

plot.fit.weighted.by.ftpstatus <- function(DT, label) {
    dplot <- subset(DT, MODEL == label) |>
        prettify_labels()

    .h <- function(lab, intercept) {
        if (label == lab) {
            geom_hline(yintercept = intercept, linetype = "dashed", color = "grey80")
        }
    }

    ggplot(dplot, aes(x = AGEYRS, y = M, ymin = CL, ymax = CU, fill = SEX_LAB, color = SEX_LAB)) +
        .h(lab = "run-gp-supp-hiv", intercept = .95^3) +
        .h(lab = "run-gp-supp-hiv", intercept = .95^2) +
        geom_ribbon(alpha = .2, color = NA) +
        geom_line() +
        facet_grid(ROUND_LAB ~ LOC_LAB) +
        scale_color_manual(values = palettes$sex) +
        scale_fill_manual(values = palettes$sex) +
        scale_y_percentage +
        scale_x_continuous(expand = c(0, 0)) +
        theme_default() +
        my_labs(y = model_dict[label]) +
        NULL
}

plot.estimated.number.viraemic.among.census.eligible <- function(DT) {
    ggplot(data = DT, aes(x = AGE_LABEL, y = M, color = SEX_LAB, fill = SEX_LAB)) +
        geom_col(position = position_dodge()) +
        facet_grid(ROUND_LAB ~ LOC_LAB, scales = "free_y") +
        theme_default() +
        scale_color_manual(values = palettes$sex) +
        scale_fill_manual(values = palettes$sex) +
        scale_y_continuous(expand = expansion(mult = c(0, .1))) +
        my_labs(y = "Estimated number of viraemic individuals among census eligible") +
        NULL
}

check_median_contr_approx1 <- function(DT){
    cond <- DT[, abs(sum(M) - 1) , by=c('MODEL', 'ROUND', 'LOC')][, sum(V1) < .5]
    if(! cond)
        warning("Unexpected value for aggregated contributions")
}



plot.estimated.contribution.viraemic.among.census.eligible <- function(DT) 
{
    ggplot(data = DT, aes(x = AGE_LABEL, y = P, color = SEX_LAB, fill = SEX_LAB)) +
        geom_line() +
        geom_point() +
        facet_grid(ROUND_LAB ~ LOC_LAB) +
        theme_default() +
        scale_color_manual(values = palettes$sex) +
        scale_fill_manual(values = palettes$sex) +
        scale_y_percentage +
        my_labs(y = "Estimated proportion of viraemic individuals among census eligible") +
        NULL
}

# from Melodie's: I don't think I can use the function as is but it is pretty nice
find_summary_output <- function(samples, output, vars, transform = NULL, standardised.vars = NULL, names = NULL, operation = NULL, save_output = F)
{
    
    # summarise outputs by period
    ps <- c('M'=0.5, 'CL'=0.025, 'CU'=0.975)

    draws = as.data.table( reshape2::melt(samples[[output]]) )
    if(!is.null(names)){
        setnames(draws, 2:(length(names) + 1), names)
    }else if( draws[, max(Var2)] == df_age[, max(INDEX_AGE)] ){
        setnames(draws, 2:4, c('INDEX_AGE', 'INDEX_DIRECTION', 'INDEX_TIME'))
    }else{
        setnames(draws, 2:4, c('INDEX_DIRECTION', 'INDEX_TIME', 'INDEX_AGE'))
    }

    if('INDEX_AGE' %in% names(draws)){
        draws <- merge(draws, df_age, by = 'INDEX_AGE')
        draws <- merge(draws, df_age_aggregated, by = c('AGE_INFECTION.RECIPIENT', 'AGE_TRANSMISSION.SOURCE'))
    }

    if(!is.null(transform)){
        draws[, value := sapply(value, transform)]
    }

    #  sum force of infection
    if(is.null(operation)){
        draws <- draws[, list(value = sum(value)), by = c('iterations', vars)]
    } else{
        draws <- draws[, list(value = sapply(value, operation)), by = c('iterations', vars)]
    }

    # standardised
    if(!is.null(standardised.vars)){
        draws[, total_value := sum(value), by = c('iterations', standardised.vars)]
        draws[, value := value / total_value]
    }

    #summarise
    draws = draws[, list(q= quantile(value, prob=ps, na.rm = T), q_label=names(ps)), by=vars]	
    draws = dcast(draws, ... ~ q_label, value.var = "q")


    if('INDEX_DIRECTION' %in% vars)
    draws <- merge(draws, df_direction, by = 'INDEX_DIRECTION')

    if('INDEX_TIME' %in% vars){
        draws <- merge(draws, df_period, by = c('INDEX_TIME'))
        draws <- merge(draws, df_community, by = 'COMM')
    }

    if('INDEX_AGE' %in% vars)
    draws <- merge(draws, df_age, by = 'INDEX_AGE')

    if(save_output)
    {
        file = paste0(outdir.table, '-output-', output, 'by_', tolower(paste0(gsub('INDEX_', '', vars), collapse = '_')))
        if(!is.null(standardised.vars)){
            file = paste0(file, 'standardisedby_', tolower(paste0(gsub('INDEX_', '', standardised.vars), collapse = '_')))
        }

        file = paste0(file, '.rds')
        saveRDS(draws, file)
    }

    return(draws)
}

plot.agesex.contributions.by.roundcomm <- function(DT, label) {

    # DT <- copy(dcontrib); label = 'run-gp-supp-pop'
    dplot <- subset(DT, MODEL == label) |> 
        prettify_labels()

    ggplot(dplot, aes(x = AGEYRS, y = M, ymin = CL, ymax = CU, fill = SEX_LAB, color = SEX_LAB)) +
        geom_ribbon(alpha = .2, color = NA) +
        geom_line() +
        facet_grid(ROUND_LAB ~ LOC_LAB) +
        scale_color_manual(values = palettes$sex) +
        scale_fill_manual(values = palettes$sex) +
        scale_y_percentage +
        scale_x_continuous(expand = c(0, 0)) +
        theme_default() +
        my_labs(y = model_dict[label]) +
        NULL
}

plot.agesex.contributions.by.roundcomm2 <- function(DT, label) {

    # DT <- copy(dcontrib); label = 'run-gp-supp-pop'
    dplot <- subset(DT, MODEL == label) |> 
        prettify_labels()

    ggplot(dplot, aes(x = AGEGROUP, y = M, ymin = CL, ymax = CU, fill = SEX_LAB, color = SEX_LAB)) +
        geom_ribbon(alpha = .2, color = NA) +
        geom_col(position='dodge') +
        facet_grid(ROUND_LAB ~ LOC_LAB) +
        scale_color_manual(values = palettes$sex) +
        scale_fill_manual(values = palettes$sex) +
        scale_y_percentage +
        # scale_x_continuous(expand = c(0, 0)) +
        theme_default() +
        my_labs(y = model_dict[label]) +
        NULL
}
