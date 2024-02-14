.make.stan.data.gp.sharedhyper <- function(DTsd,
                               num.var = NA,
                               den.var = NA,
                               sharehyper = TRUE,
                               alpha_hyper_sd = .75,
                               rho_hyper_lower_bound = 1,
                               rho_hyper_upper_bound = 35 ) {
    stopifnot(length(num.var == 1) & length(den.var == 1))

    tmp <- seq(DTsd[, min(AGE_LABEL)], DTsd[, max(AGE_LABEL) + 1], .5)
    tmp1 <- which(tmp %% 1 == 0.5)

    stan.data <- list(
        x_predict = tmp,
        N_predict = length(tmp),
        observed_idx = tmp1,
        N_observed = length(tmp1),
        alpha_hyper_par_00 = alpha_hyper_sd,
        alpha_hyper_par_10 = alpha_hyper_sd,
        alpha_hyper_par_01 = alpha_hyper_sd,
        alpha_hyper_par_11 = alpha_hyper_sd,
        rho_hyper_lower_bound = rho_hyper_lower_bound,
        rho_hyper_upper_bound = rho_hyper_upper_bound
    )

    # compute "numerators" and "denominators" by participant type, and append to stan.data
    .standata.get.num.and.den <- function(data){
        list(
            y_observed_00 = data[SEX == 0 & LOC == 0, ..num.var][[1]],
            y_observed_10 = data[SEX == 1 & LOC == 0, ..num.var][[1]],
            y_observed_01 = data[SEX == 0 & LOC == 1, ..num.var][[1]],
            y_observed_11 = data[SEX == 1 & LOC == 1, ..num.var][[1]],
            total_observed_00 = data[SEX == 0 & LOC == 0, ..den.var][[1]],
            total_observed_10 = data[SEX == 1 & LOC == 0, ..den.var][[1]],
            total_observed_01 = data[SEX == 0 & LOC == 1, ..den.var][[1]],
            total_observed_11 = data[SEX == 1 & LOC == 1, ..den.var][[1]]
        )
    }

    # if no participant type passed, return stan.data obtained from dataframe
    if( ! "PTYPE" %in% names(DTsd) ) {
        stan.data <- c(stan.data, .standata.get.num.and.den(DTsd))
        return(stan.data)
    }

    # else get stan.data for each participant type
    for(ptype in unique(DTsd$PTYPE)){

        stan.data.type <- .standata.get.num.and.den(DTsd[PTYPE == ptype, ])
        if(ptype == 'ftp' | ptype == 0)
            names(stan.data.type) <- paste0(names(stan.data.type), '_ftp')
        stan.data <- c(stan.data, stan.data.type)
    }

    # some sanity checks
    with(stan.data, {c(
        length(y_observed_00) == N_observed,
        length(y_observed_01) == N_observed,
        length(y_observed_10) == N_observed,
        length(y_observed_11) == N_observed,
        all(y_observed_00_ftp <= y_observed_00),
        all(y_observed_01_ftp <= y_observed_01),
        all(y_observed_10_ftp <= y_observed_10),
        all(y_observed_11_ftp <= y_observed_11),
        all(total_observed_00_ftp <= total_observed_00),
        all(total_observed_01_ftp <= total_observed_01),
        all(total_observed_10_ftp <= total_observed_10),
        all(total_observed_11_ftp <= total_observed_11)
    ) |> all() |> stopifnot()})

    with(stan.data, {
        any(y_observed_00 != y_observed_00_ftp) |> stopifnot()
        any(y_observed_01 != y_observed_01_ftp) |> stopifnot()
        any(y_observed_10 != y_observed_10_ftp) |> stopifnot()
        any(y_observed_11 != y_observed_11_ftp) |> stopifnot()
    })

    return(stan.data)
}

get.draws.wholepop <- function(DRAWS, DPART){

    dpartrates <- copy(DPART)
    draws_wholepop <- dcast.data.table(DRAWS,
        SEX + LOC + SEX_LABEL + LOC_LABEL + AGE_LABEL + .draw ~ PTYPE,
        value.var = "P")
    setnames(dpartrates, 
        c("AGEYRS", "SEX", "LOC"),
        c("AGE_LABEL", "SEX_LABEL", "LOC_LABEL"))
    draws_wholepop <- merge(draws_wholepop, dpartrates, by=c("LOC_LABEL", "SEX_LABEL", "AGE_LABEL"))
    draws_wholepop[, P := PARTRATE*all + (1-PARTRATE)*ftp]
    return(draws_wholepop)
}

plot.quantiles.wholepop <- function(DQUANT, ylim=NA, ylab){
    ggplot(DQUANT, aes(x=AGE_LABEL, y=M, ymin=CL, ymax=CU, color=SEX_LABEL,fill=SEX_LABEL)) + 
        geom_ribbon(alpha=.3, color=NA) +
        geom_line() +
        facet_grid(. ~ LOC_LABEL , labeller=labeller(LOC_LABEL=community_dictionary$long)) + 
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(labels=scales::label_percent(), limits=c(0,ylim), expand=c(0,0)) + 
        scale_color_manual(values=palettes$sex, labels=sex_dictionary2) +
        scale_fill_manual(values=palettes$sex, labels=sex_dictionary2) +
        my_labs(y=ylab) +
        theme_default()
}

vl.prevalence.by.gender.loc.age.gp.cmdstan.hyper <- function(
    DT=vla, 
    refit = FALSE,
    vl.out.dir. = vl.out.dir,
    alpha_hyper = .75
) {

    cat("\n\n--- Analysing HIV+ Prevalence ---\n\n")
    vla <- copy(DT)

    # Stan file locations
    file.stan <- file.path(gitdir.stan, "vl_binomial_gp_sharedhyper.stan")

    .fit.stan.and.plot.by.round <- function(DT) {

        #  DT <- copy(vla[ROUND == 16]); vl.out.dir.=vl.out.dir
        round <- DT[, unique(ROUND)]
        stopifnot(length(round) == 1)

        stan.data <- .make.stan.data.gp.sharedhyper(DT,
            num.var = "HIV_N",
            den.var = "N",
            alpha_hyper_sd = alpha_hyper
        )

        # file paths for fit and standata
        .prefix <- file.path(vl.out.dir., "cmd_hivprevalence_gp_stanfit_round")
        filename <- paste0(.prefix, round, ".rds")
        filename_standata <- paste0(.prefix, round, "_standata.rds")

        sprintf("Fitting stan model for round %s, with filename:\n %s\n",
            round, filename) |> cat()

        if (file.exists(filename) & refit == FALSE) {

            cat("Loading previously run HMC... \n")
            fit <- readRDS(filename)
            stan.data <- readRDS(filename_standata)

        } else {

            stan.model <- cmdstan_model(stan_file = file.stan, cpp_options=list(stan_threads = TRUE))
            stan.args <- yaml::read_yaml(path.stan.config)

            fit <- stan.model$sample(
                data = stan.data,
                seed = stan.args$seed,
                iter_sampling = stan.args$iter,
                iter_warmup = stan.args$warmup,
                chains = stan.args$chains,
                max_treedepth = stan.args$control$max_treedepth,
                adapt_delta = stan.args$control$adapt_delta,
                parallel_chains = stan.args$chains,
                threads_per_chain = 1
            )

            # save
            fit$save_object(file = filename)
            saveRDS(object=stan.data, file = filename_standata)

        }

        catn("Analyse posterior")
        # _______________________

        re <- fit$draws(format="df") |> as.data.table()
        names_vars <- dimnames(re)[[2]]
        group_codes <- unique(DT[, .(SEX, SEX_LABEL, LOC, LOC_LABEL)])

        # quantiles vectors
        ps <- c(CL=0.025, IL=0.25, M=0.5, IU=0.75, CU=0.975)

        # Extract summary excluding hyper parameters
        dsum <- .get.summary.without.hyperparams(FIT=fit, verbose=TRUE, vars=names_vars)

        # extract hyperparams rho
        prev.hiv.gp.pars <- extract.stan.hyperparams.rho(
            re=re, 
            encoding=group_codes
        )
        tmp <- .get.prior.ranges(
            stan.data,
            DT = group_codes,
            shape = re$rho_hyper_par_shape2[1],
            scale = re$rho_hyper_par_scale2[1]
        )
        p <- .plot.gp.hyperparameters(prev.hiv.gp.pars, tmp)
        filename <- paste0("hivprevalence_gppars_round", round, ".pdf")
        ggsave2(p, file = filename, w = 6, h = 3)

        catn("make use stan.data for PPC.")
        # _________________________________

        ppDT <- copy(DT)
        cols <- c("M", "CU", "CL")
        ppDT[, (cols) := binconf(HIV_N, N, return.df = T)]

        catn("make prevalence plot by age")
        # _________________________________

        q <- c("M"=.5, "CL"=.025, "CU"=.975)
        cols <- names(re) %which.like% '^p_predict_'
        tmp <- re[, lapply(.SD, posterior::quantile2, probs=q), .SDcols =cols]
        tmp[, quantile := names(q)]
        tmp <- melt( tmp, id.vars = "quantile")


        prev.hiv.by.age <- dcast.data.table(tmp, variable ~ quantile, value.var = "value")
        prev.hiv.by.age[ , `:=` (
            AGE_LABEL = .stan.brackets.to.age(variable, .stan.data=stan.data),
            variable = .stan.remove.brackets(variable)
        ) ]
        prev.hiv.by.age <- .stan.get.sex.and.loc(prev.hiv.by.age, 'variable', codes = group_codes)

        plots_ftp <- .plot.stan.fit(
            prev.hiv.by.age[PTYPE == "ftp"],
            DT2=ppDT[PTYPE == "ftp"],
            ylims = c(0,.75),
            ylab="HIV prevalence (95% credible intervals)"
        )
        filenames <- paste0("hivprevalence_vs_age_by_gender_fishinland_ftp_",c("", "data_"),"gp_round", round, ".pdf")
        ggsave2(plots_ftp[[1]], file = filenames[[1]], w = 6, h = 5)
        ggsave2(plots_ftp[[2]], file = filenames[[2]], w = 6, h = 5)

        plots_all <- .plot.stan.fit(
            prev.hiv.by.age[PTYPE == "all"],
            DT2=ppDT[PTYPE == "all"],
            ylims = c(0,.75),
            ylab="HIV prevalence (95% credible intervals)"
        )
        filenames <- paste0("hivprevalence_vs_age_by_gender_fishinland_all_",c("", "data_"),"gp_round", round, ".pdf")
        ggsave2(plots_all[[1]], file = filenames[[1]], w = 6, h = 5)
        ggsave2(plots_all[[2]], file = filenames[[2]], w = 6, h = 5)


        plots_comparison <- ggplot(prev.hiv.by.age, aes(x=AGE_LABEL, y=M, ymin=CL, ymax=CU, color=PTYPE, fill=PTYPE)) + 
            geom_ribbon(alpha=.3, color=NA) +
            geom_line() +
            facet_grid(LOC_LABEL ~ SEX_LABEL) +
            scale_fill_manual(values=palettes$ptype) +
            scale_color_manual(values=palettes$ptype) +
            my_labs() 
        filename <- paste0("hivprevalence_vs_age_by_gender_fishinland_ptype_gp_round", round, ".pdf")
        ggsave2(plots_comparison, file = filename, w = 6, h = 9)


        catn("extract basic prevalence estimates")
        # ________________________________________

        # ps <- c(0.025, 0.5, 0.975)
        draws <- subset(re, select=names(re) %like% '^p_predict|\\.draw')
        draws <- melt( draws,  id.vars='.draw', value.name='P')
        draws[ , `:=` (
            AGE_LABEL = .stan.brackets.to.age(variable, .stan.data=stan.data),
            variable = .stan.remove.brackets(variable)
        ) ]
        draws <- .stan.get.sex.and.loc(draws, 'variable', codes=group_codes)
        draws[, variable := NULL ]

        draws_wholepop <- get.draws.wholepop(DRAWS=draws, DPART=dpartrates[ROUND==round])
        quants_wholepop <- draws_wholepop[, quantile2(P), by=c("AGE_LABEL", "SEX_LABEL", "LOC_LABEL")] 
        p <- plot.quantiles.wholepop(quants_wholepop, ylab="HIV prevalence (95% credible intervals)")
        filename <- paste0("hivprevalence_vs_age_by_gender_fishinland_joint_round", round, ".pdf")
        cmd <- ggsave2(p, file = filename, w = 9, h = 6)

        rp <- merge(
            draws, DT[, .(SEX, LOC, PTYPE,AGE_LABEL, N )],
            by=c('SEX', 'LOC', "PTYPE", "AGE_LABEL"))
        rp <- rp[, .(P = sum(P * N) / sum(N)), by = c("LOC", "SEX", "PTYPE", ".draw")]

        prev.hiv.by.sex.loc <- rp[, quantile2(P, ps = q) , by=c("LOC", "SEX", "PTYPE")]
        prev.hiv.by.sex.loc [, LABEL := prettify_cell(M*100, CL*100, CU*100, percent = TRUE) ]
        prev.hiv.by.sex.loc <- merge(prev.hiv.by.sex.loc, group_codes, by=c("LOC", "SEX"))


        catn("extract prevalence ratio female:male and male:female")
        # __________________________________________________________

        rp <- merge(group_codes, rp, by = c("LOC", "SEX"))
        rp <- dcast.data.table(rp, LOC_LABEL + PTYPE + .draw ~ SEX_LABEL, value.var = "P")
        rp[,  `:=` (PR_FM = F / M, PR_MF = M / F, M = NULL, F = NULL)]

        rp <- melt(rp, id.vars = c("LOC_LABEL", "PTYPE",".draw"))
        rp <- rp[, quantile2(value, ps = q), by = c("LOC_LABEL", "PTYPE","variable")]
        rp[, LABEL := prettify_cell(M*100, CL*100, CU*100, percent = TRUE) ]

        prevratio.hiv.by.loc <- copy(rp)

        catn("plot prevalence ratio F:M and M:F by age")
        # ______________________________________________


        rp <- dcast.data.table(draws, 
            LOC + LOC_LABEL +  PTYPE + .draw + AGE_LABEL ~ SEX_LABEL, value.var = "P")
        rp[, `:=` (PR_FM = F / M, PR_MF = M / F, M = NULL, F = NULL)]
        rp <- melt( rp, 
            id.vars = c("LOC_LABEL", "PTYPE", "AGE_LABEL", ".draw"),
            measure.vars = c("PR_FM", "PR_MF"))
        prevratio.hiv.by.loc.age <- rp[, 
            quantile2(value),
            by = c("LOC_LABEL", "PTYPE", "AGE_LABEL", "variable")]

        p <- .plot.stan.fit.ratio(
            prevratio.hiv.by.loc.age[variable == "PR_FM"],
            ylab = "female to male HIV prevalence ratio\n(95% credible interval)\n"
        ) + facet_grid( LOC_LABEL ~ PTYPE, scales="free_y")
        filename <- paste0("fit_hivprevalenceratio_vs_age_by_fishinland_stan_round", round, ".pdf")
        ggsave2(p, file = filename, w = 8, h = 8)

        catn("extract if diff in F:M prevalence risk ratio in fish vs inland")
        # ____________________________________________________________________

        by_cols <- c("LOC", "SEX", "PTYPE", "AGE_LABEL")
        rp <- merge(subset(DT, select = c(by_cols, 'N') ), draws, by = by_cols )
        rp <- rp[, 
            .(P = sum(P * N) / sum(N)),
            by = c("LOC_LABEL","SEX_LABEL", "PTYPE",".draw")]
        rp <- dcast.data.table(rp, 
            LOC_LABEL + PTYPE + .draw ~ SEX_LABEL, value.var = "P")

        rp[, `:=` (PR_FM = F / M, PR_MF = M / F, M = NULL, F = NULL)]

        rp <- dcast.data.table(rp, PTYPE + .draw ~ LOC_LABEL , value.var = "PR_FM")
        rp[, {
            PR_FM_D <- inland - fishing
            quantile2(PR_FM_D, ps = q)
        }, by="PTYPE"]

        # Save diff. qttities
        filename <- paste0("hivprevalence_round", round, ".rda")
        save(DT, stan.data, re, prev.hiv.by.age, prevratio.hiv.by.loc,
            prev.hiv.by.sex.loc, prevratio.hiv.by.loc.age,
            file = file.path(vl.out.dir., filename)
        )

        cat("Round", round, ": done.\n")
        return(TRUE)
    }

    if (parallelise) {
        # Use parallel backend if run locally
        foreach(
            r = vla[, unique(ROUND)],
            .combine = "c"
        ) %dopar% {
            cat("Running Round", r, "\n")
            .fit.stan.and.plot.by.round(vla[ROUND == r, ])
        } -> tmp
        return(tmp)
    } else {
        # else do not
        for (r in args$round) {
            .fit.stan.and.plot.by.round(vla[ROUND == r, ])
        }
    }
}

vl.suppofinfected.by.gender.loc.age.gp.cmdstan.hyper <- function(
    DT=vla,
    refit = FALSE,
    vl.out.dir. = vl.out.dir,
    alpha_hyper = .75
) {

    cat("\n\n--- Analyse suppressed among infected ---\n\n")

    # DT <- copy(dall); refit=FALSE; vl.out.dir.=vl.out.dir
    vla <- copy(DT)

    # Stan file locations
    file.stan <- file.path(gitdir.stan, "vl_binomial_gp_sharedhyper.stan")

    .fit.stan.and.plot.by.round <- function(DT) {

        # DT <- copy(vla[ROUND == 16] )
        round <- DT[, unique(ROUND)]
        stopifnot(length(round) == 1)

        # Exclude people with unknown Viral load from analyses
        DT[, HIV_N := HIV_N - VLNA_N ]
        stopifnot(DT[, all(HIV_N>0)])
        DT[, VLSUP_N := HIV_N - VLNS_N, ]
        stan.data <- .make.stan.data.gp.sharedhyper(DT, 
            num.var = "VLSUP_N",
            den.var = "HIV_N",
            alpha_hyper_sd = alpha_hyper
        )

        .prefix <- file.path(vl.out.dir., "cmd_notsuppAmongInfected_gp_stan_round")
        filename <- paste0(.prefix, round, ".rds")
        filename_standata <- paste0(.prefix, round, "_standata.rds")

        sprintf("Fitting stan model for round %s, with filename:\n %s\n",
            round, filename) |> cat()

        if (file.exists(filename) & refit == FALSE) {

            cat("Loading previously run HMC... \n")
            fit <- readRDS(filename)
            stan.data <- readRDS(filename_standata)

        } else {

            stan.model <- cmdstan_model(stan_file = file.stan, cpp_options=list(stan_threads = TRUE))
            stan.args <- yaml::read_yaml(path.stan.config)

            fit <- stan.model$sample(
                data = stan.data,
                seed = stan.args$seed,
                iter_sampling = stan.args$iter,
                iter_warmup = stan.args$warmup,
                chains = stan.args$chains,
                max_treedepth = stan.args$control$max_treedepth,
                adapt_delta = stan.args$control$adapt_delta,
                parallel_chains = stan.args$chains,
                threads_per_chain = 1
            )

            # save
            fit$save_object(file = filename)
            saveRDS(object=stan.data, file = filename_standata)
        }

        catn("Analyse posterior")
        # _______________________

        re <- fit$draws(format="df") |> as.data.table()
        names_vars <- dimnames(re)[[2]]
        group_codes <- unique(DT[, .(SEX, SEX_LABEL, LOC, LOC_LABEL)])

        # quantiles vectors
        ps <- c(CL=0.025, IL=0.25, M=0.5, IU=0.75, CU=0.975)

        # Extract summary excluding hyper parameters
        dsum <- .get.summary.without.hyperparams(FIT=fit, verbose=TRUE, vars=names_vars)

        # extract hyperparams rho
        prev.hiv.gp.pars <- extract.stan.hyperparams.rho(
            re=re, 
            encoding=group_codes
        )
        tmp <- .get.prior.ranges(
            stan.data,
            DT = group_codes,
            shape = re$rho_hyper_par_shape2[1],
            scale = re$rho_hyper_par_scale2[1]
        )
        p <- .plot.gp.hyperparameters(prev.hiv.gp.pars, tmp)
        filename <- paste0("cmd_notsuppAmongInfected_gppars_round", round, ".pdf")
        ggsave2(p, file = filename, w = 6, h = 3)

        catn("make use stan.data for PPC.")
        # _________________________________

        ppDT <- copy(DT)
        cols <- c("M", "CU", "CL")
        ppDT[, (cols) := binconf(HIV_N - VLNS_N, HIV_N, return.df = T)]

        catn("make prevalence plot by age")
        # _________________________________
        
        ps <- c("CL"=.025, "IL"=.25, "M"=.5, "IU"=.75, "CU"=.975)
        cols <- names(re) %which.like% '^p_predict_'
        tmp <- re[, lapply(.SD, posterior::quantile2, probs=ps), .SDcols=cols]
        tmp[, quantile := names(ps)]
        tmp <- melt(tmp, id.vars = 'quantile')

        nsinf.by.age <- dcast.data.table(tmp, variable ~ quantile, value.var = "value")
        nsinf.by.age[ , `:=` (
            AGE_LABEL = .stan.brackets.to.age(variable, .stan.data=stan.data),
            variable = .stan.remove.brackets(variable)
        ) ]
        # nsinf.by.age[, sapply(.SD, function(x) any(is.na(x)))]
        nsinf.by.age <- .stan.get.sex.and.loc(nsinf.by.age, 'variable', codes=group_codes)

        plots_ftp <- .plot.stan.fit(
            nsinf.by.age[PTYPE == "ftp"],
            DT2=ppDT[PTYPE == "ftp"],
            ylims = c(0,1),
            ylab = "HIV+ individuals with suppressed viral load\n(95% credible interval)\n"
        )
        filenames <- paste0("suppAmongInfected_vs_age_by_gender_fishinland_ftp_",c("", "data_"),"gp_round", round, ".pdf")
        ggsave2(plots_ftp[[1]], file = filenames[[1]], w = 6, h = 5)
        ggsave2(plots_ftp[[2]], file = filenames[[2]], w = 6, h = 5)

        plots_all <- .plot.stan.fit(
            nsinf.by.age[PTYPE == "all"],
            DT2=ppDT[PTYPE == "all"],
            ylims = c(0,1),
            ylab = "HIV+ individuals with suppressed viral load\n(95% credible interval)\n"
        )
        filenames <- paste0("suppAmongInfected_vs_age_by_gender_fishinland_all_",c("", "data_"),"gp_round", round, ".pdf")
        ggsave2(plots_all[[1]], file = filenames[[1]], w = 6, h = 5)
        ggsave2(plots_all[[2]], file = filenames[[2]], w = 6, h = 5)


        catn("extract basic not supp estimates")
        # ______________________________________

        q <- c("CL"=0.025, "M"=0.5, "CU"=0.975)

        draws <- subset(re, select=names(re) %like% '^p_predict|\\.draw')
        draws <- melt( draws,  id.vars='.draw', value.name='P')
        draws[ , `:=` (
            AGE_LABEL = .stan.brackets.to.age(variable, .stan.data=stan.data),
            variable = .stan.remove.brackets(variable)
        ) ]
        draws <- .stan.get.sex.and.loc(draws, 'variable', codes=group_codes)
        draws[, variable := NULL ]

        draws_wholepop <- get.draws.wholepop(DRAWS=draws, DPART=dpartrates[ROUND==round])
        quants_wholepop <- draws_wholepop[, quantile2(P), by=c("AGE_LABEL", "SEX_LABEL", "LOC_LABEL")] 
        p <- plot.quantiles.wholepop(quants_wholepop, ylim=1,ylab="HIV+ individuals with suppressed viral load\n(95% credible interval)\n")
        filename <- paste0("suppAmongInfected_vs_age_by_gender_fishinland_joint_round", round, ".pdf")
        cmd <- ggsave2(p, file = filename, w = 9, h = 6)

        rp <- merge(
            draws, DT[, .(SEX, LOC, AGE_LABEL, PTYPE ,N )],
            by=c('SEX', 'LOC', "PTYPE","AGE_LABEL"))
        rp <- rp[, .(P = sum(P * N) / sum(N)), by = c("LOC", "SEX", "PTYPE",".draw")]

        nsinf.by.sex.loc <- rp[, quantile2(P, ps = q) , by=c("LOC", "SEX", "PTYPE")]
        nsinf.by.sex.loc [, LABEL := prettify_cell(M*100, CL*100, CU*100, percent = TRUE) ]
        nsinf.by.sex.loc <- merge(nsinf.by.sex.loc, group_codes, by=c('LOC', 'SEX'))


        catn("extract risk ratio of suppressed VL female:male and male:female")
        # _____________________________________________________________________

        rp <- merge(
            draws, DT[, .(SEX_LABEL, LOC_LABEL, AGE_LABEL, N )],
            by=c('SEX_LABEL', 'LOC_LABEL', "AGE_LABEL"))
        rp <- rp[, .(P = sum(P * N) / sum(N)), by = c("LOC_LABEL", "SEX_LABEL", ".draw")]

        rp <- dcast.data.table(rp, LOC_LABEL + .draw ~ SEX_LABEL, value.var = "P")
        rp <- rp[, .(PR_FM = F / M, PR_MF = M / F), by = c("LOC_LABEL", ".draw")]
        rp <- melt(rp, id.vars = c("LOC_LABEL", ".draw"))
        rp <- rp[, quantile2(value, ps = q), by = c("LOC_LABEL", "variable")]
        rp[, LABEL := .write.CIs(M, CL, CU, d = 2)]
        nsinf.by.loc <- copy(rp)

        catn("extract risk ratio of unsuppressed VL female:male and male:female by age")
        # ______________________________________________________________________________

        rp <- dcast.data.table(draws, 
            LOC_LABEL + .draw + AGE_LABEL ~ SEX_LABEL, 
            value.var = "P")
        rp[, `:=` (PR_FM = F / M, PR_MF = M / F, M = NULL, F = NULL)]
        rp <- melt(rp,
            id.vars = c("LOC_LABEL", "AGE_LABEL", ".draw"),
            measure.vars = c("PR_FM", "PR_MF"))
        rp <- rp[, quantile2(value), by=c("LOC_LABEL", "AGE_LABEL", "variable")]
        rp[, LABEL := .write.CIs(M, CL, CU, d = 2)]
        nsinf.ratio.by.loc.age <- copy(rp)


        catn("extract if F:M risk ratio of unsupp VL is diff in fishing vs inland")
        # _________________________________________________________________________

        rp <- merge(
            draws, DT[, .(SEX_LABEL, LOC_LABEL, AGE_LABEL, N )],
            by=c('SEX_LABEL', 'LOC_LABEL', "AGE_LABEL"))
        rp <- rp[, .(P = sum(P * N) / sum(N)), by = c("LOC_LABEL", "SEX_LABEL", ".draw")]
        rp <- dcast.data.table(rp, 
            LOC_LABEL + .draw ~ SEX_LABEL, value.var = "P")
        rp[, `:=` (PR_FM = F / M, PR_MF = M / F, M = NULL, F = NULL)]
        rp <- dcast.data.table(rp, .draw ~ LOC_LABEL, value.var = "PR_FM")
        rp[, PR_FM_D := fishing - inland]
        rp <- rp[, quantile2(PR_FM_D, ps)]

        filename <- paste0("suppAmongInfected_round", round, ".rda")
        save(DT, stan.data, re, 
            nsinf.by.age, nsinf.by.sex.loc,
            nsinf.by.loc, nsinf.ratio.by.loc.age,
            file = file.path(vl.out.dir., filename)
        )


        catn("make table version suppressed")
        # ___________________________________

        nsinf.by.age[, LABEL := prettify_cell(M*100, CL*100,CU*100, percent = TRUE, precision=1)]
        nsinf.ratio.by.loc.age[, LABEL2 := prettify_cell(M, CL, CU)]

        vec_ages <- c(20.5, 25.5, 30.5, 35.5, 40.5, 45.5)

        dt <- subset(nsinf.by.age, AGE_LABEL %in% vec_ages)
        dt <- dcast.data.table(dt, LOC_LABEL + AGE_LABEL ~ SEX_LABEL, value.var = "LABEL")
        .f <- function(var){
            tmp <- subset( nsinf.ratio.by.loc.age,
                variable == var & AGE_LABEL %in% vec_ages,
                select=c('LOC_LABEL', 'AGE_LABEL', 'LABEL2')
            )
            tmp <- setnames(tmp, "LABEL2", var)
            tmp
        }
        tmp_fm <- .f("PR_FM"); tmp_mf <- .f("PR_MF")
        dt <- merge(dt, tmp_fm, by = c("LOC_LABEL", "AGE_LABEL"))
        dt <- merge(dt, tmp_mf, by = c("LOC_LABEL", "AGE_LABEL"))

        filename <- paste0("tab_suppamonginfected_round", round, ".csv")
        fwrite(dt, row.names = FALSE, file = file.path(vl.out.dir., filename))

        TRUE
    }

    if (parallelise) {
        # Use parallel backend if run locally
        foreach(
            r = vla[, unique(ROUND)],
            .combine = "c"
        ) %dopar% {
            cat("Running Round", r, "\n")
            .fit.stan.and.plot.by.round(vla[ROUND == r, ])
        } -> tmp
        return(tmp)
    } else {
        # else do not
        for (r in args$round) {
            .fit.stan.and.plot.by.round(vla[ROUND == r, ])
        }
    }
}


vl.suppofpop.by.gender.loc.age.gp.cmdstan.hyper <- function(
    DT=vla, 
    refit = FALSE,
    vl.out.dir. = vl.out.dir,
    alpha_hyper = .75
)
{

    cat("\n\n--- Analyse suppression among participants ---\n\n")
    vla <- copy(DT)

    # Stan file locations
    file.stan <- file.path(gitdir.stan, "vl_binomial_gp_sharedhyper.stan")
    .fit.stan.and.plot.by.round <- function(DT) {

        #  DT <- copy(vla[ROUND == 16]); refit=FALSE
        round <- DT[, unique(ROUND)]
        stopifnot(length(round) == 1)

        # Remove people with unknwon Viral load,
        # else equivalent to assuming they are suppressed or negative
        DT[, N := N - VLNA_N]
        stan.data <- .make.stan.data.gp.sharedhyper(DT,
            num.var = "VLNS_N",
            den.var = "N",
            alpha_hyper_sd = alpha_hyper
        )

        # outputs paths
        .prefix <- file.path(vl.out.dir., "cmd_suppAmongPop_gp_stan_round")
        filename <- paste0(.prefix, round, ".rds")
        filename_standata <- paste0(.prefix, round, "_standata.rds")

        sprintf("Fitting stan model for round %s, with filename:\n %s\n",
            round, filename) |> cat()

        if (file.exists(filename) & refit == FALSE) {

            cat("Loading previously run HMC... \n")
            fit <- readRDS(filename)
            stan.data <- readRDS(filename_standata)

        } else {
            
            stan.model <- cmdstan_model(stan_file = file.stan, cpp_options=list(stan_threads = TRUE))
            stan.args <- yaml::read_yaml(path.stan.config)

            fit <- stan.model$sample(
                data = stan.data,
                seed = stan.args$seed,
                iter_sampling = stan.args$iter,
                iter_warmup = stan.args$warmup,
                chains = stan.args$chains,
                max_treedepth = stan.args$control$max_treedepth,
                parallel_chains = stan.args$chains,
                threads_per_chain = 1
            )

            # save
            fit$save_object(file = filename)
            saveRDS(object=stan.data, file = filename_standata)
        }

        catn("Analyse posterior")
        # _______________________

        re <- fit$draws(format="df") |> as.data.table()
        names_vars <- dimnames(re)[[2]]
        group_codes <- unique(DT[, .(SEX, SEX_LABEL, LOC, LOC_LABEL)])

        # quantiles vectors
        ps <- c(CL=0.025, IL=0.25, M=0.5, IU=0.75, CU=0.975)

        # Extract summary excluding hyper parameters
        dsum <- .get.summary.without.hyperparams(FIT=fit, verbose=TRUE, vars=names_vars)

        # extract hyperparams rho
        nspop.gp.pars <- extract.stan.hyperparams.rho(
            re=re,
            encoding=group_codes)
        tmp <- .get.prior.ranges(
            stan.data,
            DT = group_codes,
            shape = re$rho_hyper_par_shape2[1],
            scale = re$rho_hyper_par_scale2[1]
        )
        p <- .plot.gp.hyperparameters(nspop.gp.pars, tmp)
        filename <- paste0("notsuppAmongPop_gppars_round", round, ".pdf")
        ggsave2(p, file = filename, w = 6, h = 3)


        catn("make use stan.data for PPC.")
        # _________________________________

        ppDT <- copy(DT)
        cols <- c("M", "CU", "CL")
        ppDT[, (cols) := binconf(VLNS_N, N, return.df = T)]

        catn("make prevalence plot by age and ptype")
        # _________________________________

        q <- c("M"=.5, "CL"=.025, "CU"=.975)
        cols <- names(re) %which.like% '^p_predict_'
        tmp <- re[, lapply(.SD, posterior::quantile2, probs=q), .SDcols =cols]
        tmp[, quantile := names(q)]
        tmp <- melt( tmp, id.vars = "quantile")

        nspop.byage.ptype <- dcast.data.table(tmp, variable ~ quantile, value.var = "value")
        nspop.byage.ptype[ , `:=` (
            AGE_LABEL = .stan.brackets.to.age(variable, .stan.data=stan.data),
            variable = .stan.remove.brackets(variable)
        ) ]
        nspop.byage.ptype <- .stan.get.sex.and.loc(nspop.byage.ptype, 'variable', codes=group_codes)

        plots_ftp <- .plot.stan.fit(
            nspop.byage.ptype[PTYPE == "ftp"],
            DT2=ppDT[PTYPE == "ftp"],
            ylims = c(0,.4),
            ylab = "population with unsuppressed viral load\n(95% credible interval)\n"
        )

        filenames <- paste0("notsuppAmongPop_vs_age_by_gender_fishinland_ftp_",c("", "data_"),"gp_round", round, ".pdf")
        ggsave2(plots_ftp[[1]], file = filenames[[1]], w = 6, h = 5)
        ggsave2(plots_ftp[[2]], file = filenames[[2]], w = 6, h = 5)

        plots_all <- .plot.stan.fit(
            nspop.byage.ptype[PTYPE == "all"],
            DT2=ppDT[PTYPE == "all"],
            ylims = c(0,.4),
            ylab = "population with unsuppressed viral load\n(95% credible interval)\n"
        )
        filenames <- paste0("notsuppAmongPop_vs_age_by_gender_fishinland_all_",c("", "data_"),"gp_round", round, ".pdf")
        ggsave2(plots_all[[1]], file = filenames[[1]], w = 6, h = 5)
        ggsave2(plots_all[[2]], file = filenames[[2]], w = 6, h = 5)

        catn("extract basic not supp estimates")
        # ______________________________________

        draws <- subset(re, select=names(re) %like% '^p_predict|\\.draw')
        draws <- melt( draws,  id.vars='.draw', value.name='P')
        draws[ , `:=` (
            AGE_LABEL = .stan.brackets.to.age(variable, .stan.data=stan.data),
            variable = .stan.remove.brackets(variable)
        ) ]
        draws <- .stan.get.sex.and.loc(draws, 'variable', codes=group_codes)
        draws[, variable := NULL ]

        draws_wholepop <- get.draws.wholepop(DRAWS=draws, DPART=dpartrates[ROUND==round])
        quants_wholepop <- draws_wholepop[, quantile2(P), by=c("AGE_LABEL", "SEX_LABEL", "LOC_LABEL")] 
        p <- plot.quantiles.wholepop(quants_wholepop, ylim=.4,ylab="population with unsuppressed viral load\n(95% credible interval)\n")
        filename <- paste0("notsuppAmongPop_vs_age_by_gender_fishinland_joint_round", round, ".pdf")
        cmd <- ggsave2(p, file = filename, w = 9, h = 6)

        rp <- merge(
            draws, DT[, .(SEX, LOC, AGE_LABEL, N )],
            by=c('SEX', 'LOC', "AGE_LABEL"))
        rp <- rp[, .(P = sum(P * N) / sum(N)), by = c("LOC", "SEX", ".draw")]

        nspop.by.sex.loc <- rp[, quantile2(P, ps = q) , by=c("LOC", "SEX")]
        nspop.by.sex.loc [, LABEL := prettify_cell(M*100, CL*100, CU*100, percent = TRUE) ]
        nspop.by.sex.loc <- merge(nspop.by.sex.loc, group_codes, by=c('LOC', 'SEX'))

        catn("extract risk ratio of suppressed VL female:male and male:female")
        # _____________________________________________________________________

        rp <- merge(group_codes, rp, by = c("LOC", "SEX"))
        rp <- dcast.data.table(rp, LOC_LABEL + .draw ~ SEX_LABEL, value.var = "P")
        rp[,  `:=` (PR_FM = F / M, PR_MF = M / F, M = NULL, F = NULL)]

        rp <- melt(rp, id.vars = c("LOC_LABEL", ".draw"))
        rp <- rp[, quantile2(value, ps = q), by = c("LOC_LABEL", "variable")]
        rp[, LABEL := prettify_cell(M*100, CL*100, CU*100, percent = TRUE) ]

        nspop.by.loc <- copy(rp)

        catn("extract risk ratio of unsuppressed VL F:M and M:F by age")
        # ______________________________________________________________

        rp <- dcast.data.table(draws, LOC_LABEL + .draw + AGE_LABEL ~ SEX_LABEL, value.var = "P")
        rp[,  `:=` (PR_FM = F / M, PR_MF = M / F, M = NULL, F = NULL)]
        rp <- melt(rp, 
            id.vars = c("LOC_LABEL", "AGE_LABEL", ".draw"),
            measure.vars = c("PR_FM", "PR_MF"))
        rp <- rp[, quantile2(value, ps = q), by = c("LOC_LABEL", "AGE_LABEL","variable")]
        rp[, LABEL := prettify_cell(M, CL, CU)]
        nspop.ratio.by.loc.age <- copy(rp)


        catn("extract if diff in F:M riskratio of unsuppr VL is different in fish vs inland")
        # ___________________________________________________________________________________

        rp <- merge(draws, DT[, .(LOC_LABEL, SEX_LABEL, AGE_LABEL, N)], by=c('LOC_LABEL', "SEX_LABEL", "AGE_LABEL"))
        rp <- rp[, list(P = sum(P * N) / sum(N)), by = c("LOC_LABEL", "SEX_LABEL", ".draw")]
        rp <- dcast.data.table(rp, LOC_LABEL + .draw ~ SEX_LABEL, value.var = "P")
        rp[,  `:=` (PR_FM = F / M, PR_MF = M / F, M = NULL, F = NULL)]
        rp <- dcast.data.table(rp, .draw ~ LOC_LABEL, value.var = "PR_FM")
        rp[, PR_FM_D := fishing - inland]
        rp[, quantile2(PR_FM_D, q)]

        filename <- paste0("suppAmongPop_round", round, ".rda")
        save(DT, stan.data, re, nspop.byage.ptype,
            nspop.by.sex.loc, nspop.ratio.by.loc.age,
            file = file.path(vl.out.dir., filename)
        )


        catn("make table version suppressed")
        # ___________________________________

        nspop.byage.ptype[, LABEL := .write.CIs(M, CL, CU, percent = T, d = 1)]
        setnames(nspop.ratio.by.loc.age, "LABEL", "LABEL2")

        vec_ages <- c(20.5, 25.5, 30.5, 35.5, 40.5, 45.5)

        dt <- subset(nspop.byage.ptype, AGE_LABEL %in% vec_ages) |>
            dcast.data.table(LOC_LABEL + AGE_LABEL ~ SEX_LABEL, value.var = "LABEL")
        .f <- function(var){
            tmp <- subset(nspop.ratio.by.loc.age,
                variable == var & AGE_LABEL %in% vec_ages,
                select=c('LOC_LABEL', 'AGE_LABEL', 'LABEL2')
            )
            tmp <- setnames(tmp, "LABEL2", var)
            tmp
        }
        tmp_fm <- .f("PR_FM"); tmp_mf <- .f("PR_MF")
        dt <- merge(dt, tmp_fm, by = c("LOC_LABEL", "AGE_LABEL"))
        dt <- merge(dt, tmp_mf, by = c("LOC_LABEL", "AGE_LABEL"))

        filename <- paste0("tab_suppAmongPop_round", round, ".csv")
        fwrite(dt, row.names = FALSE, file = file.path(vl.out.dir., filename))

        cat("Round", round, ": done.\n")
        return(TRUE)
    }

    if (parallelise) {
        # Use parallel backend if run locally
        foreach(
            r = vla[, unique(ROUND)],
            .combine = "c"
        ) %dopar% {
            cat("Running Round", r, "\n")
            .fit.stan.and.plot.by.round(vla[ROUND == r, ])
        } -> tmp
        return(tmp)
    } else {
        # else do not
        for (r in args$round) {
            cat("Running Round", r, "\n")
            .fit.stan.and.plot.by.round(vla[ROUND == r, ])
        }
    }
}
