.stan.get.sex.and.loc <- function(DT, var='variable', var_overwrite=TRUE, codes=NULL){
    .s <- function(reg, x) gsub("^(.*)_([0-9])([0-9]).*$", reg, x)
    DT[, SEX := as.integer(.s("\\2", get(var)))]
    DT[, LOC := as.integer(.s("\\3", get(var)))]
    if( DT[[var]] %like% '_ftp$' |> any() ){
        DT[, PTYPE := fifelse( get(var) %like% '_ftp$', yes="ftp", no="all") ]
    }
    if(var_overwrite){
        DT[, (var) := .s("\\1", get(var))]
    }
    if(! is.null(codes)){
        DT <- merge(DT, codes, by=c("SEX", "LOC"))
    }
    return(DT)
}


extract.stan.hyperparams.rho <- function(re, encoding){

    
    ps=c(CL=.025, IL=.25,M=.5, IU=.75,CU=.975)
    .q <- function(x) transpose(as.data.table(quantile(x, probs = ps)))

    names_pars <- names(re)
    if( is.null(names_pars) ){
        names_pars <- dimnames(re)$variable
    }

    names_pars <- grep("rho_[0-9]+$|alpha_[0-9]+$", names_pars, value = TRUE)
    pars_quantiles <- lapply(as.list(re)[names_pars], .q)
    pars_quantiles <- rbindlist(pars_quantiles, idcol = "GP_hyper_par")
    names(pars_quantiles) <- c('GP_hyper_par', names(ps))


    # .stan.get.sex.and.loc
    .s <- function(reg, x) gsub("^([a-z]+)_([0-9])([0-9])", reg, x)
    pars_quantiles[, SEX := as.integer(.s("\\2", GP_hyper_par))]
    pars_quantiles[, LOC := as.integer(.s("\\3", GP_hyper_par))]
    pars_quantiles[, GP_hyper_par := .s("\\1", GP_hyper_par)]

    out <- merge(encoding, pars_quantiles, by = c("SEX", "LOC"))
    out [, col := "Posterior"]
    return(out)

}

.get.prior.ranges <- function(stan.data, DT, shape, scale) {
    # Extracts alpha and rho from stan.data
    # Then computes 95% intervals assuming inv.gamma and normal prior distributions
    cols <- grep("alpha_hyper", names(stan.data), value = TRUE)
    tmp <- unlist(as.vector(stan.data[cols]))
    tmp <- data.table(par = names(tmp), val = as.vector(tmp))

    .f <- function(reg, x) {
        as.integer(gsub("^.*?_([0-9])([0-9])$", reg, x))
    }

    tmp[, `:=`(
        SEX = .f("\\2", par),
        LOC = .f("\\1", par),
        GP_hyper_par = gsub("^.*?(alpha|rho).*?$", "\\1", par)
    )]

    tmp <- rbind(
        tmp,
        data.table(
            par = c("rho_hyper_scale", "rho_hyper_shape"),
            val = c(scale, shape), SEX = NA,
            LOC = NA, GP_hyper_par = "rho"
        )
    )
    tmp

    ps <- c(0.025, 0.25, 0.5, 0.75, 0.975)
    cols <- c("CL", "IL", "M", "IU", "CU")

    # tmp[GP_hyper_par=='rho', (cols) := as.list(invgamma::qinvgamma(p=ps, shape=val, scale=val)), by=val]

    qhalfnorm <- function(p, sd) {
        x <- abs(rnorm(1000000, mean = 0, sd = sd))
        as.list(quantile(x, ps))
    }

    qinvgamma <- function(p, sh, sc) {
        x <- MCMCpack::rinvgamma(1000000, shape = sh, scale = sc)
        # CHECK ps2 <- c(0.01, 0.99)
        as.list(quantile(x, ps))
    }

    tmp[GP_hyper_par == "alpha", (cols) := qhalfnorm(p = ps, sd = val), by = val]
    tmp[GP_hyper_par == "rho", (cols) := qinvgamma(p = ps, sh = shape, sc = scale)]

    tmp1 <- tmp[GP_hyper_par == "rho", list(
        GP_hyper_par = "rho",
        SEX = c(0, 0, 1, 1),
        LOC = c(0, 1, 0, 1),
        CL = CL[1],
        CU = CU[1]
    )]

    cols <- names(tmp1)
    tmp <- rbind(tmp1, tmp[GP_hyper_par == "alpha", ..cols])

    tmp1 <- unique(DT[, .(SEX, SEX_LABEL, LOC, LOC_LABEL)])
    tmp <- merge(tmp1, tmp, by = c("SEX", "LOC"))
    tmp[, col := "prior"]
    tmp
}

.make.stan.data.gp <- function(DTsd,
                               num.var = NA,
                               den.var = NA,
                               alpha_hyper_sd = .75,
                               rho_hyper_lower_bound = 1, rho_hyper_upper_bound = 35 ) {
    stopifnot(length(num.var == 1) & length(den.var == 1))

    tmp <- seq(DTsd[, min(AGE_LABEL)], DTsd[, max(AGE_LABEL) + 1], .5)
    tmp1 <- which(tmp %% 1 == 0.5)

    # num.var='HIV_N'
    # den.var='N'
    stan.data <- list(
        x_predict = tmp,
        N_predict = length(tmp),
        observed_idx = tmp1,
        N_observed = length(tmp1),
        y_observed_00 = DTsd[SEX == 0 & LOC == 0, ..num.var][[1]],
        y_observed_10 = DTsd[SEX == 1 & LOC == 0, ..num.var][[1]],
        y_observed_01 = DTsd[SEX == 0 & LOC == 1, ..num.var][[1]],
        y_observed_11 = DTsd[SEX == 1 & LOC == 1, ..num.var][[1]],
        total_observed_00 = DTsd[SEX == 0 & LOC == 0, ..den.var][[1]],
        total_observed_10 = DTsd[SEX == 1 & LOC == 0, ..den.var][[1]],
        total_observed_01 = DTsd[SEX == 0 & LOC == 1, ..den.var][[1]],
        total_observed_11 = DTsd[SEX == 1 & LOC == 1, ..den.var][[1]],
        alpha_hyper_par_00 = alpha_hyper_sd,
        alpha_hyper_par_10 = alpha_hyper_sd,
        alpha_hyper_par_01 = alpha_hyper_sd,
        alpha_hyper_par_11 = alpha_hyper_sd,
        rho_hyper_lower_bound = rho_hyper_lower_bound,
        rho_hyper_upper_bound = rho_hyper_upper_bound
    )

    stan.data
}

.get.summary.without.hyperparams <- function(FIT, verbose=TRUE, vars)
{
    # weird bug whereby names_vars was not recognised, so using vars instead.

    summ <- {
        if( "CmdStanFit" %in% class(FIT) ){
            FIT$summary(vars %which.not.like% "rho_hyper_par|L_cov|^\\.")
        }
    }  |> as.data.table()

    .msg <- summ[, sprintf(
            "The minimum effective sample size is: %s.\nOut of %s parameters:\n- %s had ess_tail < 4000,\n- %s had ess_tail < 1000\n",
            min(ess_tail, na.rm = T), .N, 
            sum(ess_tail < 4000, na.rm = TRUE),
            sum(ess_tail < 1000, na.rm = TRUE)
        ) ]

    if(verbose){
        cat(.msg)
    }

    return(summ)
}

.stan.remove.brackets <- function(x){
    x <- gsub( "\\[.*\\]", "\\1", x)
    return(x)
}

.stan.brackets.to.age <- function(x, .stan.data) {

    x <- gsub( "^.*\\[([0-9]*)\\]$", "\\1", x)

    dict <- unique( .stan.data$x_predict)
    names(dict) <- as.character(seq_along(.stan.data$x_predict))

    out <- unname(dict[x])

    if( any(is.na(out)))
        warning("Some ages were not found in the dictionary")

    return(out)
}

.plot.gp.hyperparameters <- function(GP, PR) {

    p <- ggplot(GP, aes(
        colour = col,
        x = paste0(GP_hyper_par, " ", LOC_LABEL, " ", SEX_LABEL)
    )) +
        geom_errorbar(data = PR, aes(ymin = CL, ymax = CU)) +
        geom_linerange(aes(ymin = CL, ymax = CU)) +
        geom_point(aes(y = M)) +
        coord_flip() +
        theme_bw() +
        theme(legend.position = "bottom") +
        labs(
            x = "GP hyperparameter\n", y = "", colour = "",
            title = "95% Credible Intervals"
        )
    p
}

.plot.stan.fit <- function(DT, DT2=ppDT, ylab, ylims=NA, round=NA) {

    tmp <- prettify_labels(DT)
    tmp2 <- prettify_labels(DT2)

    if( is.na(round) ){
        facet.formula <- formula( . ~ .  )
    }else{
        facet.formula <- formula( ROUND ~ . )
    }

    ALPHA = .3

    p_without <- ggplot(tmp, aes(x = AGE_LABEL, colour = SEX_LAB, fill = SEX_LAB, y=M, ymin=CL, ymax=CU)) +
        geom_ribbon(alpha = ALPHA, col = NA) +
        geom_line() +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(labels = scales:::percent, expand = c(0, 0)) +
        scale_colour_manual(values = palettes$sex, labels=sex_dictionary2 ) + 
        scale_fill_manual(values = palettes$sex, labels=sex_dictionary2 ) +
        coord_cartesian(ylim = ylims, expand = FALSE) +
        facet_wrap(~LOC_LAB, ncol = 2) +
        theme_default() +
        my_labs(
            x = "age at visit (years)",
            y = ylab,
        )

    p_with <- p_without + 
        geom_point(data = DT2, aes(y = M, shape = SEX_LAB, size = N)) +
        scale_size(range = c(0, 3)) +
        scale_shape_manual(values=shapes$sex, labels=sex_dictionary2) + 
        labs(pch = "Gender", size = "Number of individuals")

    return(list(p_without, p_with))
}

.plot.stan.fit.ratio <- function(DT, ylab = NA_character_){

    ALPHA = .3

    ggplot(DT, aes(x = AGE_LABEL, ymin = CL, ymax = CU, group = LOC_LABEL)) +
        geom_hline(yintercept = 1, linetype = 2) +
        geom_ribbon(alpha = ALPHA) +
        geom_line(aes(x = AGE_LABEL, y = M)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_log10(expand = c(0, 0)) +
        coord_cartesian(ylim = c(.5, 50)) +
        facet_wrap(~LOC_LABEL, ncol = 2) +
        theme_default() +
        my_labs(
            x = "age at visit (years)",
            y = ylab
        )
}

.write.CIs <- function(m, l, u, d = 1, percent = F) {
    .r <- function(x) round(x * multiplier, digits = d)
    if (percent) {
        multiplier <- 100
        out <- paste0(.r(m), "% (", .r(l), "% - ", .r(u), "%)")
    } else {
        multiplier <- 1
        paste0(.r(m), " [", .r(l), "-", .r(u), "]")
    }
}

date2numeric <- function(x) {
    if (!class(x) %in% c("Date", "character")) {
        return(x)
    }
    x <- as.POSIXlt(x)
    tmp <- x$year + 1900
    x <- tmp + round(x$yday / ifelse((tmp %% 4 == 0 & tmp %% 100 != 0) | tmp %% 400 == 0, 366, 365), d = 3)
    x
}

.preprocess.ds.oli <- function(DT, rm.na.vl = TRUE) {
    # DT <- copy(dall)
    DT <- subset(DT, AGEYRS <= 50)

    if (rm.na.vl) {
        DT <- DT[, {
            idx_posnovl <- HIV_STATUS == 1 & HIV_AND_VL == 0
            np <- sum(HIV_STATUS)
            n <- sum(idx_posnovl)
            sprintf("removing %d out of %d out of HIV positive individuals without VL.", n, np)
            .SD[!idx_posnovl]
        }, ]
    }

    # consider only ARVMED for infected
    set(DT, DT[, which(ARVMED == 1 & HIV_STATUS == 0)], "ARVMED", 0)

    # define VL_COPIES for uninfected
    set(DT, NULL, "VLC", DT$VL_COPIES)
    set(DT, DT[, which(HIV_STATUS == 0)], "VLC", 0)

    # define undetectable VL (machine-undetectable)
    # define suppressed VL (according to WHO criteria)
    DT[, `:=`(
        VLU = as.integer(VLC < VL_DETECTABLE),
        VLS = as.integer(VLC < VIREMIC_VIRAL_LOAD),
        VLD = as.integer(VLC >= VL_DETECTABLE),
        VLNS = as.integer(VLC >= VIREMIC_VIRAL_LOAD)
    )]
    DT[, HIV_AND_VLD := as.integer(VLD == 1 & HIV_AND_VL == 1)]

    # reset VLC below machine detectable to 0
    DT[which(HIV_AND_VL == 1 & VLU == 1), VLC := 0]

    setkey(DT, ROUND, FC, SEX, AGEYRS)

    DT
}

vl.vlprops.by.comm.gender.loc <- function(DT, write.csv = FALSE) {
    # DT <- copy(dall)
    DT <- .preprocess.ds.oli(DT)
    # merge two communities that fully overlap, so we have 40 communities in the end
    DT[COMM_NUM == 22, COMM_NUM := 1]

    # get.community types
    dcomm <- .get.dcomm()

    # calculate HIV seroprevalence and proportion not suppressed of HIV+ by community and gender
    .f <- function(x, y) {
        as.vector(unname(binconf(sum(x), length(y))))
    }

    vlc <- DT[,
        {
            z <- .f(HIV_STATUS == 1, HIV_STATUS)
            z2 <- .f(VLNS == 1, VLNS)
            z3 <- .f(VLNS == 1, which(HIV_STATUS == 1))
            list(
                FC = FC[1],
                N = length(HIV_STATUS),
                PHIV_MEAN = z[1],
                PHIV_CL = z[2],
                PHIV_CU = z[3],
                PVLNS_MEAN = z2[1],
                PVLNS_CL = z2[2],
                PVLNS_CU = z2[3],
                PVLNSofHIV_MEAN = z3[1],
                PVLNSofHIV_CL = z3[2],
                PVLNSofHIV_CU = z3[3],
                VLC_MEAN = mean(VLC)
            )
        },
        by = c("ROUND", "COMM_NUM", "SEX")
    ]

    .f <- function(m, l, u) {
        .r <- function(x) round(x * 100, digits = 1)
        paste0(.r(m), " [", .r(l), "-", .r(u), "]")
    }

    vlc[, PHIV_L := .f(PHIV_MEAN, PHIV_CL, PHIV_CU)]
    vlc[, PVLNS_L := .f(PVLNS_MEAN, PVLNS_CL, PVLNS_CU)]
    vlc[, PVLNSofHIV_L := .f(PVLNSofHIV_MEAN, PVLNSofHIV_CL, PVLNSofHIV_CU)]

    setkey(vlc, SEX, PHIV_MEAN)
    vlc[, SEX := factor(SEX, levels = c("M", "F"), labels = c("men", "women"))]
    names(dcomm)
    vlc <- merge(vlc, dcomm[, .(COMM_NUM, FC2 = TYPE)], by = "COMM_NUM", all.x = TRUE)


    p_inf <- ggplot(vlc, aes(x = PHIV_MEAN, y = PVLNSofHIV_MEAN)) +
        scale_x_continuous(labels = scales:::percent) +
        scale_y_continuous(labels = scales:::percent) +
        geom_errorbar(aes(ymin = PVLNSofHIV_CL, ymax = PVLNSofHIV_CU), alpha = 0.2) +
        geom_errorbarh(aes(xmin = PHIV_CL, xmax = PHIV_CU), alpha = 0.2) +
        geom_point(aes(colour = FC2)) +
        geom_text(aes(label = COMM_NUM), size = 2) +
        facet_grid(ROUND ~ SEX) +
        # scale_colour_manual(values=palettes$comm) +
        scale_colour_manual(values = palettes$comm2) +
        theme_default() +
        my_labs(color = "Community type")

    filename <- "220729_hivnotsuppofhiv_vs_hivprev_by_round_gender_fishinland.pdf"
    ggsave2(p_inf, file = filename, LALA = glm.out.dir, w = 9, h = 12)


    p_pop <- ggplot(vlc, aes(x = PHIV_MEAN, y = PVLNS_MEAN)) +
        scale_x_continuous(labels = scales:::percent) +
        scale_y_continuous(labels = scales:::percent) +
        geom_errorbar(aes(ymin = PVLNS_CL, ymax = PVLNS_CU), alpha = 0.2) +
        geom_errorbarh(aes(xmin = PHIV_CL, xmax = PHIV_CU), alpha = 0.2) +
        geom_point(aes(y = PVLNS_MEAN, colour = FC2)) +
        geom_text(aes(y = PVLNS_MEAN, label = COMM_NUM), size = 2) +
        facet_grid(ROUND ~ SEX) +
        scale_colour_manual(values = palettes$comm2) +
        theme_default() +
        my_labs(colour = "Community type")

    filename <- "220729_hivnotsuppofpop_vs_hivprev_by_round_gender_fishinland.pdf"
    ggsave2(p_pop, file = filename, LALA = glm.out.dir, w = 9, h = 12)

    if (write.csv) {
        # write results to file
        filename <- file.path(
            glm.out.dir,
            "220729_hivnotsuppofhiv_vs_hivprev_by_round_gender_fishinland.csv"
        )
        fwrite(vlc, file = filename)
    }

    list(DT = vlc, p_among_pop = p_pop, p_among_inf = p_inf)
}

vl.prevalence.by.gender.loc.age.gp.cmdstan <- function(
    DT=vla, 
    refit = FALSE,
    vl.out.dir. = vl.out.dir,
    alpha_hyper = .75
) {

    cat("\n\n--- Analysing HIV+ Prevalence ---\n\n")

    vla <- copy(DT)

    # Stan file locations
    file.stan <- file.path(gitdir.stan, "vl_binomial_gp2.stan")

    .fit.stan.and.plot.by.round <- function(DT) {

        #  DT <- copy(vla[ROUND == 16]); vl.out.dir.=vl.out.dir
        round <- DT[, unique(ROUND)]
        stopifnot(length(round) == 1)

        stan.data <- .make.stan.data.gp(
            DTsd = DT,
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

        plots <- .plot.stan.fit(
            prev.hiv.by.age, 
            DT2=ppDT,
            ylims = c(0,.75),
            ylab="HIV seroprevalence (95% credible intervals)")

        filenames <- paste0("hivprevalence_vs_age_by_gender_fishinland_",c("", "data_"),"gp_round", round, ".pdf")

        ggsave2(plots[[1]], file = filenames[[1]], w = 6, h = 5)
        ggsave2(plots[[2]], file = filenames[[2]], w = 6, h = 5)

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

        rp <- merge(
            draws, DT[, .(SEX, LOC, AGE_LABEL, N )],
            by=c('SEX', 'LOC', "AGE_LABEL"))
        rp <- rp[, .(P = sum(P * N) / sum(N)), by = c("LOC", "SEX", ".draw")]

        prev.hiv.by.sex.loc <- rp[, quantile2(P, ps = q) , by=c("LOC", "SEX")]
        prev.hiv.by.sex.loc [, LABEL := prettify_cell(M*100, CL*100, CU*100, percent = TRUE) ]
        prev.hiv.by.sex.loc <- merge(prev.hiv.by.sex.loc, group_codes, by=c('LOC', 'SEX'))


        catn("extract prevalence ratio female:male and male:female")
        # __________________________________________________________

        rp <- merge(group_codes, rp, by = c("LOC", "SEX"))
        rp <- dcast.data.table(rp, LOC_LABEL + .draw ~ SEX_LABEL, value.var = "P")
        rp[,  `:=` (PR_FM = F / M, PR_MF = M / F, M = NULL, F = NULL)]

        rp <- melt(rp, id.vars = c("LOC_LABEL", ".draw"))
        rp <- rp[, quantile2(value, ps = q), by = c("LOC_LABEL", "variable")]
        rp[, LABEL := prettify_cell(M*100, CL*100, CU*100, percent = TRUE) ]

        prevratio.hiv.by.loc <- copy(rp)

        catn("plot prevalence ratio F:M and M:F by age")
        # ______________________________________________


        rp <- dcast.data.table(draws, 
            LOC + LOC_LABEL + .draw + AGE_LABEL ~ SEX_LABEL, value.var = "P")
        rp[, `:=` (PR_FM = F / M, PR_MF = M / F, M = NULL, F = NULL)]
        rp <- melt( rp, 
            id.vars = c("LOC_LABEL", "AGE_LABEL", ".draw"),
            measure.vars = c("PR_FM", "PR_MF"))
        prevratio.hiv.by.loc.age <- rp[, 
            quantile2(value),
            by = c("LOC_LABEL", "AGE_LABEL", "variable")]

        p <- .plot.stan.fit.ratio(
            prevratio.hiv.by.loc.age[variable == "PR_FM"],
            ylab = "female to male HIV seroprevalence ratio\n(95% credible interval)\n"
        )

        filename <- paste0("fit_hivprevalenceratio_vs_age_by_fishinland_stan_round", round, ".pdf")
        ggsave2(p, file = filename, w = 8, h = 5)

        catn("extract if diff in F:M prevalence risk ratio in fish vs inland")
        # ____________________________________________________________________

        by_cols <- c("LOC", "SEX", "AGE_LABEL")
        rp <- merge(subset(DT, select = c(by_cols, 'N') ), draws, by = by_cols )
        rp <- rp[, 
            .(P = sum(P * N) / sum(N)),
            by = c("LOC_LABEL", "SEX_LABEL", ".draw")]
        rp <- dcast.data.table(rp, 
            LOC_LABEL + .draw ~ SEX_LABEL, value.var = "P")

        rp[, `:=` (PR_FM = F / M, PR_MF = M / F, M = NULL, F = NULL)]

        rp <- dcast.data.table(rp, .draw ~ LOC_LABEL, value.var = "PR_FM")
        rp[, {
            PR_FM_D <- inland - fishing
            quantile2(PR_FM_D, ps = q)
        }]

        # Save diff. qttities
        filename <- paste0("hivprevalence_round", round, ".rda")
        save(DT, stan.data, re, prev.hiv.by.age, prevratio.hiv.by.loc,
            prev.hiv.by.sex.loc, prevratio.hiv.by.loc.age,
            file = file.path(vl.out.dir., filename)
        )

        catn("make table version suppressed")

        # prev.hiv.by.age[, LABEL := .p1(M, CL, CU)]
        prev.hiv.by.age[, LABEL := prettify_cell(M*100, CL*100, CU*100, percent=TRUE) ]
        prevratio.hiv.by.loc.age[, LABEL2 := prettify_cell(M, CL, CU)]

        vec_ages <- c(20.5, 25.5, 30.5, 35.5, 40.5, 45.5)

        dt <- subset(prev.hiv.by.age, AGE_LABEL %in% vec_ages) |> 
            dcast.data.table(
                LOC_LABEL + AGE_LABEL ~ SEX_LABEL,
                value.var = "LABEL" )

        .f <- function(var){
            tmp <- subset( prevratio.hiv.by.loc.age,
                variable == var & AGE_LABEL %in% vec_ages,
                select=c('LOC_LABEL', 'AGE_LABEL', 'LABEL2')
            )
            tmp <- setnames(tmp, "LABEL2", var)
            tmp
        }
        tmp_fm <- .f("PR_FM"); tmp_mf <- .f("PR_MF")

        dt <- merge(dt, tmp_fm, by = c("LOC_LABEL", "AGE_LABEL"))
        dt <- merge(dt, tmp_mf, by = c("LOC_LABEL", "AGE_LABEL"))

        filename <- file.path(vl.out.dir., paste0("tab_hivprevalence_round", round, ".csv"))
        fwrite(dt, row.names = FALSE, file = filename)

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

vl.suppofinfected.by.gender.loc.age.gp.cmdstan <- function(
    DT=vla,
    refit = FALSE,
    vl.out.dir. = vl.out.dir,
    alpha_hyper = .75
) {

    cat("\n\n--- Analyse suppressed among infected ---\n\n")

    vla <- copy(DT)

    file.stan <- file.path(gitdir.stan, "vl_binomial_gp2.stan")

    .fit.stan.and.plot.by.round <- function(DT) {

        # DT <- copy(vla[ROUND == 16] )
        round <- DT[, unique(ROUND)]
        stopifnot(length(round) == 1)

        # Exclude people with unknown Viral load from analyses
        DT[, HIV_N := HIV_N - VLNA_N ]
        stopifnot(DT[, all(HIV_N>=0)])
        DT[, VLSUP_N := HIV_N - VLNS_N, ]
        stan.data <- .make.stan.data.gp(DT, 
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

        catn("compare to self-report")
        # ____________________________

        arv_bool <- DT[, any(ARV_N) > 0]
        if (arv_bool) {

            DT[, HIV_NARV_N := HIV_N - ARV_N]
            stan.data.2 <- .make.stan.data.gp(DT,
                num.var = "HIV_NARV_N",
                den.var = "HIV_N",
                alpha_hyper_sd = alpha_hyper
            )

            
            # file paths for fit and standata
            .prefix <- file.path(vl.out.dir., "cmd_notARVAmongInfected_gp_stan_round")
            filename.2 <- paste0(.prefix, round, ".rds")
            filename_standata.2 <- paste0(.prefix, round, "_standata.rds")

            if (file.exists(filename.2) & refit == FALSE) {

                cat("Loading previously run HMC... \n")
                fit2 <- readRDS(filename.2)
                stan.data.2 <- readRDS(filename_standata.2)

            } else {

                stan.model <- cmdstan_model(stan_file = file.stan, cpp_options=list(stan_threads = TRUE))
                stan.args <- yaml::read_yaml(path.stan.config)

                fit2 <- stan.model$sample(
                    data = stan.data.2,
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
                fit2$save_object(file = filename.2)
                saveRDS(object=stan.data.2, file = filename_standata.2)
            }
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
        nsinf.by.age <- .stan.get.sex.and.loc(nsinf.by.age, 'variable', codes=group_codes)

        plots <- .plot.stan.fit(
            nsinf.by.age,
            DT2=ppDT,
            ylims = c(0,1),
            ylab = "HIV+ individuals with suppressed viral load\n(95% credible interval)\n"
        )

        filenames <- paste0("suppAmongInfected_vs_age_by_gender_fishinland_",c("", "data_"),"gp_round", round, ".pdf")
        ggsave2(plots[[1]], file = filenames[[1]], w = 6, h = 5)
        ggsave2(plots[[2]], file = filenames[[2]], w = 6, h = 5)

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

        rp <- merge(
            draws, DT[, .(SEX, LOC, AGE_LABEL, N )],
            by=c('SEX', 'LOC', "AGE_LABEL"))
        rp <- rp[, .(P = sum(P * N) / sum(N)), by = c("LOC", "SEX", ".draw")]

        nsinf.by.sex.loc <- rp[, quantile2(P, ps = q) , by=c("LOC", "SEX")]
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
        if (!arv_bool) {
            nainf.by.age <- re2 <- data.table()
        }
        save(DT, stan.data, re, re2, nainf.by.age,
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

vl.suppofpop.by.gender.loc.age.gp.cmdstan <- function(
    DT, 
    refit = FALSE,
    vl.out.dir. = vl.out.dir,
    alpha_hyper = .75
){

    cat("\n\n--- Analyse suppression among participants ---\n\n")

    # DT <- copy(dall); refit=FALSE
    DT <- .preprocess.ds.oli(DT)
    tmp <- c("N", "HIV_N", "VLNS_N", "ARV_N")
    vla <- .preprocess.make.vla(DT, select = tmp)

    # Stan file locations
    file.stan <- file.path(gitdir.stan, "vl_binomial_gp2.stan")

    .fit.stan.and.plot.by.round <- function(DT) {

        #  DT <- copy(vla[ROUND == 16]); refit=FALSE
        round <- DT[, unique(ROUND)]
        stopifnot(length(round) == 1)

        # Remove people with unknwon Viral load,
        # else equivalent to assuming they are suppressed or negative
        DT[, N := N - VLNA_N]

        stan.data <- .make.stan.data.gp(DT,
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

        catn("make prevalence plot by age")
        # _________________________________

        q <- c("M"=.5, "CL"=.025, "CU"=.975)
        cols <- names(re) %which.like% '^p_predict_'
        tmp <- re[, lapply(.SD, posterior::quantile2, probs=q), .SDcols =cols]
        tmp[, quantile := names(q)]
        tmp <- melt( tmp, id.vars = "quantile")
        nspop.by.age <- dcast.data.table(tmp, variable ~ quantile, value.var = "value")
        nspop.by.age[ , `:=` (
            AGE_LABEL = .stan.brackets.to.age(variable, .stan.data=stan.data),
            variable = .stan.remove.brackets(variable)
        ) ]
        nspop.by.age <- .stan.get.sex.and.loc(nspop.by.age, 'variable', codes=group_codes)

        plots <- .plot.stan.fit(
            nspop.by.age, 
            DT2=ppDT,
            ylims = c(0,.4),
            ylab = "population with unsuppressed viral load\n(95% credible interval)\n")

        filenames <- paste0("hivprevalence_vs_age_by_gender_fishinland_",c("", "data_"),"gp_round", round, ".pdf")

        ggsave2(plots[[1]], file = filenames[[1]], w = 6, h = 5)
        ggsave2(plots[[2]], file = filenames[[2]], w = 6, h = 5)

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
        save(DT, stan.data, re, nspop.by.age,
            nspop.by.sex.loc, nspop.ratio.by.loc.age,
            file = file.path(vl.out.dir., filename)
        )


        catn("make table version suppressed")
        # ___________________________________

        nspop.by.age[, LABEL := .write.CIs(M, CL, CU, percent = T, d = 1)]
        setnames(nspop.ratio.by.loc.age, "LABEL", "LABEL2")

        vec_ages <- c(20.5, 25.5, 30.5, 35.5, 40.5, 45.5)

        dt <- subset(nspop.by.age, AGE_LABEL %in% vec_ages) |>
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

vl.vldistribution.by.gender.loc <- function() {

    # DT <- copy(dall)
    DT <- .preprocess.ds.oli(DT)

    tmp <- seq.int(min(DT$AGEYRS), max(DT$AGEYRS))
    tmp1 <- DT[, sort(unique(ROUND))]

    # plot proportion of population with viral load > x
    x <- seq(log(1), log(max(DT$VLC)), length.out = 1e3)
    x <- c(0, exp(x))

    vld <- as.data.table(expand.grid(
        ROUND = tmp1,
        X = x,
        SEX = c("M", "F"),
        FC = c("fishing", "inland"),
        HIV_AND_VLD = c(0, 1)
    ))

    ans <- vld[,
        {
            n <- sum(DT$ROUND == ROUND, DT$SEX == SEX & DT$FC == FC & DT$HIV_AND_VLD >= HIV_AND_VLD)
            k <- sum(DT$ROUND == ROUND, DT$SEX == SEX & DT$FC == FC & DT$HIV_AND_VLD >= HIV_AND_VLD & X < DT$VLC)
            z <- as.vector(unname(binconf(k, n)))
            list(N = n, K = k, P_M = z[1], P_CL = z[2], P_CU = z[3])
        },
        by = names(vld)
    ]

    set(ans, NULL, "HIV_AND_VLD", factor(ans[, HIV_AND_VLD], levels = c(0, 1), labels = c("all study participants", "infected study participants\nwith detectable viral load")))
    set(ans, NULL, "SEX", ans[, factor(SEX, levels = c("M", "F"), labels = c("men", "women"))])

    ans <- subset(ans, !(HIV_AND_VLD == "infected study participants\nwith detectable viral load" & X < VL_DETECTABLE))
    ans <- subset(ans, !(HIV_AND_VLD == "all study participants" & X < VL_DETECTABLE))

    p <- ggplot(ans) +
        geom_line(aes(x = X, y = P_M, group = interaction(FC, SEX), colour = SEX, linetype = FC)) +
        scale_x_log10() +
        scale_y_continuous(labels = scales:::percent, expand = c(0, 0)) +
        scale_colour_manual(values = c("men" = "royalblue3", "women" = "deeppink2")) +
        geom_text(aes(x = 1e3, y = P_M * 1.03, label = "")) +
        theme_bw() +
        facet_wrap(~HIV_AND_VLD, scales = "free", ncol = 2) +
        labs(
            x = "\nviral load\n(copies / ml)",
            y = "proportion of individuals with larger viral load\n",
            colour = "gender", linetype = "location"
        )

    # TODO: add `_round?`
    filename <- "220729_vldistribution_by_gender_fishinland.pdf"
    ggsave2(p, file = filename, w = 9, h = 5)
}

vl.keystats.by.gender.loc <- function(DT) {
    # TODO: by_round as well?

    # DT <- copy(dall)
    vl.out.dir <- file.path(vl.out.dir)
    DT <- .preprocess.ds.oli(DT)

    # entire population:
    # mean viral load, proportion with DVL
    DT[, mean(VLC)]
    # 2290.494
    binconf(length(which(DT$VLD == 1)), nrow(DT))
    # PointEst   Lower      Upper
    # 0.05717964 0.05393704 0.06060469

    # stratified by men/women inland/fishing
    # ______________________________________

    .f <- function(x, y) {
        as.vector(unname(binconf(sum(x), length(y))))
    }

    ans <- DT[,
        {
            z <- .f(HIV_STATUS == 1, HIV_STATUS)
            z2 <- .f(VLNS == 1, VLNS)
            z3 <- .f(VLNS == 1, which(HIV_STATUS == 1))

            list(
                N = length(HIV_STATUS),
                PHIV_MEAN = z[1],
                PHIV_CL = z[2],
                PHIV_CU = z[3],
                PVLNS_MEAN = z2[1],
                PVLNS_CL = z2[2],
                PVLNS_CU = z2[3],
                PVLNSofHIV_MEAN = z3[1],
                PVLNSofHIV_CL = z3[2],
                PVLNSofHIV_CU = z3[3],
                VLC_MEAN = mean(VLC)
            )
        },
        by = "SEX"
    ]

    ans[, FC := "overall"]
    tmp <- DT[,
        {
            z <- .f(HIV_STATUS == 1, HIV_STATUS)
            z2 <- .f(VLNS == 1, VLNS)
            z3 <- .f(VLNS == 1, which(HIV_STATUS == 1))

            list(
                N = length(HIV_STATUS),
                PHIV_MEAN = z[1],
                PHIV_CL = z[2],
                PHIV_CU = z[3],
                PVLNS_MEAN = z2[1],
                PVLNS_CL = z2[2],
                PVLNS_CU = z2[3],
                PVLNSofHIV_MEAN = z3[1],
                PVLNSofHIV_CL = z3[2],
                PVLNSofHIV_CU = z3[3],
                VLC_MEAN = mean(VLC)
            )
        },
        by = c("FC", "SEX")
    ]

    ans <- rbind(tmp, ans)

    set(ans, NULL, "FC", factor(ans$FC, levels = c("overall", "fishing", "inland")))
    set(ans, NULL, "SEX", factor(ans$SEX, levels = c("F", "M")))
    setkey(ans, FC, SEX)

    .f <- function(m, l, u) {
        .r <- function(x) round(x * 100, digits = 1)
        paste0(.r(m), " [", .r(l), "-", .r(u), "]")
    }

    ans[, PHIV_L := .f(PHIV_MEAN, PHIV_CL, PHIV_CU)]
    ans[, PVLNS_L := .f(PVLNS_MEAN, PVLNS_CL, PVLNS_CU)]
    ans[, PVLNSofHIV_L := .f(PVLNSofHIV_MEAN, PVLNSofHIV_CL, PVLNSofHIV_CU)]

    # TODO: change name?
    fwrite(ans, file = file.path(vl.out.dir, "220729_keystats_by_gender_fishinland.csv"))
}

vl.vlratio.by.loc <- function(DT) {
    # DT <- copy(dall)
    vl.out.dir <- file.path(vl.out.dir)
    DT <- .preprocess.ds.oli(DT)

    ans <- as.data.table(expand.grid(BS = 1:1e3, FC = c("fishing", "inland")))
    set.seed(42)

    ans <- ans[,
        {
            zm <- which(DT$FC == FC & DT$SEX == "M")
            zf <- which(DT$FC == FC & DT$SEX == "F")
            zm <- sample(zm, length(zm), replace = TRUE)
            zf <- sample(zf, length(zf), replace = TRUE)
            list(VLCM_M = mean(DT$VLC[zm]), VLCM_F = mean(DT$VLC[zf]))
        },
        by = c("FC", "BS")
    ]

    ans[, VLCR := VLCM_M / VLCM_F]

    ans <- ans[, list(
        V = quantile(VLCR, prob = c(0.5, 0.025, 0.975)),
        P = c("M", "CL", "CU")
    ), by = c("FC")]
    ans <- dcast.data.table(ans, FC ~ P, value.var = "V")

    #
    # stratified by men/women inland/fishing
    #
    .f <- function(x, y) {
        as.vector(unname(binconf(sum(x), length(y))))
    }

    ans <- DT[,
        {
            z <- .f(HIV_STATUS == 1, HIV_STATUS)
            z2 <- .f(VLNS == 1, VLNS)
            z3 <- .f(VLNS == 1, which(HIV_STATUS == 1))

            list(
                N = length(HIV_STATUS),
                PHIV_MEAN = z[1],
                PHIV_CL = z[2],
                PHIV_CU = z[3],
                PVLNS_MEAN = z2[1],
                PVLNS_CL = z2[2],
                PVLNS_CU = z2[3],
                PVLNSofHIV_MEAN = z3[1],
                PVLNSofHIV_CL = z3[2],
                VLC_MEAN = mean(VLC)
            )
        },
        by = "SEX"
    ]

    ans[, FC := "overall"]

    tmp <- DT[,
        {
            z <- .f(HIV_STATUS == 1, HIV_STATUS)
            z2 <- .f(VLNS == 1, VLNS)
            z3 <- .f(VLNS == 1, which(HIV_STATUS == 1))

            list(
                N = length(HIV_STATUS),
                PHIV_MEAN = z[1],
                PHIV_CL = z[2],
                PHIV_CU = z[3],
                PVLNS_MEAN = z2[1],
                PVLNS_CL = z2[2],
                PVLNS_CU = z2[3],
                PVLNSofHIV_MEAN = z3[1],
                PVLNSofHIV_CL = z3[2],
                PVLNSofHIV_CU = z3[3],
                VLC_MEAN = mean(VLC)
            )
        },
        by = c("FC", "SEX")
    ]

    ans <- rbind(tmp, ans)
    set(ans, NULL, "FC", factor(ans$FC, levels = c("overall", "fishing", "inland")))
    set(ans, NULL, "SEX", factor(ans$SEX, levels = c("F", "M")))
    setkey(ans, FC, SEX)
}

vl.age.gender <- function() {
    require(data.table)
    prjdir <- "~/Box Sync/OR_Work/2018/2018_RakaiViralLoad"
    infile <- file.path(prjdir, "data", "191101_data_round17_vl_gps.rda")
    load(infile)

    # drop few infecteDT with missing VL
    DT <- subset(DT, HIV_STATUS == 0 | (HIV_STATUS == 1 & !is.na(VL_COPIES)))
    # set VL for uninfected to 0, and VL with undetectable VL to 0
    set(DT, DT[, which(HIV_STATUS == 0)], "VL_COPIES", 0)
    set(DT, DT[, which(HIV_STATUS == 1 & VL_UNDETECTABLE == 1)], "VL_COPIES", 0)

    #
    # calculate proportion with VL > x among participants

    # do general by as characters
    # then determine sort index
    # then calculate empirical quantile
    DT <- DT[order(SEX, VL_COPIES), ]
    DT[VL]

    DT[, sort(unique(VL_COPIES))]
    # dv <- data.table(VL:= )
}

vl.get.data.round17 <- function() {
    require(data.table)
    prjdir <- "~/Box/OR_Work/2018/2018_RakaiViralLoad"
    infile <- "data_raw/ViralLoad_Data_Pangea_Ratmann.rda"
    load(file.path(prjdir, infile))

    # subset to survey round 17
    # _________________________

    ds <- subset(as.data.table(survey_data), visit == 17)
    # reset dates from Date format to numeric
    for (x in c("visit_date", "lastNegDate", "firstPosDate"))
    {
        set(ds, NULL, x, date2numeric(ds[[x]]))
    }
    # make all column names upper case
    setnames(ds, colnames(ds), toupper(colnames(ds)))
    # define FISHING_COMM
    ds[, FC := as.character(factor(COMM_NUM %in% c(770, 771, 774, 38), levels = c(TRUE, FALSE), labels = c("fishing", "inland")))]
    # define ARVMED
    set(ds, ds[, which(ARVMED == 8)], "ARVMED", NA_integer_)
    set(ds, NULL, "ARVMED", ds[, as.integer(as.character(factor(ARVMED, levels = c(1, 2), labels = c("1", "0"))))])

    #
    # prepare GPS coordinates
    #
    dg <- as.data.table(gpsdat)
    # bring dates into common format
    setnames(dg, colnames(dg), gsub("\\.", "_", toupper(colnames(dg))))
    tmp <- which(dg[, grepl("([0-9]+)/([0-9]+)/([0-9]+)", GPS_DATE)])
    set(dg, tmp, "GPS_DATE", dg[tmp, gsub("([0-9]+)/([0-9]+)/([0-9]+)", "\\3-\\1-\\2", GPS_DATE)])
    tmp <- which(dg[, grepl("([0-9]+)-([A-Za-z]+)-([0-9]+)", GPS_DATE)])
    set(dg, tmp, "GPS_DATE", dg[tmp, gsub("([0-9]+)-([A-Za-z]+)-([0-9]+)", "20\\3-\\2-\\1", GPS_DATE)])
    set(dg, NULL, "GPS_DATE", dg[, gsub("Nov", "11", gsub("Oct", "10", gsub("Sep", "09", gsub("Aug", "08", gsub("July", "07", gsub("Jun", "06", gsub("May", "05", GPS_DATE)))))))])
    # reset dates from character format to numeric
    set(dg, NULL, "GPS_DATE", date2numeric(dg[, GPS_DATE]))
    # make households per date unique
    dg <- unique(dg, by = c("HHID", "GPS_DATE"))

    #
    # add to surveyed individuals the GPS of their households
    #
    tmp <- unique(subset(ds, select = c(RCCS_STUDYID, VISIT_DATE, HHID)))
    tmp <- merge(tmp, dg, by = "HHID", all.x = TRUE)
    # some households do not have GPS coordinates
    ch <- subset(tmp, is.na(LATITUDE_JITTER) | is.na(LATITUDE_JITTER))
    if (nrow(ch)) {
        cat("\nNumber of households without GPS coordinates, n=", nrow(ch))
        write.csv(ch, file = file.path(prjdir, "data/check_missing_coordinates.csv"))
        # 	521 households without GPS coordinates
    }
    # for every individual, extract house closest in time
    tmp <- subset(tmp, !is.na(LATITUDE_JITTER) & !is.na(LATITUDE_JITTER))
    tmp2 <- tmp[, list(GPS_DATE = GPS_DATE[which.min(abs(GPS_DATE - VISIT_DATE))[1]]), by = c("RCCS_STUDYID", "VISIT_DATE")]
    tmp <- merge(tmp, tmp2, by = c("RCCS_STUDYID", "VISIT_DATE", "GPS_DATE"))
    stopifnot(nrow(tmp) == nrow(tmp2))
    set(tmp, NULL, c("COMM", "HOUSE"), NULL)
    ds <- merge(ds, tmp, by = c("RCCS_STUDYID", "VISIT_DATE", "HHID"), all.x = TRUE)

    #
    # extract viral loads from round 17
    #
    dvl <- subset(as.data.table(viralLoads), visit == 17)
    setnames(dvl, colnames(dvl), toupper(colnames(dvl)))
    setnames(dvl, c("DATE", "COPIES", "DONEBY"), c("VL_DATE", "VL_COPIES", "VL_DONEBY"))
    set(dvl, NULL, "VL_DATE", date2numeric(dvl[, VL_DATE]))
    stopifnot(!nrow(subset(dvl, is.na(VL_COPIES))))
    # check if one viral load measurement per person
    tmp <- dvl[, list(N_VL = length(VL_COPIES)), by = "RCCS_STUDYID"]
    stopifnot(!nrow(subset(tmp, N_VL > 1)))
    # merge with main data
    set(dvl, NULL, "VISIT", dvl[, as.numeric(VISIT)])
    set(dvl, NULL, "VL_DONEBY", NULL)
    ds <- merge(ds, dvl, by = c("RCCS_STUDYID", "VISIT"), all.x = TRUE)

    # check if viral load for all infected
    ch <- subset(ds, HIV_STATUS == 1 & is.na(VL_COPIES))
    if (nrow(ch)) {
        cat("\nFound infected individuals without VL measurement, n=", nrow(ch))
        write.csv(ch, file = file.path(prjdir, "data/check_missing_viralloads.csv"))
        # 13 HIV+ individuals without VL
    }
    ds[, HIV_AND_VL := as.integer(HIV_STATUS == 1 & !is.na(VL_COPIES))]


    save(ds, file = file.path(prjdir, "data", "191101_data_round17_vl_gps.rda"))
}

prop.dectectable.viraemia <- function() {
    require(data.table)
    require(rgdal)
    require(rgeos)
    library(raster)
    require(RColorBrewer) # Map colours

    # load data
    infile <- "~/Box Sync/OR_Work/2018/2018_RakaiViralLoad/data/merged_round17_vl_gps.rda"
    load(infile)

    tmp <- ds[, list(
        HIV_POS = sum(HIV_STATUS == 1),
        HIV_NEG = sum(HIV_STATUS == 0),
        HIV_PREV = sum(HIV_STATUS == 1) / length(HIV_STATUS)
    ), by = "COMM_NUM"]

    thr <- 1e3
    tmp2 <- ds[HIV_STATUS == 1]
    tmp2 <- tmp2[, list(
        VL_D = sum(VL_COPIES > thr),
        VL_U = sum(VL_COPIES <= thr),
        VL_DP = sum(VL_COPIES > thr) / length(VL_COPIES)
    ), by = "COMM_NUM"]
    tmp <- merge(tmp, tmp2, by = "COMM_NUM")
    tmp[, POP_VL_DP := HIV_PREV * VL_DP]
    ggplot(tmp, aes(y = COMM_NUM, x = POP_VL_DP)) +
        geom_point()

    tmp3 <- subset(ds, HIV_STATUS == 1 & COMM_NUM == 38)
    ggplot(tmp3, aes(x = VL_COPIES)) +
        geom_histogram() +
        facet_grid(~ARVMED)

    tmp3 <- subset(ds, HIV_STATUS == 1 & COMM_NUM == 38 & ARVMED == 2)
    tmp3[, VL_COPIES_C := cut(VL_COPIES, breaks = c(0, 1, 10, 100, 1000, 1e4, 1e5, 1e6, 1e7, 1e10), right = FALSE)]
    tmp3[, table(VL_COPIES_C)]

    tmp3 <- subset(ds, HIV_STATUS == 1 & COMM_NUM == 38 & ARVMED == 1)
    tmp3[, VL_COPIES_C := cut(VL_COPIES, breaks = c(0, 1, 10, 100, 1000, 1e4, 1e5, 1e6, 1e7, 1e10), right = FALSE)]
    tmp3[, table(VL_COPIES_C)]

    tmp3 <- subset(ds, HIV_STATUS == 1 & COMM_NUM == 38)
    tmp3[, table(VL_COPIES > 1)]
}

make.map.190129 <- function(DT) {
    require(data.table)
    require(rgdal)
    require(rgeos)
    library(raster)
    require(RColorBrewer) # Map colours

    # DT <- copy(dall)
    DT <- .preprocess.ds.oli(DT)
    infile <- file.path(indir.deepsequence.data, "RCCS_R15_R18", "Rakai_community_geography_R15.rda")

    # Get longitude and latitude
    tmp <- new.env()
    load(infile, envir = tmp)
    ds <- as.data.table(tmp$comgps)
    cols <- c("latitude", "longitude")
    ds[, (cols) := lapply(.SD, unlist), .SDcols = cols]
    ds[, COMM_NUM := as.integer(COMM_NUM)]
    setnames(ds, cols, paste0(toupper(cols), "_JITTER"))

    # merge with data
    DT <- merge(DT, ds, all.x = T, by = "COMM_NUM")
    DT[is.na(LATITUDE_JITTER), uniqueN(COMM_NUM)] -> tmp
    DT[, cat("missing geoloc for", tmp, "communities\n")]

    # convert the data into a data table
    setnames(DT, "VLU", "VL_UNDETECTABLE")
    DT <- DT[, .(STUDY_ID, SEX, AGEYRS, HIV_STATUS, LATITUDE_JITTER, LONGITUDE_JITTER, VL_COPIES, VL_UNDETECTABLE)]
    # set the NA VL to 0
    DT[is.na(VL_COPIES), VL_COPIES := 0]
    DT[, VL_DETECTABLE := as.numeric(VL_COPIES >= 1000)]
    DT[, RCCS_STUDYID2 := seq_len(.N)]

    ##############################
    # Load in Uganda Shape files #
    ##############################

    uganda1 <- raster::getData("GADM", country = "UGA", level = 1) # Admin unit 1
    uganda3 <- raster::getData("GADM", country = "UGA", level = 3)

    rakai1 <- subset(uganda1, NAME_1 == "Rakai")
    rakai3 <- subset(uganda3, NAME_1 == "Rakai")
    masaka1 <- subset(uganda1, NAME_1 == "Masaka")
    # Create a smaller Rakai for plotting (not current Rakai region no longer includes kabula subdistrict 3)
    # minirak<-rakai3[which(rakai3$NAME_2!="Kabula" | rakai3$NAME_3=="Lyantonde Tc" | rakai3$NAME_3=="Lyantonde"),]
    minirak <- rakai3[which(rakai3$NAME_2 != "Kabula" | rakai3$NAME_3 == "Lyantonde Tc"), ]

    ##################################
    # set the coordinates of the data #
    ##################################

    DT <- DT[!is.na(LATITUDE_JITTER) & !is.na(LONGITUDE_JITTER)]
    coordinates(DT) <- ~ LONGITUDE_JITTER + LATITUDE_JITTER
    # set coordinate system to match uganda files
    proj4string(DT) <- proj4string(uganda1)

    # convert to m in order to build a 30x30m grid
    newcrs <- CRS("+proj=robin +datum=WGS84")
    DTnew <- spTransform(DT, newcrs)
    rakai1trans <- spTransform(rakai1, newcrs)
    miniraktrans <- spTransform(minirak, newcrs)
    masaka1trans <- spTransform(masaka1, newcrs)

    #######################################
    # Combine rakai1trans and masaka1trans #
    #######################################

    outline <- union(rakai1trans, masaka1trans)
    # find the extent of the data
    exnew <- extent(DTnew)
    # extent of the maps
    exmap <- extent(outline)

    # chose extent to cover all the data and rakai district

    # With a 30m grid, I think the same individuals are usually entering calculations for a large number of grid points
    # Do we really need a 30m grid? Why not 100m?

    grid <- raster(xmn = min(exnew[1], exmap[1]), xmx = exnew[2], ymn = exmap[3], ymx = exnew[4], res = 100)
    grid[] <- 1:ncell(grid) # No longer needed

    # set the coordinate reference system to match
    proj4string(grid) <- proj4string(DTnew)

    # restrict grid to map
    gridmask <- mask(grid, outline) # Restrict the map after
    # plot(gridmask)

    # consider the grid points in a data frame
    id <- as.data.table(1:ncell(gridmask))
    setnames(id, "V1", "ID")
    griddf <- as.data.table(SpatialPoints(grid))
    griddf <- data.table(id, griddf)
    setnames(griddf, gsub("y", "LAT_GRID", gsub("x", "LONG_GRID", colnames(griddf))))

    bw <- 3000
    bw2 <- bw * bw
    # require(mvtnorm)
    # dmvnorm( c(3.84,0) )	# ~ 9.996634e-05
    threshold <- bw * 3.84 # cut if density is < 1e-4
    threshold <- threshold * threshold # square the threshold, to avoid sqrt calculations in loop
    norm.const <- 1 / (2 * pi * bw2)

    tmp <- griddf[1:1e4, ]
    anst <- system.time({
        ans <- tmp[,
            {
                z1 <- LONG_GRID - DTnew@coords[, "LONGITUDE_JITTER"]
                z2 <- LAT_GRID - DTnew@coords[, "LATITUDE_JITTER"]
                z1 <- z1 * z1 + z2 * z2 # square distance
                z2 <- which(z1 < threshold) # avoid sqrt on 2e4 entries
                w <- norm.const * exp(-0.5 * z1 / bw2) # now with correct normalising constant
                # 	Xiayue
                # z3 <-  z1*z1 + z2*z2
                # z4 <- which(z3<threshold)
                # z <- cbind(matrix(z1[z4],ncol=1),matrix(z2[z4],ncol=1))
                # OR: the source code in Boom seems quite slow, with Cholesky decomposition etc. DIY faster?
                # w <- dmvn(z,mu=c(0,0),bw^2*diag(2))
                # z2 <- z4
                # 	olli
                # z1	<- z1*z1 + z2*z2 		# square distance
                # z2	<- which(z1<threshold)	# avoid sqrt on 2e4 entries
                # # to avoid very large output data, calculate directly all smooths here
                # z1	<- sqrt(z1[z2])			# sqrt on few entries
                # w	<- dnorm(z1, mean=0, sd=bw) # OR: I agree the normalising constant is not right
                # code assumes @coords and @data has same order.
                list(
                    HIV_STATUS_MEAN = mean(DTnew@data$HIV_STATUS[z2]), # no weighting by distance
                    HIV_STATUS_KERNEL = sum(DTnew@data$HIV_STATUS[z2] * w) / sum(w), # Gaussian kernel
                    VL_COPIES_KERNEL_GEOMMEAN = exp(sum(w * log(DTnew@data$VL_COPIES[z2] + 1)) / sum(w)) - 1, # Geometric Mean Kernel
                    VL_DETECTABLE_KERNEL = sum(DTnew@data$VL_DETECTABLE[z2] * w) / sum(w) # Detectable Prevelance
                )
            },
            by = c("ID", "LONG_GRID", "LAT_GRID")
        ]
    })

    grid[] <- ans[, VL_DETECTABLE_KERNEL]
    gridmask <- mask(grid, outline)

    # Breaks chosen by looking at data - need refining
    plot(gridmask,
        breaks = c(0, 0.025, 0.05, 0.075, 0.1, 0.5),
        col = brewer.pal(11, "RdYlGn")[c(10, 9, 5, 4, 3)],
        axes = FALSE, box = FALSE, ylim = c(exmap[3], -6000), legend = FALSE
    )

    plot(outline, add = TRUE)
    par(xpd = TRUE)
    legend("right",
        legend = c("0-2.5", "2.5-5", "5-7.5", "7.5-10", ">10"),
        fill = brewer.pal(11, "RdYlGn")[c(10, 9, 5, 4, 3)],
        horiz = FALSE, inset = -0.175,
        title = "Prevelence of \n Detectable \n Viremia (%)",
        cex = 0.8, box.lty = 0
    )

    grid[] <- ans[, VL_COPIES_KERNEL_GEOMMEAN]
    gridmask <- mask(grid, outline)
    plot(gridmask, breaks = c(0, 0.8, 1.5, 2.5, 3, 145), col = brewer.pal(11, "RdYlGn")[c(10, 9, 5, 4, 3)], axes = FALSE, box = FALSE, ylim = c(exmap[3], -6000), legend = FALSE)
    plot(outline, add = TRUE)
    par(xpd = TRUE)

    legend("right",
        legend = c("0-0.8", "0.8-1.5", "1.5-2.5", "2.5-3", ">3"),
        fill = brewer.pal(11, "RdYlGn")[c(10, 9, 5, 4, 3)],
        horiz = FALSE, inset = -0.175,
        title = "Geometric Mean \n VL (Copies/ml)",
        cex = 0.8, box.lty = 0
    )
}

glm_random_effects <- function(formula_glm) {
    glm_reffs <- glmer(
        data = dglm,
        formula = formula_glm,
        weights = N_HIV,
        control = glmerControl(optimizer = "Nelder_Mead"),
        family = "binomial"
    )
    ilink <- family(glm_reffs)$linkinv

    dpreds <- copy(dglm)
    dpreds[, ROUND := as.factor(ROUND)]
    dpreds[, PRED := ilink(predict(glm_reffs))]
    dpreds[, `:=`(
        CL = qbinom(p = .25, size = N_HIV, prob = PRED) / N_HIV,
        CU = qbinom(p = .975, size = N_HIV, prob = PRED) / N_HIV
    ), by = PRED]

    cat(
        "The proportion of observations falling within prediction intervals are:\n",
        dpreds[, mean(PVLNSofHIV_MEAN <= CU & PVLNSofHIV_MEAN >= CL)], "\n"
    )

    g <- ggplot(dpreds, aes(x = ordered(COMM_NUM, levels = comm_lvls), color = ROUND, pch = FC2)) +
        geom_linerange(aes(ymin = CL, ymax = CU), position = position_dodge(width = .4)) +
        geom_point(aes(y = PVLNSofHIV_MEAN, size = N_HIV), position = position_dodge(width = .4)) +
        geom_vline(xintercept = 6.5, linetype = "dotted") +
        geom_vline(xintercept = 36.5, linetype = "dotted") +
        facet_grid(SEX ~ .) +
        theme_bw()
    print(g)

    glm_reffs
}

glm_no_random_effects <- function(formula_glm) {
    glm_base <- glm(
        data = dglm,
        formula = formula_glm,
        weights = N_HIV,
        family = "binomial"
    )
    ilink <- family(glm_base)$linkinv

    dpreds <- copy(dglm)
    dpreds[, ROUND := as.factor(ROUND)]
    dpreds[, PRED := ilink(predict(glm_base))]
    dpreds[, `:=`(
        CL = qbinom(p = .25, size = N_HIV, prob = PRED) / N_HIV,
        CU = qbinom(p = .975, size = N_HIV, prob = PRED) / N_HIV
    ), by = PRED]

    cat(
        "The proportion of observations falling within prediction intervals are:\n",
        dpreds[, mean(PVLNSofHIV_MEAN <= CU & PVLNSofHIV_MEAN >= CL)], "\n"
    )

    g <- ggplot(dpreds, aes(x = ordered(COMM_NUM, levels = comm_lvls), pch = FC2, color = ROUND)) +
        geom_linerange(aes(ymin = CL, ymax = CU), position = position_dodge(width = .4)) +
        geom_point(aes(y = PVLNSofHIV_MEAN, size = N_HIV), position = position_dodge(width = .4)) +
        geom_vline(xintercept = 6.5, linetype = "dotted") +
        geom_vline(xintercept = 36.5, linetype = "dotted") +
        facet_grid(SEX ~ .) +
        theme_bw()
    print(g)

    glm_base
}

glm_random_effects_stan <- function(formula_glm) {
    glm_reffs <- stan_glmer(
        data = dglm,
        formula = formula_glm,
        weights = N_HIV,
        control = glmerControl(optimizer = "Nelder_Mead"),
        family = "binomial"
    )
    ilink <- family(glm_reffs)$linkinv

    dpreds <- copy(dglm)
    dpreds[, ROUND := as.factor(ROUND)]
    dpreds[, PRED := ilink(predict(glm_reffs))]
    dpreds[, `:=`(
        CL = qbinom(p = .25, size = N_HIV, prob = PRED) / N_HIV,
        CU = qbinom(p = .975, size = N_HIV, prob = PRED) / N_HIV
    ), by = PRED]

    cat(
        "The proportion of observations falling within prediction intervals are:\n",
        dpreds[, mean(PVLNSofHIV_MEAN <= CU & PVLNSofHIV_MEAN >= CL)], "\n"
    )

    g <- ggplot(dpreds, aes(x = ordered(COMM_NUM, levels = comm_lvls), color = ROUND, pch = FC2)) +
        geom_linerange(aes(ymin = CL, ymax = CU), position = position_dodge(width = .4)) +
        geom_point(aes(y = PVLNSofHIV_MEAN, size = N_HIV), position = position_dodge(width = .4)) +
        geom_vline(xintercept = 6.5, linetype = "dotted") +
        geom_vline(xintercept = 36.5, linetype = "dotted") +
        facet_grid(SEX ~ .) +
        theme_bw()
    print(g)

    glm_reffs
}


mcmc_intervals_2 <- function(x,
                             pars = character(),
                             regex_pars = character(),
                             transformations = list(),
                             ...,
                             prob = 0.5,
                             prob_outer = 0.9,
                             point_est = c("median", "mean", "none"),
                             outer_size = 0.5,
                             inner_size = 2,
                             point_size = 4,
                             rhat = numeric(),
                             color = NA) {
    .gs <- function(x) as.integer(gsub("[A-z]|\\)|\\(|:", "", x))
    p_mint$data$COMM_NUM <- .gs(p$data$parameter)
    tmp <- unique(dglm[, .(COMM_NUM, FC2), ])
    p_mint$data <- merge(p$data, tmp)

    out <- ggplot(
        p_mint$data,
        aes(x = m, y = parameter, xmin = ll, xmax = hh, col = FC2)
    ) +
        geom_point(size = 1.5) +
        geom_errorbarh(height = 0, size = 0.5) +
        geom_errorbarh(aes(xmin = l, xmax = h), height = 0, size = 1)

    out
}

plot.first.participant.estimates <- function(type="supp-pop"){

    type <- match.arg(type, c("supp-pop", "prevalence", 'supp-hiv'))
    if( type == 'supp-hiv'){
        .expr <- expr(Hmisc::binconf(x=sum(VL_COPIES[z] >= VIREMIC_VIRAL_LOAD), n=sum(HIV_STATUS[z] == 1), return.df =T ))
        .ylab <- 'Proportion of PLHIV with unsuppressed VL'
    }else if(type == 'prevalence'){
        .expr <- expr(Hmisc::binconf(x=sum(HIV_STATUS==1, n=length(HIV_STATUS)), return.df =T ))
        .ylab <- 'Prevalence of HIV'
    }else{
        .expr <- expr(Hmisc::binconf(x=sum(VL_COPIES[z] >= VIREMIC_VIRAL_LOAD), n=length(HIV_STATUS), return.df =T))
        .ylab <- 'Proportion of pop with unsuppressed VL'
    }

    tmp <- copy(dall)
    tmp[, AGEGROUP := split.agegroup(AGEYRS)]
    tmp <- tmp[ FIRST_PARTICIPATION == 1 & ROUND %in% c(16 ,19), {
        z <- ! is.na(VL_COPIES); eval(.expr)
    }, by=.(FC,ROUND, SEX, AGEGROUP)] 
    prettify_labels(tmp)
    ggplot(tmp, aes(x=AGEGROUP, y=PointEst, ymin=Lower, ymax=Upper, color=SEX_LAB)) +
        geom_point(position = position_dodge(width = .4)) +
        geom_linerange(position = position_dodge(width = .4)) +
        facet_grid(FC_LAB ~ ROUND_LAB) +
        scale_color_manual(values=palettes$sex, labels=sex_dictionary2) +
        my_labs(y=.ylab)
}


.reconstruct_ppDT <- function(standata){
    nms <- which(lapply(standata, length) == 35) |> names()
    DT <- as.data.table(standata[nms])
    DT[, `:=` (AGE_LABEL = observed_idx/2 + 14, observed_idx = NULL)]
    DT <- melt(DT, id.vars="AGE_LABEL")
    DT[, `:=` (
        TYPE = fifelse(variable %like% "y_observed", yes="y", no="N"),
        PTYPE = fifelse(variable %like% "ftp", yes="ftp", no="all"),
        SEX_LABEL = sub("^.*([0-1])([0-1]).*$", "\\1", variable),
        LOC_LABEL = sub("^.*([0-1])([0-1]).*$", "\\2", variable),
        variable = NULL
    )]
    DT <- dcast(DT, AGE_LABEL + SEX_LABEL + LOC_LABEL + PTYPE ~ TYPE, value.var = "value" )

    cols <- c("M", "CU", "CL")
    return(DT[, (cols) := binconf(y, N, return.df = T)])
}

plot_single_posterior_fit <- function(fit_rds, standata_rds, model, round){

    fit <- readRDS(fit_rds)
    stan.data <- readRDS(standata_rds)

    re <- fit$draws(format="df") |> as.data.table()
    names_vars <- dimnames(re)[[2]]

    # reconstruct PT
    ppDT <- .reconstruct_ppDT(standata=stan.data)

    # make plots

    q <- c("M"=.5, "CL"=.025, "CU"=.975)
    cols <- names(re) %which.like% '^p_predict_'
    tmp <- re[, lapply(.SD, posterior::quantile2, probs=q), .SDcols =cols]
    tmp[, quantile := names(q)]
    tmp <- melt( tmp, id.vars = "quantile")

    prop.by.age.sex.loc.type <- dcast.data.table(tmp, variable ~ quantile, value.var = "value")
    prop.by.age.sex.loc.type[ , `:=` (
        AGE_LABEL = .stan.brackets.to.age(variable, .stan.data=stan.data),
        variable = .stan.remove.brackets(variable)
    ) ]
    tmp1 <- .stan.get.sex.and.loc(prop.by.age.sex.loc.type, 'variable')

    prop.by.age.sex.loc.type <- merge(x=tmp1, y=group_codes, by=c('SEX', 'LOC'))
    # .stan.get.sex.and.loc(prop.by.age.sex.loc.type, 'variable', codes=unique(group_codes[, -"PTYPE"]))
    # prop.by.age.sex.loc.type <- .stan.get.sex.and.loc(nspop.byage.ptype, 'variable', codes=group_codes)

    # plot settings specific to model type
    .ylims <- if(model %like% "prevl"){
        c(0, .75) 
    }else if (model %like% "supp-hiv"){
        c(0, 1)
    } else {
        c(0,.5)
    }

    prettify_labels(ppDT)
    prettify_labels(prop.by.age.sex.loc.type)

    plot_ftp <- .plot.stan.fit(
        prop.by.age.sex.loc.type[PTYPE == "ftp"],
        DT2=ppDT[PTYPE == "ftp"],
        ylims = .ylims,
        ylab = NULL
    )[[2]]
    plot_all <- .plot.stan.fit(
        prop.by.age.sex.loc.type[PTYPE == "all"],
        DT2=ppDT[PTYPE == "all"],
        ylims = .ylims,
        ylab = NULL
    )[[2]]

    force(plot_ftp); force(plot_all)
    gg_list[[paste0(model, "_", round, "_ftp")]] <- plot_ftp
    gg_list[[paste0(model, "_", round, "_all")]] <- plot_all
    gg_list <<- gg_list

    list(ftp=plot_ftp, all=plot_all)
}
