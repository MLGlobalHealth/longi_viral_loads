question_dict <- list(
    rhivever2 = "Have you ever received your HIV results from anywhere?",
    RHIVEVER2 = "Have you ever received your HIV results from anywhere?",
    hivperiod2 = "How long ago did you last receive your HIV results?(years)",
    HIVPERIOD2 = "How long ago did you last receive your HIV results?(years)",
    hivrslt2 = "What was the results of this last HIV test?",
    HIVRSLT2 = "What was the results of this last HIV test?",
    oghivts2 = "From whom did you receive the test?",
    OGHIVTS2 = "From whom did you receive the test?"
)

answers_dict <- list(
    rhivever = c(
        "1" = "1.Yes",
        "2" = "2.No",
        "8" = "8.Not Applicable",
        "9" = "9.NR"
    ),
    hivperiod = c(
        "1" = "1.(<1)",
        "2" = "2.(1-2)",
        "3" = "3.(3-4)",
        "4" = "4.(>4)",
        "8" = "NA?"
    ),
    hivrslt = c(
        "1" = "1.Negative",
        "2" = "2.Positive",
        "3" = "3.Indeterminate",
        "7" = "7.Don't know",
        "8" = "8. ?",
        "9" = "9.NR"
    ),
    oghivts1 = c(
        "1" = "1. Rakai Project",
        "2" = "2. RAIN",
        "3" = "3. Other NGO (specify)",
        "4" = "4. Gov. Doctor/nurse",
        "5" = "5. Priv. Doctor/nurse/HW",
        "6" = "6. Other (specify)",
        "7" = "7.97.98.998. ?",
        "88" = "88. No additional response",
        "97" = "7.97.98.998. ?",
        "98" = "7.97.98.998. ?",
        "998" = "7.97.98.998. ?"
    )
)

stan_dicts <- list(
    SEXtoINT = c(
        `F` = 0L,
        `M` = 1L,
        `Female` = 0L,
        `Male` = 1L,
        `female` = 0L,
        `male` = 1L
    ),
    INTtoSEX = c(
        `0` = "F",
        `1` = "M"
    ),
    LOCtoINT = c(
        `inland` = 0L,
        `fishing` = 1L
    ),
    INTtoLOC = c(
        `0` = "inland",
        `1` = "fishing"
    )
)

community_dictionary <- list(
    short = c(
        "F" = "Fishing",
        "I" = "Inland",
        "fishing" = "Fishing",
        "inland" = "Inland"
    ),
    long = c(
        "F" = "Fishing communities",
        "I" = "Inland communities",
        "fishing" = "Fishing communities",
        "inland" = "Inland communities",
        "Fishing" = "Fishing communities",
        "Inland" = "Inland communities",
        `Fishing, Round 16` = "Fishing, Round 16\nJul 2013 to Jan 2015" ,
        `Fishing, Round 19` = "Fishing, Round 19\nJun 2018 to May 2019" ,
        `Inland, Round 16` = "Inland, Round 16\nJul 2013 to Jan 2015",
        `Inland, Round 19` = "Inland, Round 19\nJun 2018 to May 2019"
    ),
    longest = c(
        "F" = "Fishing communities with high HIV seroprevalence",
        "I" = "Inland communities with typical HIV seroprevalence",
        "fishing" = "Fishing communities with high HIV seroprevalence",
        "inland" = "Inland communities with typical HIV seroprevalence",
        "Fishing" = "Fishing communities with high HIV seroprevalence",
        "Inland" = "Inland communities with typical HIV seroprevalence"
    ),
    longestn = c(
        "F" = "Fishing communities\nwith high HIV seroprevalence",
        "I" = "Inland communities\nwith typical HIV seroprevalence",
        "fishing" = "Fishing communities\nwith high HIV seroprevalence",
        "inland" = "Inland communities\nwith typical HIV seroprevalence",
        "Fishing" = "Fishing communities\nwith high HIV seroprevalence",
        "Inland" = "Inland communities\nwith typical HIV seroprevalence"
    ),
    longest2 = c(
        "I" = "Inland communities with HIV seroprevalence typical of rural and semi-urban East Africa",
        "F" = "Hyperendemic Lake Victoria fishing communities",
        "inland" = "Inland communities with HIV seroprevalence typical of rural and semi-urban East Africa",
        "fishing" = "Hyperendemic Lake Victoria fishing communities",
        "Inland" = "Inland communities with HIV seroprevalence typical of rural and semi-urban East Africa",
        "Fishing" = "Hyperendemic Lake Victoria fishing communities"
    ),
    longest2n = c(
        "I" = "Inland communities with HIV seroprevalence\ntypical of rural and semi-urban East Africa",
        "F" = "Hyperendemic Lake Victoria\nfishing communities",
        "inland" = "Inland communities with HIV seroprevalence typical\nof rural and semi-urban East Africa",
        "fishing" = "Hyperendemic Lake Victoria\nfishing communities",
        "Inland" = "Inland communities with HIV seroprevalence typical\nof rural and semi-urban East Africa",
        "Fishing" = "Hyperendemic Lake Victoria\nfishing communities"
    ),
    none = c(
        "F" = "",
        "I" = "",
        "fishing" = "",
        "inland" = "",
        "Fishing" = "",
        "Inland" = ""
    ),
    none2 = c(
        "F" = NULL,
        "I" = NULL,
        "fishing" = NULL,
        "inland" = NULL,
        "Fishing" = NULL,
        "Inland" = NULL
    )

    # I would probably try "Fishing communities with high HIV seroprevalence", "Inland communities with typical/more moderate HIV seroprevalence"
)

ptype_dict = c(
    "all" = "All participants",
    "ftp" = "First-time participants"
)

sex_dictionary <- c(M = "Male", F = "Female", `0` = "Female", `1` = "Male", `Total`="Total")

sex_dictionary2 <- c(
    M="Men",
    F="Women",
    Male="Men",
    Female="Women",
    `0`="Women",
    `1`="Men",
    `Total`= "Total",
    `Round 16, Male` = "Men, Round 16\nJul 2013 to Jan 2015" ,
    `Round 19, Male` = "Men, Round 19\nJun 2018 to May 2019" ,
    `Round 16, Female` = "Women, Round 16\nJul 2013 to Jan 2015",
    `Round 19, Female` = "Women, Round 19\nJun 2018 to May 2019"
)

loc_dictionary <- c(
    inland = "Inland",
    fishing = "Fishing",
    trading = "Trading",
    agrarian = "Agrarian",
    `0` = "Inland",
    `1` = "Fishing"
)

model_dictionary <- c(
    "run-gp-supp-hiv" = "Prevalence of viral suppression among HIV positive individuals",
    "run-gp-supp-pop" = "Prevalence of viraemia",
    "run-gp-prevl" = "HIV seroprevalence"
)

round_dictionary <- c(
    "R016" = "Round 16",
    "R017" = "Round 17",
    "R018" = "Round 18",
    "R019" = "Round 19",
    "R020" = "Round 20",
    "16" = "Round 16",
    "17" = "Round 17",
    "18" = "Round 18",
    "19" = "Round 19",
    "20" = "Round 20",
    "AGGR" = "Aggregate over rounds"
)

drounds <- data.table(
    ROUND = 16:19,
    LABS = paste0("Round ", 16:19, "\n"),
    START = c("07/2013", "02/2015", "10/2016", "06/2018"),
    END = c("01/2015", "09/2016", "05/2018", "05/2019"),
    START2 = c("Jul 2013", "Feb 2015", "Oct 2016", "Jun 2018"),
    END2 = c("Jan 2015", "Sep 2016", "May 2018", "May 2019")
)
drounds[, LABS := paste0(LABS, START2, " to ", END2)]
with(drounds, {
    round_labs <<- rep(LABS, 2)
    names(round_labs) <<- c(ROUND, paste("Round", ROUND))
    round_labs2 <<- gsub("\n", ":  ", round_labs)
})

round_labs3 <- gsub("Jul ", "07/", round_labs2) |>
    gsub("Jan ", "01/", x=_) |>
    gsub("Feb ", "02/", x=_) |>
    gsub("Mar ", "03/", x=_) |>
    gsub("Apr ", "04/", x=_) |>
    gsub("May ", "05/", x=_) |>
    gsub("Jun ", "06/", x=_) |>
    gsub("Sep ", "09/", x=_) |>
    gsub("Oct ", "10/", x=_) |>
    gsub(" to ", "-", x=_) |>
    gsub("201", "1", x=_) |>
    gsub(":  ", "\n", x=_)

round_labs4 <- gsub("Jul ", "July ", round_labs2) |>
    gsub("Jan ", "January ", x=_) |>
    gsub("Feb ", "Februay ", x=_) |>
    gsub("Mar ", "March ", x=_) |>
    gsub("Apr ", "April ", x=_) |>
    gsub("May ", "May ", x=_) |>
    gsub("Jun ", "June ", x=_) |>
    gsub("Sep ", "September ", x=_) |>
    gsub("Oct ", "October ", x=_) |>
    gsub(":  ", "\n", x=_)

dall_dictionaries <- list(
    TEST_YEAR_AGO = c(
        "0.5" = "<1",
        "1.5" = "1-2",
        "3" = "3-4",
        "5" = ">5",
        `NA` = "never"
    ),
    TEST_EVER = c(),
    TEST_LAST_RESULT = c()
)

postproc_dictionaries <- list(
    VIR_STATUS = c(
        "UNSUPP" = "unsuppressed virus",
        "SUPP" = "suppressed virus",
        "NEG" = "HIV negative"
    ),
    PARTICIPATION_STATUS = c(
        "IN" = "in-study",
        "OUT" = "out-of-study"
    )
)

model_dict <- c(
    `run-gp-prevl` = "HIV seroprevalence",
    `run-gp-supp-hiv` = "Prevalence of suppression among HIV positive",
    `run-gp-supp-pop` = "Prevalence of viraemia"
)

# to work together with my_labs() function
my_labs_dictionary <- c(
    AGE = "Age",
    AGEYRS = "Age",
    AGEGROUP = "Age group",
    AGE_LABEL = "Age",
    AGE_LAB = "Age",
    COMM = "Community type",
    COMM_NUM = "Community number",
    FIRST_VISIT = "First participant",
    FC = "Community type",
    FC_LAB = "Community type",
    FC2 = "Community type",
    FTP_LAB = "Participants type",
    HIV = "HIV status",
    HIV_LAB = "HIV status",
    LOC = "Community type",
    LOC_LAB = "Community type",
    LOESS_SPAN = "Loess span",
    PARTRATE = "Participation rate",
    PART_RATE = "Participation rate",
    PHIV_MEAN = "HIV seroprevalence",
    PTYPE = "Participant type",
    PVLNS_MEAN = "Prevalence of viraemia among population",
    PVLNSofHIV_MEAN = "Prevalence of viraemia among HIV infected",
    ROUND = "Survey round",
    ROUND_LAB = "Survey round",
    ROUND_LABEL = "Survey round",
    SEX = "Gender",
    SEX_LAB = "Gender",
    SEX_LABEL = "Gender",
    VL_COPIES = "Viral Load Copies",
    NULL
)

dict_stan_params <- c(
    `alpha_00` = expression(paste("s.deviation:  ", alpha[0][0])),
    `alpha_01` = expression(paste("s.deviation:  ", alpha[0][1])),
    `alpha_10` = expression(paste("s.deviation:  ", alpha[1][0])),
    `alpha_11` = expression(paste("s.deviation:  ", alpha[1][1])),
    `rho_00` = expression(paste("lengthscale:  ", rho[0][0])),
    `rho_01` = expression(paste("lengthscale:  ", rho[0][1])),
    `rho_10` = expression(paste("lengthscale:  ", rho[1][0])),
    `rho_11` = expression(paste("lengthscale:  ", rho[1][1])),
    sex_0_loc_0 = expression(paste("baseline:  ", nu[0][0])),
    sex_0_loc_1 = expression(paste("baseline:  ", nu[0][1])),
    sex_1_loc_0 = expression(paste("baseline:  ", nu[1][0])),
    sex_1_loc_1 = expression(paste("baseline:  ", nu[1][1])),
    NULL
)

dict_stan_params2 <- with(stan_dicts, c(
    `alpha_00` = expression(paste("variance:  ", alpha["F"]["-fi"])),
    `alpha_01` = expression(paste("variance:  ", alpha["F"]["-in"])),
    `alpha_10` = expression(paste("variance:  ", alpha["M"]["-fi"])),
    `alpha_11` = expression(paste("variance:  ", alpha["M"]["-in"])),
    `rho_00` = expression(paste("lengthscale:  ",  rho["F"]["-fi"])),
    `rho_01` = expression(paste("lengthscale:  ",  rho["F"]["-in"])),
    `rho_10` = expression(paste("lengthscale:  ",  rho["M"]["-fi"])),
    `rho_11` = expression(paste("lengthscale:  ",  rho["M"]["-in"])),
    sex_0_loc_0 = expression(paste("baseline:  ",   nu["F"]["-fi"])),
    sex_0_loc_1 = expression(paste("baseline:  ",   nu["F"]["-in"])),
    sex_1_loc_0 = expression(paste("baseline:  ",   nu["M"]["-fi"])),
    sex_1_loc_1 = expression(paste("baseline:  ",   nu["M"]["-in"])),
    NULL
))

# same but with names = t(parname)
dict_stan_params_t <- copy(dict_stan_params)
names(dict_stan_params_t) <- paste0( "t(", names(dict_stan_params_t), ")")

plabels <- list(
    
    sex = c(
        M = "Men",
        Male = "Men",
        F = "Women",
        W = "Women",
        Female = "Women"
    )
)

dict_table_names <- list(
    percent_reduction_old = c(
        LOC_LAB = "Location",
        SEX_LAB = "Gender",
        ROUND_LAB = "Round",
        N_ELIGIBLE = "Number of census\neligible individuals",
        CELL_HIV_N = "Estimated number\nof PLHIV",
        CELL_HIV_P = "Estimated percent reduction\nin number of PLHIV",
        CELL_UNSUPP_N = "Estimated number\nof unsuppressed PLHIV",
        CELL_UNSUPP_P = "Estimated percent reduction\nin number of unsuppressed PLHIV",
        N_HIV = "Number of PLHIV (ftp)",
        N_UNSUPP = "Number of unsuppressed (ftp)",
        P_HIV_REDUCTION = "Percent reduction\nin number of PLHIV (ftp)",
        P_UNSUPP_REDUCTION = "Percent reduction\nin number of unsuppressed (ftp)"
    ),

    percent_reduction = c(
        LOC_LAB = "Location\n\n\n\n\n\n",
        SEX_LAB = "Gender\n\n\n\n\n\n",
        ROUND_LAB = "Round\n\n\n\n\n\n",
        N_ELIGIBLE = "Census-eligible\nindividuals\n\n\n\n(n)\n",
        CELL_HIV_N = "People with HIV \n\n\n\n\n(posterior median\nestimate, (95%CrI))",
        CELL_UNSUPP_N = "People with unsuppressed HIV\n\n\n\n\n(posterior median estimate,\n(95%CrI)) ",
        CELL_SUPPHIV = "Proportion of\npeople exhibiting viramia\namong PLHIV",
        VIR_P = "Proportion of\npeople \nexhibiting viraemia",
        CELL_UNSUPP_P = "Percent reduction in the absolute number \nof people with \nunsuppressed HIV relative to\nround 16\n\n(posterior median estimate,\n(95%CrI))",
        CELL_MF_GAP = "Male-female difference in the\nproportion of people\nwith HIV who have unsuppressed vireamia\n\n\n(posterior median estimate,\n(95%CrI))",
        NULL
    ),

    mean_ages = c(
        LOC_LAB = "Location",
        SEX_LAB = "Gender",
        ROUND_LAB = "Round",
        `AGEMEAN_run-gp-prevl` = "Mean age of PLHIV",
        `AGESTD_run-gp-prevl` = "Standard deviation of PLHIV",
        `AGEMEAN_run-gp-supp-pop` = "Mean age of unsuppressed population",
        `AGESTD_run-gp-supp-pop` = "Standard deviation of unsuppressed population",
        `AGE25_run-gp-prevl` = "First quartile for age of individuals with HIV",
        `AGE50_run-gp-prevl` = "Median age of individuals with HIV",
        `AGE75_run-gp-prevl` = "Third quartile for age of individuals with HIV",
        `AGE25_run-gp-supp-pop` = "First quartile for age of individuals with unsuppressed virus",
        `AGE50_run-gp-supp-pop` = "Median age of individuals with unsuppressed virus",
        `AGE75_run-gp-supp-pop` = "Third quartile for age of individuals with unsuppressed virus"
    ),

    eligible_participants = c(
        FC_LAB = "Location",
        SEX_LAB = "Gender",
        # ROUND_LAB = "Survey Round",
        ROUND= "Survey\nRound",
        ROUND_LAB = "Survey\nRound",
        EL = "Number of census eligible individuals",
        N_EL = "Number of census\neligible individuals",
        PART = "Number of participants",
        N_PAR = "Number of\nparticipants",
        PERC = "Participation\nrate",
        RATE = "Participation\nrate"
    ),

    first_time_participants = c(
        COMM_LAB = "Location",
        COMM = "Location",
        SEX_LAB = "Gender",
        # ROUND_LAB = "Survey Round",
        ROUND= "Survey\nRound",
        ROUND_LAB = "Survey\nRound",
        N = "First-time\nparticipants",
        N_HIV =  "First-time\nparticipants\nwith HIV",
        N_UV = "First-time\nparticipants\nwith\nunsuppressed virus",
        R_UV = "Ratio of\nfirst-time\nparticipants\nwih HIV\nrelative to\nround 16",
        R_HIV = "Ratio of\nfirst-time\nparticipants\nwith\nunsuppressed\nvirus\nrelative to\nround 16",
        NULL
    ),

    NULL
)
