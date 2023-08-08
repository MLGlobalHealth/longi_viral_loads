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
        "I" = "Inland communitiesd",
        "fishing" = "Fishing communities",
        "inland" = "Inland communities",
        "Fishing" = "Fishing communities",
        "Inland" = "Inland communities"
    ),
    longest = c(
        "F" = "Fishing communities with high HIV prevalence",
        "I" = "Inland communities with typical HIV prevalence",
        "fishing" = "Fishing communities with high HIV prevalence",
        "inland" = "Inland communities with typical HIV prevalence",
        "Fishing" = "Fishing communities with high HIV prevalence",
        "Inland" = "Inland communities with typical HIV prevalence"
    )

    # I would probably try "Fishing communities with high HIV prevalence", "Inland communities with typical/more moderate HIV prevalence"
)

sex_dictionary <- c(M = "Male", F = "Female", `0` = "Female", `1` = "Male", `Total`="Total")

sex_dictionary2 <- c(M="Men", F="Women", Male="Men", Female="Women", `0`="Women", `1`="Men", `Total`= "Total")

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
    "run-gp-prevl" = "HIV prevalence"
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
        "UNSUPP" = "unsuppressed viral load",
        "SUPP" = "suppressed viral load",
        "NEG" = "HIV negative"
    ),
    PARTICIPATION_STATUS = c(
        "IN" = "in-study",
        "OUT" = "out-of-study"
    )
)

model_dict <- c(
    `run-gp-prevl` = "HIV prevalence",
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
    FC2 = "Community type",
    FTP_LAB = "Participants type",
    HIV = "HIV status",
    HIV_LAB = "HIV status",
    LOC = "Community type",
    LOC_LAB = "Community type",
    LOESS_SPAN = "Loess span",
    PARTRATE = "Participation rate",
    PART_RATE = "Participation rate",
    PHIV_MEAN = "HIV prevalence",
    PVLNS_MEAN = "Prevalence of viraemia among population",
    PVLNSofHIV_MEAN = "Prevalence of viraemia among HIV infected",
    ROUND = "Interview round",
    ROUND_LAB = "Interview round",
    ROUND_LABEL = "Interview round",
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
