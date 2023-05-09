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

community_dictionary <- list(
    short = c(
        "F" = "Fishing",
        "I" = "Inland",
        "fishing" = "Fishing",
        "inland" = "Inland"
    ),
    long = c(
        "F" = "Fishing communities",
        "I" = "Inlan communitiesd",
        "fishing" = "Fishing communities",
        "inland" = "Inland communities"
    )
)

sex_dictionary <- c(M = "Male", F = "Female")

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
    LABS = paste0("Round ", 16:19, ":\n"),
    START = c("07/2013", "02/2015", "10/2016", "06/2018"),
    END = c("01/2015", "09/2016", "05/2018", "05/2019")
)
drounds[, LABS := paste0(LABS, START, " to ", END)]

round_labs <- drounds$LABS
names(round_labs) <- drounds$ROUND


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

# to work together with my_labs() function
my_labs_dictionary <- c(
    AGE = "Age",
    AGEYRS = "Age",
    AGE_LABEL = "Age",
    AGE_LAB = "Age",
    COMM = "Community type",
    COMM_NUM = "Community number",
    FC = "Community type",
    FC2 = "Community type",
    HIV = "HIV status",
    HIV_LAB = "HIV status",
    PHIV_MEAN='HIV prevalence',
    PVLNS_MEAN='Prevalence of viraemia among population',
    PVLNSofHIV_MEAN='Prevalence of viraemia among HIV infected',
    SEX = "Gender",
    SEX_LAB = "Gender",
    VL_COPIES = "Viral Load Copies",
    NULL
)
