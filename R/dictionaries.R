question_dict <- list(
    rhivever2= 'Have you ever received your HIV results from anywhere?',
    RHIVEVER2= 'Have you ever received your HIV results from anywhere?',
    hivperiod2= 'How long ago did you last receive your HIV results?(years)',
    HIVPERIOD2= 'How long ago did you last receive your HIV results?(years)',
    hivrslt2= 'What was the results of this last HIV test?',
    HIVRSLT2= 'What was the results of this last HIV test?',
    oghivts2= 'From whom did you receive the test?',
    OGHIVTS2= 'From whom did you receive the test?'
)

answers_dict <- list(
    rhivever = c(
        '1'='1.Yes',
        '2'='2.No',
        '8'='8.Not Applicable',
        '9'='9.NR'
    ),
    hivperiod= c(
        '1'= '1.(<1)',
        '2'= '2.(1-2)',
        '3'= '3.(3-4)',
        '4'= '4.(>4)',
        '8'= 'NA?'),
    hivrslt = c(
        '1'='1.Negative',
        '2'='2.Positive',
        '3'='3.Indeterminate',
        '7'="7.Don't know",
        '8'='8. ?',
        '9'='9.NR'),
    oghivts1 = c(
        '1'='1. Rakai Project',
        '2'='2. RAIN',
        '3'='3. Other NGO (specify)',
        '4'='4. Gov. Doctor/nurse',
        '5'='5. Priv. Doctor/nurse/HW',
        '6'='6. Other (specify)',
        '7'='7.97.98.998. ?',
        '88'='88. No additional response',
        '97'='7.97.98.998. ?',
        '98'='7.97.98.998. ?',
        '998'='7.97.98.998. ?'
    )
)

sex_dictionary <- c(M='Male', F='Female')


round_dictionary = c(
    'R016'="Round 16", 
    'R017'="Round 17", 
    'R018'="Round 18", 
    'R019'="Round 19", 
    'R020'="Round 20", 
    '16'="Round 16", 
    '17'="Round 17", 
    '18'="Round 18", 
    '19'="Round 19", 
    '20'="Round 20", 
    'AGGR' = "Aggregate over rounds"
)

dall_dictionaries <- list( 
    TEST_YEAR_AGO=c(
        '0.5' = "<1",
        '1.5' = "1-2",
        '3' = "3-4",
        '5' = '>5',
        `NA`= 'never'),
    TEST_EVER = c(),
    TEST_LAST_RESULT = c(),
)

# extract_from_dall_dicts <- function(dict, keys)
# {
#     NA2textNA <- function(x)
#     {
#         x[is.na(x)] <- "NA"
#         return(x)
#     }
# 
#     keys <- as.character(keys) |> NA2textNA()
#     # return(keys)
#     return(dall_dictionaries[[dict]][keys])
# }

