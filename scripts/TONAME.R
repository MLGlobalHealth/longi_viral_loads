# AIMS:


################
# DEPENDENCIES #
################

library(data.table)
library(ggplot2)
library(foreign)
library(nnet)
library(stargazer)

################
#    PATHS     #
################

usr <- Sys.info()[['user']]
if(usr == 'andrea')
{

}else{

}


################
#    HELPERS   #
################


################
#     MAIN     #
################

# Load envs

# Look below for instructions:
# https://www.princeton.edu/~otorres/LogitR101.pdf 

dclass_16 <- dclass_16[, .(study_id, class)]
dclass_16
table(dclass_16[, class])

dvl[]

dvl_16[, uniqueN(comm) == 1, by='study_id'][, all(V1)]
tmp <- participant_comm2mode(dvl_16, quietly=TRUE)
tmp <- tmp[, .(comm=comm[1], sex=sex[1], lastvisit=hivdate[.N]), by='study_id']

tmp <- merge(tmp, dbirth, by=c('study_id', 'sex'))
tmp[, age_lastvisit := round((lastvisit-birthdate)/365.25)]
tmp[, `:=` (lastvisit=NULL, birthdate=NULL)]
tmp <- merge(dclass_16, tmp)

# refactor ARV as consistenly, never, interrupted?
dclass_16

darv_16 <- merge(dvl_16[, .(study_id, round)], darv, by=c('study_id', 'round'))
darv_16[, arv_class := NA_character_]
darv_16[is.na(arv_class),
        arv_class := ifelse(all(cuarvmed == TRUE),
                           'always',
                           NA),
     by='study_id']
darv_16[is.na(arv_class),
        arv_class := ifelse(all(cuarvmed == FALSE),
                           'never',
                           NA),
     by='study_id']

darv_16[is.na(arv_class),
        arv_class := ifelse(!is.unsorted(cuarvmed),
                        'eventually',
                        NA),
        by='study_id']

# So does that mean that no one interrupted treatment???
# it seems that only 4 individuals did not report ARV at the last visit.
# TODO: check that this is not a bug in my code
# For example, load crude data and merge to dvl_16...
idx <- darv_16[, cuarvmed[.N] == FALSE ,study_id][V1 == TRUE, study_id]
darv_16[study_id %in% idx]

stopifnot(darv_16[is.na(arv_class), .N == 0])
tmp1 <- unique(darv_16[, .(study_id, arv_class)])
tmp <- merge(tmp, tmp1)
tmp

tmp[, unique(class)]
tmp[, factor(class, 
             levels= c("durably_suppressed", "newly_suppressed", "intermittently_viremic", "newly_viremic", "durably_viremic")
             )]
model_multi0 <- multinom(data=tmp, 
                         class ~ sex + comm + arv_class)

summary(model_multi0)
stargazer(model_multi0, type='text')

# how would we check the quality of such a model?
# read as:
# "If you are (MALE) you are more likely to be in (DURABLY VIREMIC) than (DURABLY SUPPRESSED) as compared to a (FEMALE)" 

dclass_16[, table(class)]

dbirth
