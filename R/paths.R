gitdir.stan <- file.path(gitdir, 'stan')
gitdir.functions <- file.path(gitdir, 'functions')
gitdir.scripts <- file.path(gitdir, 'scripts')
gitdir.src <- file.path(gitdir, 'src')
gitdir.data <- file.path(gitdir, 'data')
gitdir.R <- file.path(gitdir, 'R')

#####################
# check indir paths #
#####################

usr <- Sys.info()[['user']]

if(usr=='andrea') {
    # Local machine
    indir.deepdata <- '/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live'
    indir.deepsequence_analyses <- '/home/andrea/HPC/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI'
    OUTDIR <- '/home/andrea/HPC/ab1820/home/projects/2022/longvl'
}else if(usr == 'ab1820') {
    # HPC 
    indir.deepdata <- '/rds/general/project/ratmann_pangea_deepsequencedata/live'
    indir.deepsequence_analyses <- '/rds/general/project/ratmann_xiaoyue_jrssc2022_analyses/live/PANGEA2_RCCS1519_UVRI'
}else {
    stop("Need to specify input and output directory in R/paths.R")
}

dir.exists(c(indir.deepdata, indir.deepsequence_analyses)) |> 
    all() |> 
    stopifnot("Input data directories could not be found in paths.R"=_)


# simplify sourcing of all R/*.R helpers
R_scripts <- list.files( gitdir.R, pattern = "\\.R$", full.names = TRUE)
R_scripts <- R_scripts[! R_scripts %like% 'local_cores_parallelisation.R|paths.R']
for (path in R_scripts)
    source(file=path)


##############################
# proceed to all definitions #
##############################

indir.deepdata.r1520 <- file.path(indir.deepdata, 'RCCS_R15_R20')
indir.deepdata.r1518 <- file.path(indir.deepdata, 'RCCS_R15_R18')
indir.deepdata.r0914 <- file.path(indir.deepdata, 'RCCS_R9_R14')

# primary files:
# ______________

path.viral.loads1 <- file.path(indir.deepdata.r1520, "Quest_R015_R020_VOIs_May062022.csv")
path.viral.loads2 <- file.path(indir.deepdata.r1520, "Allpcr_data_for_R015_R020_study_ids.xlsx")
path.quest_r1520_220830 <- file.path(indir.deepdata.r1520, 'Quest_R015_R020_VOIs_August302022.csv')
path.quest_r1519_221207 <- file.path(indir.deepdata.r1520, 'quest_R15_R19_VoIs_Dec072022.csv')
path.community.types <- file.path(indir.deepdata.r1518, 'community_names.csv')
path.community.idx <- file.path(indir.deepdata.r1518, 'community_id_index.csv')
path.negatives.r1520 <- file.path(indir.deepdata.r1520, 'R016_R020_Data_for_HIVnegatives.csv')
path.hivres.r0914 <- file.path(indir.deepdata.r0914, "HIV_R09_R14.csv")

# for first participants: 
path.participation <- file.path(indir.deepdata.r1520, "RCCS_1strndquest_participation.dta")

# for get_census_eligible_count.R
path.flow.r1518 <- file.path(indir.deepdata.r1518, 'FlowR15_R18_VoIs_221118.csv')
path.flow.r19 <- file.path(indir.deepdata.r1520, 'flowR19_VOIs.dta')

# Pre-processed files
# ___________________

# obtained from process_data.R (also have a 230502 version)
path.hivstatusvl.r1520 <- file.path(indir.deepdata.r1520, 'all_participants_hivstatus_vl_230515.csv')
path.viralloads.processed.r1520 <- file.path(indir.deepdata.r1520, 'viral_loads_r15r20_processed_230502.csv')

# processed files: 
path.processed.hivstatus.r0920 <- file.path(indir.deepdata.r1520,"RCCS_processed_participants_hivstatus_230328.rds")
path.processed.testing.r1519 <- file.path(indir.deepdata.r1520,"RCCS_processed_participants_testing_230328.rds")

# I believe this comes from melodie?
# path.census.eligible <- file.path(gitdir.data, 'RCCS_census_eligible_individuals_221209.csv')
path.census.eligible <- file.path(gitdir.data, 'census_eligible_individuals_230514.csv')

path.participation.rates <- file.path(gitdir.data, "participation_rates_230517.rds")
