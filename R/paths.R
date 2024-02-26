gitdir.stan      <- file.path(gitdir, 'stan')
gitdir.functions <- file.path(gitdir, 'functions')
gitdir.scripts   <- file.path(gitdir, 'scripts')
gitdir.src       <- file.path(gitdir, 'src')
gitdir.data      <- file.path(gitdir, 'data')
gitdir.R         <- file.path(gitdir, 'R')

#####################
# check indir paths #
#####################

indir.zenodo <- "PATH-TO-DOWNLOADED-ZENODO-DATA"
indir.deepdata <- "TODO-IF-YOU-HAVE-ACCESS-TO-CONFIDENTIAL-DATA"
indir.deepsequence_analyses <- "TODO-IF-YOU-HAVE-ACCESS-TO-CONFIDENTIAL-DATA"

usr <- Sys.info()[['user']]

if(usr=='andrea') {
    # Local machine
    indir.deepdata <- '/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live'
    indir.deepsequence_analyses <- '/home/andrea/HPC/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI'
    indir.zenodo <- "/home/andrea/OneDrive/Imperial/zenodo-longivl/"
    OUTDIR <- '/home/andrea/HPC/ab1820/home/projects/2022/longvl'
}else if(usr == 'ab1820') {
    # HPC 
    indir.deepdata <- '/rds/general/project/ratmann_pangea_deepsequencedata/live'
    indir.deepsequence_analyses <- '/rds/general/project/ratmann_xiaoyue_jrssc2022_analyses/live/PANGEA2_RCCS1519_UVRI'
    indir.zenodo <- "TODO"
}

# simplify sourcing of all R/*.R helpers
R_scripts <- list.files( gitdir.R, pattern = "\\.R$", full.names = TRUE)
R_scripts <- R_scripts[! R_scripts %like% 'test_palettes.R|local_cores_parallelisation.R|paths.R']
for (path in R_scripts)
source(file=path)

if( args$confidential ){
    dir.exists(c(indir.deepdata, indir.deepsequence_analyses)) |> 
        all() |> 
        stopifnot("Input data directories could not be found in paths.R"=_)
}

##############################
# proceed to all definitions #
##############################

if ( exists("indir.deepdata") && dir.exists(indir.deepdata) ){

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
    path.community.gps <- file.path(indir.deepdata.r1518, 'Rakai_community_geography_R15.rda')
    path.negatives.r1520 <- file.path(indir.deepdata.r1520, 'R016_R020_Data_for_HIVnegatives.csv')
    path.hivres.r0914 <- file.path(indir.deepdata.r0914, "HIV_R09_R14.csv")

    # for first participants: 
    path.participation <- file.path(indir.deepdata.r1520, "RCCS_1strndquest_participation.dta")

    # for get_census_eligible_count.R
    path.flow.r1518 <- file.path(indir.deepdata.r1518, 'FlowR15_R18_VoIs_221118.csv')
    path.flow.r19 <- file.path(indir.deepdata.r1520, 'flowR19_VOIs.dta')

    # Pre-processed files
    # ___________________

    # participants used in Nature Medicine paper
    # ABl: rincp from l.192 of get_table_characteristics_pop.R
    path.nm.participants <- file.path(indir.deepdata.r1518, "participants_for_Andrea.RDS")

    # obtained from process_data.R (also have a 230502 version)
    path.hivstatusvl.r1520.old2 <- file.path(indir.deepdata.r1520, 'all_participants_hivstatus_vl_220729.csv')
    path.hivstatusvl.r1520.old <- file.path(indir.deepdata.r1520, 'all_participants_hivstatus_vl_230515.csv')
    path.hivstatusvl.r1520 <- file.path(indir.deepdata.r1520,'all_participants_hivstatus_vl_240214.csv')
    path.viralloads.processed.r1520 <- file.path(indir.deepdata.r1520, 'viral_loads_r15r20_processed_230502.csv')

    # processed files: 
    path.processed.hivstatus.r0920 <- file.path(indir.deepdata.r1520,"RCCS_processed_participants_hivstatus_230328.rds")
    path.processed.testing.r1519 <- file.path(indir.deepdata.r1520,"RCCS_processed_participants_testing_230328.rds")
    path.census.eligible <- file.path(gitdir.data, 'census_eligible_individuals_230514.csv')

    # OLD NAMES
    # path.aggregated.nums.denoms.r1619 <- file.path(indir.deepdata.r1520,"RCCS_aggregated_nums_denoms.csv")
    # path.participation.rates <- file.path(gitdir.data, "participation_rates_230517.rds")
    # path.participation.rates <- file.path(gitdir.data,  "participation_rates_240214.rds")
    # path.census.eligible.aggregated <- file.path(gitdir.data, "RCCS_census_eligible_aggregated.csv")
    # path.comm.censsize <- file.path(gitdir.data, 'censsize_by_community.rds')
}


##################
# files to share #
##################

if (indir.zenodo != "PATH-TO-DOWNLOADED-ZENODO-DATA") {
    indir.zenodo.data <- file.path(indir.zenodo, 'data')
    indir.zenodo.results <- file.path(indir.zenodo, 'results')

    path.aggregated.nums.denoms.r1619 <- file.path(indir.zenodo.data, "RCCS_aggregated_nums_denoms.csv")
    path.participation.rates <- file.path(indir.zenodo.data,  "participation_rates_240214.rds")
    path.census.eligible.aggregated <- file.path(indir.zenodo.data,  "RCCS_census_eligible_aggregated.csv")
    path.comm.censsize <- file.path(indir.zenodo.data,  'censsize_by_community.rds')
}

#####################
#     stan args     #
#####################

path.stan.config <- file.path(gitdir.stan, 'binomial_gp_config.yml')
if( usr == "andrea" & interactive() ){
    path.stan.config  <- file.path(gitdir.stan, 'binomial_gp_config_local.yml')
}
