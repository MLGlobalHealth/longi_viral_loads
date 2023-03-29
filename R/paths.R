gitdir.stan <- file.path(gitdir, 'stan')
gitdir.functions <- file.path(gitdir, 'functions')
gitdir.scripts <- file.path(gitdir, 'scripts')

# simplify sourcing of all R/*.R helpers
gitdir.R <- file.path(gitdir, 'R')
R_scripts <- list.files(
    file.path(gitdir, 'R'),
    pattern = "\\.R$", full.names = TRUE)

usr <- Sys.info()[['user']]

if(usr=='andrea')
{
    indir.deepdata <- '/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live/'
    # indir.deepsequence_analyses <- '/home/andrea/HPC/project/ratmann_xiaoyue_jrssc2022_analyses/live/PANGEA2_RCCS1519_UVRI'
}

indir.deepdata.r1520 <- file.path(indir.deepdata, 'RCCS_R15_R20')
indir.deepdata.r1518 <- file.path(indir.deepdata, 'RCCS_R15_R18')
indir.deepdata.r0914 <- file.path(indir.deepdata, 'RCCS_R9_R14')

# primary files:
path.quest_r1520_220830 <- file.path(indir.deepdata.r1520, 'Quest_R015_R020_VOIs_August302022.csv')
path.quest_r1519_221207 <- file.path(indir.deepdata.r1520, 'quest_R15_R19_VoIs_Dec072022.csv')
path.community.types <- file.path(indir.deepdata.r1518, 'community_names.csv')
path.hivstatusvl.r1520 <- file.path(indir.deepdata.r1520, 'all_participants_hivstatus_vl_220729.csv')
path.negatives.r1520 <- file.path(indir.deepdata.r1520, 'R016_R020_Data_for_HIVnegatives.csv')
path.hivres.r0914 <- file.path(indir.deepdata.r0914, "HIV_R09_R14.csv")
path.flow.r19 <- file.path(indir.deepdata.r1520, 'flowR19_VOIs.dta')

# processed files:
path.processed.hivstatus.r0920 <- file.path(indir.deepdata.r1520,"RCCS_processed_participants_hivstatus_230328.rds")
path.processed.testing.r1519 <- file.path(indir.deepdata.r1520,"RCCS_processed_participants_testing_230328.rds")

