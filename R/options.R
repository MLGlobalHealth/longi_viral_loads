# Define all possible options here . 
# Not all scripts require the same exact options, so it is possible to subset the option_list as we want in each script. 

library(optparse)

option_list <- list(
    make_option(
        "--viremic-viral-load",
        type = "numeric",
        default = 1000,
        help = "Duration of acute phase, in years [Defaults to 2 months]", 
        dest = 'viremic.viral.load'
    ),
    make_option(
        "--vl-detectable",
        type = "numeric",
        default = 150,
        help = "Duration of acute phase, in years [Defaults to 2 months]", 
        dest = 'detectable.viral.load'
    ),
    make_option(
        "--outdir",
        type = "character",
        default = '/home/andrea/HPC/ab1820/home/projects/2022/longvl',
        help = "Path to output directory where to store results [Defaults to my directory]", 
        dest = 'out.dir.prefix'
    ),
    make_option(
        "--indir",
        type = "character",
        default = '/home/andrea/HPC/ab1820/home/projects/2022/longvl',
        help = "Path to input directory where results are stored [Defaults to my directory]", 
        dest = 'indir'
    ),
    make_option(
        "--round",
        type = "numeric",
        default = c(16,17,18,19),
        help = "Interview Rounds to include in the analysis.", 
        dest = 'round'
    ),
    make_option(
        "--jobname",
        type = "character",
        default = NA_character_,
        help = "Jobname used to indicate model run",
        dest= "jobname"
    )
)

args <-  parse_args(OptionParser(option_list = option_list))
