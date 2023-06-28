# Define all possible options here . 
# Not all scripts require the same exact options, so it is possible to subset the option_list as we want in each script. 

library(optparse)

option_list <- list(
    make_option(
        "--refit",
        type = "logical",
        default = FALSE,
        help = "Flag on whether to re-run the stan models even if already exist [Defaults to FALSE]",
        dest = "refit"
    ),
    make_option(
        "--run-gp-prevl",
        type = "logical",
        default = FALSE,
        help = "Flag on whether to run the Stan model for prevalence [Defaults to FALSE] ",
        dest = "run.gp.prevl"
    ),
    make_option(
        "--run-icar-mean-vl",
        type = "logical",
        default = FALSE,
        help = "Flag on whether to run the Stan model for population level viral loads [Defaults to FALSE] ",
        dest = "run.icar.mean.vl"
    ),
    make_option(
        "--run-gp-supp-hiv",
        type = "logical",
        default = FALSE,
        help = "Flag on whether to run the Stan model for suppression among HIV positive [Defaults to FALSE] ",
        dest = "run.gp.supp.hiv"
    ),
    make_option(
        "--run-gp-supp-pop",
        type = "logical",
        default = FALSE,
        help = "Flag on whether to run the Stan model for suppression among population [Defaults to FALSE] ",
        dest = "run.gp.supp.pop"
    ),
    make_option(
        "--run-comm-analysis",
        type = "logical",
        default = FALSE,
        help = "Flag on whether to run the Stan GLM model for community-level suppression [Defaults to FALSE] ",
        dest = "run.comm.analysis"
    ),
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
        "--iter_warmup", 
        type = "integer", 
        default = 5e2,
        help = "HMC warmup iterations [default %default]",
        dest = 'iter.warmup'),
    make_option(
        "--iter_sampling",
        type = "integer",
        default = 10e3,
        help = "HMC of sampling iterations iterations [default %default]",
        dest = 'iter.sampling'),
    make_option(
        "--chains",
        type = "integer",
        default = 4,
        help = "Number of MCMC chains",
        dest = 'chains'),
    make_option(
        "--outdir-prefix",
        type = "character",
        default = '/home/andrea/HPC/ab1820/home/projects/2022/longvl',
        help = "Path to prefix output directory where to store results [Defaults to my directory]", 
        dest = 'out.dir.prefix'
    ),
    make_option(
        "--outdir",
        type = "character",
        default = NA_character_,
        help = "Path to exact output directory where to store results [Defaults to NA]", 
        dest = 'out.dir.exact'
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
        "--firstparticipants",
        type = "logical",
        default = FALSE,
        help = "Subsets the analysis to only first-time participants", 
        dest = 'only.firstparticipants'
    ),
    make_option(
        "--jobname",
        type = "character",
        default = NA_character_,
        help = "Jobname used to identify model run",
        dest= "jobname"
    )
)

args <-  parse_args(OptionParser(option_list = option_list))
