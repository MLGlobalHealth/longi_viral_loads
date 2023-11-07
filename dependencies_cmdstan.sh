#!/bin/sh

ANACONDA_PATH="$HOME/anaconda3"

# load Anaconda
if [ -d $ANACONDA_PATH ]; then
    echo -e "\nanaconda already installed\n\n";
else
    # authomatic installation of anaconda on Imperial HPC.
    echo -e "\ninstalling anaconda3/personal\n\n";
    anaconda-setup
fi


# Create new conda environment called "longivl_cmdstan"
if [ -d $ANACONDA_PATH/envs/longivl_cmdstan ]; then
    echo "###############################################"
    echo -e "\nlongivl_cmdstan conda environment is already present"
    echo -e "\nIf you wish you re-install please remove the conda environment first with:"
    echo -e "\tconda remove -n longivl_cmdstan --all -y"
    echo -e "\n\n###############################################"
    exit 1
else
    echo -e "\nCreating Conda environment for R packages: longivl_cmdstan"
    conda create -n longivl_cmdstan -y
    source activate longivl_cmdstan
    # export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
fi


# Install base R
echo "\n\n=========================================\n\n
Installing compilers and base R\n\n"
conda install -c conda-forge compilers r r-base r-essentials r-devtools # udunits2 libgdal


# Install initial dependencies
echo -e "\nInstalling Dependencies: R packages via conda"

echo "\n\n=========================================\n\n
Installing R dependencies\n\n"

R -e '
r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)
options(unzip = "internal"); 
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
install.packages(
c(
    "DiagrammeR",
    "DiagrammeRsvg",
    "Hmisc",
    "MCMCpack",
    "RColorBrewer",
    "bayesplot",
    "bayestestR",
    "bh",
    "data.table",
    "dplyr", 
    "foreach",
    "foreign",
    "geomtextpath",
    "gganimate",
    "ggpubr",
    "ggtext",
    "ggthemes",
    # "ggsn",
    "haven",
    "here",
    "htmltools",
    "knitr",
    "lme4",
    "loo",
    "lubridate",
    "mvtnorm",
    "nnet",
    "optparse",
    "osmdata",
    "patchwork",
    "raster",
    "readxl",
    "rgdal",
    "rgeos",
    "rnaturalearth",
    "rsvg",
    "scales",
    "sf",
    "stargazer",
    "tidybayes",
    "tint",
    "yaml"
))
'

echo "=========================================\n\n
longivl_cmdstan: completed installation.\n
For next steps see\n
https://github.com/TODO
"
