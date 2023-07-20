module load anaconda3/personal
if [ -d $HOME/anaconda3 ]; then
    echo -e "\nanaconda3/personal already installed\n\n";
else
    echo -e "\ninstalling anaconda3/personal\n\n";
    anaconda-setup
fi


# Create new conda environment called "longivl_cmdstan"
if [ -d $HOME/anaconda3/envs/longivl_cmdstan ]; then
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
conda install -c conda-forge compilers r r-base r-essentials r-devtools


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
    "data.table", "dplyr", 
    "lubridate", "haven", "Hmisc",
    "DiagrammeR", "DiagrammeRsvg",
    "loo", "tidybayes", "bayestestR", "bayesplot", 
    "rsvg", "htmltools", "raster", "rnaturalearth", "osmdata", "sf", "foreach", "patchwork", "knitr", "foreign", "nnet", "stargazer", "tint", "scales", "geomtextpath","ggthemes", "gganimate", "ggtext", "ggpubr", "readxl", "rgdal", "rgeos", "RColorBrewer", "mvtnorm", "lme4", "optparse", "bh", "MCMCpack", "here", "yaml"
))
'

echo "=========================================\n\n
longivl_cmdstan: completed installation.\n
For next steps see\n
https://github.com/TODO
"
