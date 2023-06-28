module load anaconda3/personal
if [ -d $HOME/anaconda3 ]; then
    echo -e "\nanaconda3/personal already installed\n\n";
else
    echo -e "\ninstalling anaconda3/personal\n\n";
    anaconda-setup
fi


# Create new conda environment called "longi_vl_cmdstan"
if [ -d $HOME/anaconda3/envs/longi_vl_cmdstan ]; then
    echo "###############################################"
    echo -e "\nlongi_vl_cmdstan conda environment is already present"
    echo -e "\nIf you wish you re-install please remove the conda environment first with:"
    echo -e "\tconda remove -n longi_vl_cmdstan --all -y"
    echo -e "\n\n###############################################"
    exit 1
else
    echo -e "\nCreating Conda environment for R packages: longi_vl_cmdstan"
    conda create -n longi_vl_cmdstan -y
    source activate longi_vl_cmdstan
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
options(unzip = "internal"); 
devtools::install_github("stan-dev/cmdstanr"))
insall.packages(
c(
    "data.table", "dplyr", 
    "lubridate", "haven", "Hmisc",
    "DiagrammeR", "DiagrammeRsvg",
    "loo", "tidybayes", "bayestestR", "bayesplot", 
    "rsvg", "htmltools", "raster", "rnaturalearth", "osmdata", "sf", "foreach", "patchwork", "knitr", "foreign", "nnet", "stargazer", "tint", "scales", "ggthemes", "gganimate", "ggtext", "ggpubr", "readxl", "rgdal", "rgeos", "RColorBrewer", "mvtnorm", "lme4", "optparse", "bh", "mcmcpack", "here"
))
'

echo "=========================================\n\n
longi_vl_cmdstan: completed installation.\n
For next steps see\n
https://github.com/TODO
"
