module load anaconda3/personal
if [ -d $HOME/anaconda3 ]; then
    echo -e "\nanaconda3/personal already installed\n\n";
else
    echo -e "\ninstalling anaconda3/personal\n\n";
    anaconda-setup
fi


# Create new conda environment called "longi_vl"
if [ -d $HOME/anaconda3/envs/longi_vl ]; then
    echo "###############################################"
    echo -e "\nlongi_vl conda environment is already present"
    echo -e "\nIf you wish you re-install please remove the conda environment first with:"
    echo -e "\tconda remove -n longi_vl --all -y"
    echo -e "\n\n###############################################"
    exit 1
else
    echo -e "\nCreating Conda environment for R packages: longi_vl"
    conda create -n longi_vl -y
    source activate longi_vl
    # export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
fi

# Install initial dependencies
echo -e "\nInstalling Dependencies: R packages via conda"

conda install r r-base r-cmdstanr r-data.table r-loo r-tidybayes r-bayestestR r-bayesplot r-dplyr r-lubridate r-haven r-Hmisc r-DiagrammeR r-DiagrammeRsvg r-rsvg r-htmltools r-raster r-rnaturalearth r-osmdata r-sf r-rstan r-rstanarm r-foreach r-patchwork r-knitr r-foreign r-nnet r-stargazer r-tint r-knitr r-scales r-ggthemes r-gganimate r-ggtext r-ggpubr r-patchwork r-readxl r-rgdal r-rgeos r-RColorBrewer r-data.table r-mvtnorm r-lme4 r-optparse r-bh r-mcmcpack r-here

# also need to install wesanderson in R
