#!/bin/sh

VL="${1:=1000}" 
echo "VL set to @VL"
# VL=1000
JOBNAME="vl_$VL"
INDIR="/rds/general/user/ab1820/home/git/longi_viral_loads"
OUTDIR="/rds/general/user/ab1820/home/projects/2022/longvl"


cat > $OUTDIR/bash-$JOBNAME-postprocessing.pbs <<EOF
  
#!/bin/sh
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=10:ompthreads=1:mem=480gb
#PBS -j oe
module load anaconda3/personal

INDIR=$INDIR
OUTDIR=$OUTDIR
STAN_MODEL=$STAN_MODEL
JOBNAME=$JOBNAME
  
# main directory
CWD=\$PWD/\$JOBNAME

# directories for figure and table
mkdir \$CWD
mkdir \$CWD/figures
mkdir \$CWD/tables

Rscript \$INDIR/scripts/VL_scripts/VL_postprocessing.R -indir --viremic-viral-load $VL --outdir \$CWD --indir \$CWD

EOF

cd $OUTDIR
qsub bash-$JOBNAME-postprocessing.pbs
