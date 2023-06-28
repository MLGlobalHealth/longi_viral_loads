#!/bin/bash

# use as:
# ./run_stan.bash
# Flags:
#    VL : viremic viral load threshold
#    FTP: subset to first time participants or not
#    STANWAY : either (r)stan or cmdstan
#    INDIR : directory where github directory is located
#    OUTDIR : directory where output will be saved

for ARGUMENT in "$@"
do
   KEY=$(echo $ARGUMENT | cut -f1 -d=)

   KEY_LENGTH=${#KEY}
   VALUE="${ARGUMENT:$KEY_LENGTH+1}"

   echo "export $KEY=$VALUE"
   export "$KEY"="$VALUE"
done

# default options
VL="${VL:-1000}" 
FTP="${FTP:-FALSE}"
STANWAY="${STAN_MODEL:-stan}"
INDIR="${INDIR:-/rds/general/user/ab1820/home/git/longi_viral_loads}"
OUTDIR="${OUTDIR:-/rds/general/user/ab1820/home/projects/2022/longvl}"
# ENVNAME="${ENVNAME:-phylowSI-RakaiAgeGender}"

#################
# OPTION CHECKS #
#################

# stop if STANWAY is not stan or cmdstan
if [[ ! "$STANWAY" =~ ^(stan|cmdstan)$ ]]; then
    echo "STANWAY must be either stan (rstan) or cmdstan (cmdstanr)."
    exit 1
fi
# stop if FTP is not True, TRUE, T, False, FALSE, F
if [[ ! "$FTP" =~ ^(TRUE|FALSE)$ ]]; then
    echo "FTP must be either TRUE or FALSE."
    exit 1
fi

echo "Selected options:"
echo "VL = $VL"
echo "FTP = $FTP"
echo "STANWAY = $STANWAY"
echo "INDIR = $INDIR"
echo "OUTDIR = $OUTDIR"

# "build" jobname
JOBNAME="${STANWAY}_vl_$VL"
if [[ "$FTP" == "TRUE" ]]; then
    JOBNAME="$JOBNAME"_firstpart
fi
echo "JOBNAME = $JOBNAME"

########
# MAIN #
########

mkdir $OUTDIR/$JOBNAME

MODELS="run-gp-prevl run-gp-supp-hiv run-gp-supp-pop"
for MODEL in $MODELS
do
    cat > $OUTDIR/bash-$JOBNAME-$MODEL.pbs <<EOF

#!/bin/sh
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=9:ompthreads=1:mem=124gb
#PBS -j oe
#PBS -J 16-19
module load anaconda3/personal
source activate longivl

JOB_TEMP=\${EPHEMERAL}/\${PBS_JOBID}
mkdir -p \$JOB_TEMP
cd \$JOB_TEMP 
PWD=\$(pwd)

INDIR=$INDIR
OUTDIR=$OUTDIR
JOBNAME=$JOBNAME
ROUND=\$PBS_ARRAY_INDEX

# main directory
CWD=\$PWD/\$JOBNAME/$MODEL

mkdir -p \$CWD

Rscript \$INDIR/scripts/VL_run_$STANWAY.R \
    --viremic-viral-load $VL \
    --outdir \$CWD \
    --$MODEL TRUE \
    --round \$ROUND \
    --firstpart $FTP

cp -R --no-preserve=mode,ownership \$PWD/\$JOBNAME/. \$OUTDIR/\$JOBNAME

# cd \$OUTDIR
# qsub bash-$JOBNAME-$MODEL-postprocessing.pbs

EOF

done

cd $OUTDIR

for MODEL in $MODELS
do
    qsub bash-${JOBNAME}-${MODEL}.pbs
done
