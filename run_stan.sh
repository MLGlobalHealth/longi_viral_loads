#!/bin/sh

# STAN_MODEL=TODO
VL=1000
JOBNAME="vl_$VL"
INDIR="/rds/general/user/ab1820/home/git/longi_viral_loads"
OUTDIR="/rds/general/user/ab1820/home/projects/2022/longvl"

mkdir $OUTDIR/$JOBNAME

for MODEL in run-gp-prevl run-gp-supp-hiv run-gp-supp-pop run-icar-mean-vl 
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
OUTDIR=$OUTDIR/$JOBNAME
JOBNAME=$JOBNAME
ROUND=\$PBS_ARRAY_INDEX

# main directory
CWD=\$PWD/\$JOBNAME/$MODEL

mkdir -p \$CWD

Rscript \$INDIR/scripts/VL_run_stan.R --viremic-viral-load $VL --outdir \$CWD --$MODEL TRUE --round \$ROUND

cp -R --no-preserve=mode,ownership \$PWD/\$JOBNAME/. \$OUTDIR

# cd \$OUTDIR
# qsub bash-$JOBNAME-$MODEL-postprocessing.pbs

EOF

done

qsub bash-$JOBNAME-run-gp-supp-hiv.pbs
qsub bash-$JOBNAME-run-gp-supp-pop.pbs
qsub bash-$JOBNAME-run-icar-mean-vl.pbs
