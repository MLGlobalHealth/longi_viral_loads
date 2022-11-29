#!/bin/sh

# STAN_MODEL=TODO
JOBNAME="viremic200"
INDIR="/rds/general/user/ab1820/home/git/longi_viral_loads"
OUTDIR="/rds/general/user/ab1820/home/projects/2022/longvl"

mkdir $OUTDIR

for MODEL in run-gp-prevl run-gp-supp-hiv run-gp-supp-pop run-icar-mean-vl 
do
    cat > $OUTDIR/bash-$JOBNAME-$MODEL.pbs <<EOF

#!/bin/sh
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=10:ompthreads=1:mem=240gb
#PBS -j oe
module load anaconda3/personal
source activate longivl

JOB_TEMP=\${EPHEMERAL}/\${PBS_JOBID}
mkdir -p \$JOB_TEMP
cd \$JOB_TEMP 
PWD=\$(pwd)

INDIR=$INDIR
OUTDIR=$OUTDIR
# STAN_MODEL=$STAN_MODEL
JOBNAME=$JOBNAME

# main directory
CWD=\$PWD/\$JOBNAME-\$MODEL

mkdir \$CWD

Rsript \$INDIR/scripts/analyse_all_participants.R --viremic-viral-load 200 --outdir \$CWD --$MODEL TRUE 

cp -R --no-preserve=mode,ownership \$PWD/* \$OUTDIR

# cd \$OUTDIR
# qsub bash-$JOBNAME-$MODEL-postprocessing.pbs

EOF

done



# TODO: postprocessing
# cat > $OUTDIR/bash_$STAN_MODEL-$JOBNAME-postprocessing.pbs <<EOF
#   
# #!/bin/sh
# #PBS -l walltime=24:00:00
# #PBS -l select=1:ncpus=10:ompthreads=1:mem=480gb
# #PBS -j oe
# module load anaconda3/personal
# 
# INDIR=$INDIR
# OUTDIR=$OUTDIR
# STAN_MODEL=$STAN_MODEL
# JOBNAME=$JOBNAME
#   
# # main directory
# CWD=\$OUTDIR/\$STAN_MODEL-\$JOBNAME
# 
# # directories for figure and table
# mkdir \$CWD/figures
# mkdir \$CWD/tables
# 
# Rscript \$INDIR/scripts/postprocessing_assess_mixing.R -indir \$INDIR -outdir \$CWD -stan_model \$STAN_MODEL -jobname \$JOBNAME 
# Rscript \$INDIR/scripts/postprocessing_figures.R -indir \$INDIR -outdir \$CWD -stan_model \$STAN_MODEL -jobname \$JOBNAME 
# 
# EOF
  
cd $OUTDIR
qsub bash-$JOBNAME.pbs

