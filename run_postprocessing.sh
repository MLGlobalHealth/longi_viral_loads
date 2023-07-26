#!/bin/bash
ABSOLUTE_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
ABSOLUTE_PATH="$(dirname $ABSOLUTE_PATH)"

if [[ "$1" == "-h" || "$1" == "-help" || "$1" == "--help" ]]; then
    echo "Usage: ./run_postprocessing OPTION1=value1 OPTION2=value2 ..."
    echo "OPTIONs:"
    echo "    VL : viremic viral load threshold             [default: 1000]"
    echo "    STANWAY : either (r)stan or cmdstan           [default: stan]"
    echo "    INDIR : directory where github directory is located"
    echo "    OUTDIR : directory where output will be saved"
    exit 0
fi

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

# "build" jobname and envname
ENVNAME="longivl"
JOBNAME="vl_$VL"
if [[ "$STANWAY" == "cmdstan" ]]; then
	JOBNAME="${STANWAY}_${JOBNAME}"
	ENVNAME="${ENVNAME}_cmdstan"
fi
echo "JOBNAME = $JOBNAME"


cat > $OUTDIR/bash-$JOBNAME-joint_postprocessing.pbs <<EOF
  
#!/bin/sh
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=10:ompthreads=1:mem=480gb
#PBS -j oe
module load anaconda3/personal
source activate longivl

INDIR=$INDIR
OUTDIR=$OUTDIR
STAN_MODEL=$STAN_MODEL
JOBNAME=$JOBNAME
  
Rscript \$INDIR/scripts/VL_jointpostprocessing.R \\
    --viremic-viral-load $VL \\
    --jobname \$JOBNAME \\
    --outdir \$OUTDIR \\
    --indir \$OUTDIR

Rscript \$INDIR/src/make_paper_figures.R \\
    --viremic-viral-load $VL \\
    --jobname \$JOBNAME \\
    --outdir \$OUTDIR \\
    --indir \$OUTDIR

EOF

cd $OUTDIR
qsub bash-$JOBNAME-postprocessing.pbs
