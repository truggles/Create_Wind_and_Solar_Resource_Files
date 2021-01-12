#!/bin/bash
#
#SBATCH --export=ALL
#SBATCH -n 1
#SBATCH -p DGE
#SBATCH --mem=80000
#
#cd $SLURM_SUBMIT_DIR
#
#module load anaconda/anaconda3
module load python/3.6.0
#source activate cdat8
#conda init bash
conda deactivate
conda activate geo_stuff

export UVCDAT_ANONYMOUS_LOG=no
echo "hello world, ${YEAR}"
python step0_get_temperature_profile.py $YEAR
