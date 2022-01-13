#!/bin/bash
#######################################################################################
#SBATCH --time=72:00:00
#SBATCH --job-name=prospect-slurm
#SBATCH --out="prospect_job-%j.out"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --tasks-per-node=2
#SBATCH --mem-per-cpu=20GB
#SBATCH --mail-type=ALL
#######################################################################################

#######################################################################################
### SLURM usage example
# sbatch -w node03 -c 3 --partition compute --job-name=prospect --mail-user=sserbin@bnl.gov batch/invert_rrtm_prospect_slurm.sh
# sbatch -w node03 --partition compute --job-name=prospect --mail-user=sserbin@bnl.gov batch/invert_rrtm_prospect_slurm.sh
#######################################################################################

# which version of the script/model? 
# 5b
script_name=Rscripts/PROSPECT/Invert_leaf_refl_spectra_rrtm_PROSPECT-5b.R
# D

echo "Running ${script_name}"
echo "Starting at: `date`"
echo "Job submitted to the ${SLURM_JOB_PARTITION} partition on ${SLURM_CLUSTER_NAME}"
echo "Job name: ${SLURM_JOB_NAME}, Job ID: ${SLURM_JOB_ID}"
echo "There are ${SLURM_CPUS_ON_NODE} CPUs on compute node $(hostname)"

echo "Run started on: `date`"

Rscript --vanilla ${script_name}
