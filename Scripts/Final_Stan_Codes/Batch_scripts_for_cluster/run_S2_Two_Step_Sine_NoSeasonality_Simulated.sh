#!/bin/bash

#SBATCH --job-name=cabueico_job
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -p long
#SBATCH --requeue

# Specify the working directory
cd /exafs1/well/cooper-who/users/bdi478/

# Load required modules (adjust as necessary for your environment)
module load R/4.4.1

# Pass variable to R script
export scenario="Two_step_sine_noseasonality" 
export data_source="simulated"  # can be "simulated" or "observed"


# -------------------------------
# Check variable
# -------------------------------
echo "Bash sees scenario: $scenario"
echo "Bash sees data source: $data_source"


# Run R script 
Rscript ./CABU_EICO/model-code/S2_Two_Step_Sine.R "$scenario" "$data_source"
