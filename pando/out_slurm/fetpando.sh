#!/bin/bash
#SBATCH --job-name=fetpando
#SBATCH --partition=short
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH -o out_slurm/fetpando-%j.out
#SBATCH -e out_slurm/fetpando-%j.err
#SBATCH --nodelist=mhgcp-d02




start=$(date +%s)
echo "starting at $(date) on $(hostname)"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID"

# Run the application
Rscript run_pando_fet.R

end=$(date +%s)
echo "ended at $(date) on $(hostname). Time elapsed: $(date -u -d @$((end-start)) +'%H:%M:%S')"
exit 0
