#!/bin/bash
#SBATCH --job-name=tmp.qFz0SaE9VL
#SBATCH --partition=short
#SBATCH --ntasks-per-node=8
#SBATCH --nodes=1
#SBATCH --mem=25G
#SBATCH --time=12:00:00
#SBATCH -o out_slurm/tmp.qFz0SaE9VL-%j.out
#SBATCH -e out_slurm/tmp.qFz0SaE9VL-%j.err
#SBATCH --nodelist=mhgcp-d03




start=$(date +%s)
echo "starting at $(date) on $(hostname)"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID"

# Run the application
Rscript create_archr_project.R

end=$(date +%s)
echo "ended at $(date) on $(hostname). Time elapsed: $(date -u -d @$((end-start)) +'%H:%M:%S')"
exit 0
