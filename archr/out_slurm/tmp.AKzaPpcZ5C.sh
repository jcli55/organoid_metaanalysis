#!/bin/bash
#SBATCH --job-name=tmp.AKzaPpcZ5C
#SBATCH --partition=interactive
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --time=3-00:00:00
#SBATCH -o out_slurm/tmp.AKzaPpcZ5C-%j.out
#SBATCH -e out_slurm/tmp.AKzaPpcZ5C-%j.err





start=$(date +%s)
echo "starting at $(date) on $(hostname)"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID"

# Run the application
Rscript temp.R

end=$(date +%s)
echo "ended at $(date) on $(hostname). Time elapsed: $(date -u -d @$((end-start)) +'%H:%M:%S')"
exit 0
