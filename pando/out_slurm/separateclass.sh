#!/bin/bash
#SBATCH --job-name=separateclass
#SBATCH --partition=short
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=25G
#SBATCH --time=12:00:00
#SBATCH -o out_slurm/separateclass-%j.out
#SBATCH -e out_slurm/separateclass-%j.err
#SBATCH --nodelist=mhgcp-d01




start=$(date +%s)
echo "starting at $(date) on $(hostname)"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID"

# Run the application
Rscript separate_class.R

end=$(date +%s)
echo "ended at $(date) on $(hostname). Time elapsed: $(date -u -d @$((end-start)) +'%H:%M:%S')"
exit 0
