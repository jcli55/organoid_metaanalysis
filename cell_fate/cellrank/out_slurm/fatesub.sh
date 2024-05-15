#!/bin/bash
#SBATCH --job-name=fatesub
#SBATCH --partition=short
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --time=12:00:00
#SBATCH -o out_slurm/fatesub-%j.out
#SBATCH -e out_slurm/fatesub-%j.err





start=$(date +%s)
echo "starting at $(date) on $(hostname)"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID"

# Run the application
python infer_fate_subclass.py

end=$(date +%s)
echo "ended at $(date) on $(hostname). Time elapsed: $(date -u -d @$((end-start)) +'%H:%M:%S')"
exit 0
