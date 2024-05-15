#!/bin/bash
#SBATCH --job-name=subumap
#SBATCH --partition=gpu
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=15G
#SBATCH --time=7-00:00:00
#SBATCH -o out_slurm/subumap-%j.out
#SBATCH -e out_slurm/subumap-%j.err

#SBATCH --gres=gpu



start=$(date +%s)
echo "starting at $(date) on $(hostname)"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID"

# Run the application
python scVI_umap.py

end=$(date +%s)
echo "ended at $(date) on $(hostname). Time elapsed: $(date -u -d @$((end-start)) +'%H:%M:%S')"
exit 0
