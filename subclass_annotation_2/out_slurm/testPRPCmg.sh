#!/bin/bash
#SBATCH --job-name=testPRPCmg
#SBATCH --partition=gpu
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=12G
#SBATCH --time=7-00:00:00
#SBATCH -o out_slurm/testPRPCmg-%j.out
#SBATCH -e out_slurm/testPRPCmg-%j.err

#SBATCH --gres=gpu



start=$(date +%s)
echo "starting at $(date) on $(hostname)"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID"

# Run the application
python temp.py

end=$(date +%s)
echo "ended at $(date) on $(hostname). Time elapsed: $(date -u -d @$((end-start)) +'%H:%M:%S')"
exit 0
