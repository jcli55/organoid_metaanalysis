#!/bin/bash
#SBATCH --job-name=nrpcmg
#SBATCH --partition=gpu
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --time=7-00:00:00
#SBATCH -o out_slurm/nrpcmg-%j.out
#SBATCH -e out_slurm/nrpcmg-%j.err

#SBATCH --gres=gpu



start=$(date +%s)
echo "starting at $(date) on $(hostname)"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID"

# Run the application
python annotate_mg_NRPC.py

end=$(date +%s)
echo "ended at $(date) on $(hostname). Time elapsed: $(date -u -d @$((end-start)) +'%H:%M:%S')"
exit 0
