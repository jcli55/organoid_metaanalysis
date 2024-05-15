#!/bin/bash
#SBATCH --job-name=geneint
#SBATCH --partition=interactive
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --time=3-00:00:00
#SBATCH -o out_slurm/geneint-%j.out
#SBATCH -e out_slurm/geneint-%j.err





start=$(date +%s)
echo "starting at $(date) on $(hostname)"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID"

# Run the application
Rscript archr_gene_integration.R

end=$(date +%s)
echo "ended at $(date) on $(hostname). Time elapsed: $(date -u -d @$((end-start)) +'%H:%M:%S')"
exit 0
