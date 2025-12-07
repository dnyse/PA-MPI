#!/bin/bash
#SBATCH --job-name=3D-Solver-Serial           # Job name
#SBATCH --output=%x_%j.out          # Output file (%x=job-name, %j=job-id)
#SBATCH --error=%x_%j.err           # Error file
#SBATCH --partition=singlenode           # Partition/queue name
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks=1                  # Number of MPI tasks
#SBATCH --time=02:00:00             # Wall clock time limit (HH:MM:SS)
#SBATCH --mail-type=END,FAIL        # Email notifications
#SBATCH --mail-user=dennys.huber@fau.de  # Your email address

# Load required modules
module purge
module load intel

# Print job information
echo "Job started on $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Working directory: $(pwd)"
echo "Number of tasks: $SLURM_NTASKS"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"

make distclean
make

echo "Run Canal Serial"
srun ./exe-ICC canal.par
echo "Run Dcavity Serial"
srun ./exe-ICC dcavity.par

echo "Job finished on $(date)"
