#!/bin/bash -l
#SBATCH --job-name=bench_intra_numa
#SBATCH --output=Fritz_ICX_3D_Solver_intra_numa_%j.out
#SBATCH --error=Fritz_ICX_3D_Solver_intra_numa_%j.err
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=18
#SBATCH --time=10:00:00
#SBATCH --export=NONE
#SBATCH --cpu-freq=2400000-2400000:performance
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

unset SLURM_EXPORT_ENV

module load likwid intel intelmpi

# LIKWID settings
export LIKWID_PIN=1
export LIKWID_SILENT=1

FILENAME="result_intra_numa_scaling.csv"

# Domain configuration
DOMAIN="200x40x40"
OMEGA="1.8"

cd /home/hpc/pavl/pavl164v/exercises/PAMPI/assignment7/benchmarking

rm -f $FILENAME
touch $FILENAME
echo "Ranks,Time_seconds,Domain,Omega" >> $FILENAME

# Scale from 1 to 18 cores within one NUMA domain
for np in $(seq 1 18); do
    echo "Running with $np ranks..."
    
    # Use likwid-mpirun with pinning to NUMA domain 0 (cores 0-17)
    np_1=$(($np - 1))
    result=$(likwid-mpirun -mpi slurm -np $np -pin M0:0-$np_1 ./exe-ICX canal_memdomain.par 2>&1)
    
    time=$(echo "$result" | grep "Solution took" | grep -oP "\d+\.\d+(?=s)")
    
    if [ -z "$time" ]; then
        echo "WARNING: Could not extract time for $np ranks"
        time="NA"
    fi
    
    echo "$np,$time,$DOMAIN,$OMEGA" >> $FILENAME
    echo "Completed: $np ranks - Time: $time seconds"
done

echo "Intra-NUMA scaling benchmark completed. Results in $FILENAME"
echo "Domain: $DOMAIN, Omega: $OMEGA"
