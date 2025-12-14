#!/bin/bash -l
#SBATCH --job-name=bench_inter_numa
#SBATCH --output=Fritz_ICX_3D_Solver_inter_numa_%j.out
#SBATCH --error=Fritz_ICX_3D_Solver_inter_numa_%j.err
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
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

FILENAME="result_inter_numa_scaling.csv"

# Domain configuration
DOMAIN="300x150x150"
OMEGA="1.8"

cd /home/hpc/pavl/pavl164v/exercises/PAMPI/assignment7/benchmarking

rm -f $FILENAME
touch $FILENAME
echo "Ranks,Time_seconds,Domain,Omega" >> $FILENAME

# Scale across NUMA domains: 18, 36, 54, 72 cores
# Assuming 4 NUMA domains with 18 cores each
for np in 18 36 54 72; do
    echo "Running with $np ranks..."
    
    # np=18, nperdomain M:18
    # -> Uses 1 NUMA domain with 18 processes

    # np=36, nperdomain M:18  
    # -> Uses 2 NUMA domains with 18 each

    # np=54, nperdomain M:18
    # -> Uses 3 NUMA domains with 18 each

    # np=72, nperdomain M:18
    # -> Uses 4 NUMA domains with 18 each
    result=$(likwid-mpirun -mpi slurm -np $np -nperdomain M:18 ./exe-ICX canal_node.par 2>&1)
    
    time=$(echo "$result" | grep "Solution took" | grep -oP "\d+\.\d+(?=s)")
    
    if [ -z "$time" ]; then
        echo "WARNING: Could not extract time for $np ranks"
        time="NA"
    fi
    
    echo "$np,$time,$DOMAIN,$OMEGA" >> $FILENAME
    echo "Completed: $np ranks - Time: $time seconds"
done

echo "Inter-NUMA scaling benchmark completed. Results in $FILENAME"
echo "Domain: $DOMAIN, Omega: $OMEGA"
