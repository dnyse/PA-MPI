#!/bin/bash -l
#SBATCH --job-name=bench_internode
#SBATCH --output=Fritz_ICX_3D_Solver_internode_%j.out
#SBATCH --error=Fritz_ICX_3D_Solver_internode_%j.err
#SBATCH --partition=multinode
#SBATCH --nodes=4
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

FILENAME="result_internode_scaling.csv"

# Domain configuration
DOMAIN="600x200x200"
OMEGA="1.8"

cd /home/hpc/pavl/pavl164v/exercises/PAMPI/assignment7/benchmarking

rm -f $FILENAME
touch $FILENAME
echo "Ranks,Nodes,Time_seconds,Domain,Omega" >> $FILENAME

# Scale across nodes: 72, 144, 216, 288 cores (1, 2, 3, 4 nodes)
for nodes in 1 2 3 4; do
    np=$(($nodes * 72))
    
    echo "Running with $np ranks on $nodes node(s)..."
    
    # Use likwid-mpirun with socket-level pinning across all nodes
    # -npn specifies number of processes per node
    result=$(likwid-mpirun -mpi slurm -np $np -nperdomain N:72 ./exe-ICX canal_cluster.par 2>&1)
    
    time=$(echo "$result" | grep "Solution took" | grep -oP "\d+\.\d+(?=s)")
    
    if [ -z "$time" ]; then
        echo "WARNING: Could not extract time for $np ranks"
        time="NA"
    fi
    
    echo "$np,$nodes,$time,$DOMAIN,$OMEGA" >> $FILENAME
    echo "Completed: $np ranks ($nodes nodes) - Time: $time seconds"
done

echo "Internode scaling benchmark completed. Results in $FILENAME"
echo "Domain: $DOMAIN, Omega: $OMEGA"
