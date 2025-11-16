#!/bin/bash -l
#SBATCH --job-name=solver_scalability_study
#SBATCH --output=solver_study_output
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --time=10:00:00
#SBATCH --export=NONE
#SBATCH --cpu-freq=2400000-2400000:performance

unset SLURM_EXPORT_ENV

module load likwid intel intelmpi
# Pinning processes and fixing frequency as required by (a)
export I_MPI_PIN=1
export I_MPI_DEBUG=0

# --- Configuration ---
# Number of processes per node (assuming Intel Ice Lake / 72 logical cores)
NPM=18
# Base directory for the project
BASE_DIR="~/PA-MPI/assignment04/code/"
# Configuration file for domain size and omega
CONFIG_FILE="${BASE_DIR}/poisson.par"

# Domain sizes (imax and jmax, where imax=jmax=N)
DOMAIN_SIZES=(100 400 1000 4000)
# Number of MPI ranks to test (1, 2, 4, 8, 16, 32, 64 ranks on one node)
RANKS=(1 2 4 8 16 32 64)
# Relaxation parameters (omega) to test for stability (original is 1.1 )
OMEGAS=(1.1 1.5 1.8 1.95)
# Name of the output CSV file
FILENAME="solver_scalability_results.csv"

# --- Functions ---

_run_test() {
    local ranks=$1
    local solver_type=$2 # Can be 'STD' or 'RB'
    local omega=$3
    local N=$4
    local tag=$5

    # 1. Update poisson.par with the new omega value
   
    # Replace omg, imax, and jmax values in the config file
    sed -i "s/^omg[[:space:]]\+.*$/omg        ${omega}/" "${CONFIG_FILE}"
    sed -i "s/^imax[[:space:]]\+.*$/imax       ${N}/" "${CONFIG_FILE}"
    sed -i "s/^jmax[[:space:]]\+.*$/jmax       ${N}/" "${CONFIG_FILE}"

    local DEFS_FLAG=""
    if [ "${solver_type}" == "RB" ]; then
        DEFS_FLAG="-DRED_BLACK_SOR"
    fi

    # 3. Recompile with the required DEFINES (if any)
    # Using 'make clean' is safer but slower, 'make' usually recompiles if flags change.
    # We pass the DEFINES flag to the Makefile.
    make -C "${BASE_DIR}" DEFINES="${DEFS_FLAG} -D_GNU_SOURCE"

    # 4. Set MPI process pinning
    np_1=$(($ranks - 1))
    export I_MPI_PIN_PROCESSOR_LIST=0-$np_1

    echo "Running ${solver_type} Solver: Ranks=${ranks}, N=${N}x${N}, Omega=${omega}"

    result=$(mpirun -n ${ranks} "${BASE_DIR}/exe-ICX" "${CONFIG_FILE}" 2>&1 | \
             grep -E "Walltime|Solver took")

    # Extract Walltime and Iterations
    WALLTIME=$(echo "${result}" | grep "Walltime" | awk '{print $2}' | sed 's/s//')
    NITER=$(echo "${result}" | grep "Solver took" | awk '{print $3}')

    # Extract the convergence residuum (sqrt(res)) at termination
    CONVERGENCE=$(echo "${result}" | grep "Solver took" | awk '{print $5}')
    
    # Check if solver completed (NITER found) or failed (NITER is empty)
    if [ -z "$NITER" ]; then
        NITER="FAIL"
        WALLTIME="FAIL"
        CONVERGENCE="FAIL"
    fi

    # 5. Write results to CSV
    echo "${ranks},${N},${omega},${solver_type},${NITER},${CONVERGENCE},${WALLTIME},${tag}" >> $FILENAME
}

cd "${BASE_DIR}"
make distclean

rm $FILENAME
touch $FILENAME
echo "Ranks,DomainSize,Omega,SolverType,Iterations,Convergence,Walltime,TestTag" >>$FILENAME

# --- Part (a) Scalability Study on one node
# We will test both the standard and Red-Black versions.

echo "--- STARTING: Part (a) Scalability Study (Omega=1.1) ---"
for N in "${DOMAIN_SIZES[@]}"; do
    for ranks in "${RANKS[@]}"; do
        # Standard Solver (default)
        _run_test $ranks "STD" 1.1 $N "A_Scalability"
        
        # Red-Black SOR Solver (requires -DRED_BLACK_SOR flag)
        _run_test $ranks "RB" 1.1 $N "A_Scalability"
    done
done

echo "--- STARTING: Part (b) Stability Study (N=1000) ---"
N_STABILITY=1000
for omega in "${OMEGAS[@]}"; do
    for ranks in "${RANKS[@]}"; do
        _run_test $ranks "STD" $omega $N_STABILITY "B_Stability"
        _run_test $ranks "RB" $omega $N_STABILITY "B_Stability"
    done
done

# Resetting the config file to a known good state
sed -i "s/^omg[[:space:]]\+.*$/omg        1.1/" "${CONFIG_FILE}"
sed -i "s/^imax[[:space:]]\+.*$/imax       100/" "${CONFIG_FILE}"
sed -i "s/^jmax[[:space:]]\+.*$/jmax       100/" "${CONFIG_FILE}"

echo "--- STUDY COMPLETE. Results saved to ${FILENAME} ---"
