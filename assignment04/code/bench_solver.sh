#!/bin/bash -l
#SBATCH --job-name=solver_scalability_study
#SBATCH --output=solver_study_output
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --time=10:00:00
#SBATCH --export=NONE
#SBATCH --cpu-freq=2400000-2400000:performance
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${EMAIL} # Email address for notification

unset SLURM_EXPORT_ENV

module load likwid intel intelmpi
export I_MPI_PIN=1
export I_MPI_DEBUG=0

# --- Configuration ---
BASE_DIR="."
CONFIG_FILE="./poisson.par"
DOMAIN_SIZES=(100 200 300)
RANKS=(1 2 4 8 16 32 64 72)
OMEGAS=(1.1 1.5 1.8 1.95)
FILENAME="solver_scalability_results.csv"

# --- Executable Names ---
EXE_STD="./exe-STD"
EXE_RB="./exe-RB"

# --- Functions ---

# --- Main Execution ---
_run_test() {
    local ranks=$1
    local solver_type=$2 # Can be 'STD' or 'RB'
    local omega=$3
    local N=$4
    local tag=$5
    local executable_path
    local FULL_OUTPUT
    local result
    local NITER
    local CONVERGENCE
    local WALLTIME

    if [ "${solver_type}" == "RB" ]; then
    	executable_path="${EXE_RB}"
    else
	executable_path="${EXE_STD}"
    fi
    sed -i "s/^omg[[:space:]]\+.*$/omg        ${omega}/" "${CONFIG_FILE}"
    sed -i "s/^imax[[:space:]]\+.*$/imax       ${N}/" "${CONFIG_FILE}"
    sed -i "s/^jmax[[:space:]]\+.*$/jmax       ${N}/" "${CONFIG_FILE}"
    # Increase itermax for large domains to avoid premature failure
    sed -i "s/^itermax[[:space:]]\+.*$/itermax  500000/" "${CONFIG_FILE}"

    # 4. Set MPI process pinning
    np_1=$(($ranks - 1))
    export I_MPI_PIN_PROCESSOR_LIST=0-$np_1

    echo "Running ${solver_type} Solver: Ranks=${ranks}, N=${N}x${N}, Omega=${omega}"

# 1. Run the executable and capture ALL output (stdout + stderr) into the variable
    FULL_OUTPUT=$(mpirun -n ${ranks} "${executable_path}" "${CONFIG_FILE}" 2>&1)
	        
    # NEW: Print the full output to stdout (which is your SLURM output file)
    echo "--- Full Output for Ranks=${ranks}, N=${N} ---"
    echo "${FULL_OUTPUT}"
    echo "-----------------------------------"

			        # 2. Filter the captured output to get the result lines needed for the CSV
    result=$(echo "${FULL_OUTPUT}" | grep -E "Walltime|Solver took")
    # Extract Walltime and Iterations
    # NITER is the 3rd field: "Solver took [28172]..."
    NITER=$(echo "${result}" | grep "Solver took" | awk '{print $3}')
    
    # CONVERGENCE is the 7th field: "...to reach [0.000001]..."
    CONVERGENCE=$(echo "${result}" | grep "Solver took" | awk '{print $7}') # <-- CORRECTED FIELD NUMBER

    # WALLTIME is the 2nd field of the Walltime line, with 's' removed.
    WALLTIME=$(echo "${result}" | grep "Walltime" | awk '{print $2}' | sed 's/s//')

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

# --- STEP 1: Compilation ---

echo "--- COMPILATION: Standard Solver ---"
# Compile Standard (STD) version (no -DRED_BLACK_SOR)
make -C "${BASE_DIR}" DEFINES="-D_GNU_SOURCE"
mv "${BASE_DIR}/exe-ICC" "${EXE_STD}"

echo "--- COMPILATION: Red-Black SOR Solver ---"
# Compile Red-Black (RB) version (-DRED_BLACK_SOR)
make -C "${BASE_DIR}" DEFINES="-DRED_BLACK_SOR -D_GNU_SOURCE"
mv "${BASE_DIR}/exe-ICC" "${EXE_RB}"

# --- STEP 2: Part (a) Scalability Study (Omega=1.1) ---

echo "--- STARTING: Part (a) Scalability Study (Omega=1.1) ---"
for N in "${DOMAIN_SIZES[@]}"; do
    for ranks in "${RANKS[@]}"; do
        # Standard Solver
       _run_test $ranks "STD" 1.1 $N "A_Scalability"

        # Red-Black SOR Solver
        _run_test $ranks "RB" 1.1 $N "A_Scalability"
    done
done

# --- STEP 3: Part (b) Stability Study (N=1000) ---

echo "--- STARTING: Part (b) Stability Study (N=1000) ---"
N_STABILITY=250
for omega in "${OMEGAS[@]}"; do
    for ranks in "${RANKS[@]}"; do
        _run_test $ranks "STD" $omega $N_STABILITY "B_Stability"
        _run_test $ranks "RB" $omega $N_STABILITY "B_Stability"
    done
done

# --- STEP 4: Final Cleanup ---

# Resetting the config file to a known good state
sed -i "s/^omg[[:space:]]\+.*$/omg        1.1/" "${CONFIG_FILE}"
sed -i "s/^imax[[:space:]]\+.*$/imax       100/" "${CONFIG_FILE}"
sed -i "s/^jmax[[:space:]]\+.*$/jmax       100/" "${CONFIG_FILE}"

echo "--- STUDY COMPLETE. Results saved to ${FILENAME} ---"
