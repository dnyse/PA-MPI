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
#SBATCH --mail-user=dennys.huber@fau.de # Email address for notification

unset SLURM_EXPORT_ENV

module load likwid intel intelmpi
export I_MPI_PIN=1
export I_MPI_DEBUG=0

# --- Configuration ---
BASE_DIR="."
CONFIG_FILE="./poisson.par"
DOMAIN_SIZES=(100 400 1000 4000)
RANKS=(1 2 4 8 16 32 64)
OMEGAS=(1.1 1.5 1.8 1.95)
FILENAME="solver_scalability_results.csv"

# --- Executable Names ---
EXE_STD="./exe-STD"
EXE_RB="./exe-RB"

# --- Functions ---

_run_test() {
    local ranks=$1
    local solver_type=$2 # 'STD' or 'RB'
    local omega=$3
    local N=$4
    local tag=$5
    local executable_path # Determined by solver_type

    # 1. Select the correct executable path
    if [ "${solver_type}" == "RB" ]; then
        executable_path="${EXE_RB}"
    else
        executable_path="${EXE_STD}"
    fi

    # 2. Update poisson.par with the new omega, imax, and jmax values
    sed -i "s/^omg[[:space:]]\+.*$/omg        ${omega}/" "${CONFIG_FILE}"
    sed -i "s/^imax[[:space:]]\+.*$/imax       ${N}/" "${CONFIG_FILE}"
    sed -i "s/^jmax[[:space:]]\+.*$/jmax       ${N}/" "${CONFIG_FILE}"

    # 3. Set MPI process pinning
    np_1=$(($ranks - 1))
    export I_MPI_PIN_PROCESSOR_LIST=0-$np_1

    echo "Running ${solver_type} Solver: Ranks=${ranks}, N=${N}x${N}, Omega=${omega}"

    # 4. Run the selected executable
    result=$(mpirun -n ${ranks} "${executable_path}" "${CONFIG_FILE}" 2>&1 | \
             grep -E "Walltime|Solver took")

    # 5. Extract and Write Results (no change to this logic)
    WALLTIME=$(echo "${result}" | grep "Walltime" | awk '{print $2}' | sed 's/s//')
    NITER=$(echo "${result}" | grep "Solver took" | awk '{print $3}')
    CONVERGENCE=$(echo "${result}" | grep "Solver took" | awk '{print $5}')

    if [ -z "$NITER" ]; then
        NITER="FAIL"
        WALLTIME="FAIL"
        CONVERGENCE="FAIL"
    fi

    echo "${ranks},${N},${omega},${solver_type},${NITER},${CONVERGENCE},${WALLTIME},${tag}" >> $FILENAME
}

# --- Main Execution ---

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
N_STABILITY=1000
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
