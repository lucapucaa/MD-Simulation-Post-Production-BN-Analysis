#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=dyrk1a_array
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=2-00:00:00
#SBATCH --array=1-5                       # 🚀 Runs all 5 replicas in parallel!
#SBATCH --output=logs/slurm_%A_%a.log     # 🚀 Individual log files (%A=job, %a=array_id)

set -e

# ==============================================================================
# 1. SETUP CENTRALIZED WORKSPACE & ARRAY VARIABLES
# ==============================================================================
PROJECT_ROOT="/home/lchill/work/pproduct_analysis"
cd "$PROJECT_ROOT/pipeline"

REP_NUM="${SLURM_ARRAY_TASK_ID}"
RUN_ID="L207I_rep_${REP_NUM}"

# 🎯 DYNAMIC PATH MAPPING (Reads directly from main input folder)
RAW_PDB="7o7k_L207I.pdb"
RAW_TPR="L207I_mut_sim/L207I_rep_${REP_NUM}.tpr"
RAW_XTC="L207I_mut_sim/final_delivery/replica_${REP_NUM}/merged_raw.xtc"
LIGAND="ATP"

# 💡 UNIQUE FILE NAMES: Prevents parallel replicas from overwriting each other
CLEAN_XTC="${RUN_ID}_final.xtc"
MATCHED_PDB="matched_topology_${RUN_ID}.pdb"

# Check that the specific replica files actually exist before starting
if [ ! -f "${PROJECT_ROOT}/input/${RAW_TPR}" ] || [ ! -f "${PROJECT_ROOT}/input/${RAW_XTC}" ]; then
    echo "❌ Error: Missing input TPR or XTC for replica ${REP_NUM}."
    exit 1
fi

# ==============================================================================
# 2. ISOLATED STATUS LOG SETUP
# ==============================================================================
LOGS_DIR="${PROJECT_ROOT}/logs"
mkdir -p "$LOGS_DIR"

STATUS_LOG="${LOGS_DIR}/${RUN_ID}_pipeline_status.log"
rm -f "$STATUS_LOG" # Start fresh for this specific array task

echo "==================================================" > "$STATUS_LOG"
echo "    PARALLEL PIPELINE LEDGER FOR: ${RUN_ID}" >> "$STATUS_LOG"
echo "==================================================" >> "$STATUS_LOG"

log_phase() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ✅ $1" >> "$STATUS_LOG"
}

failure_handler() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ❌ FATAL ERROR: ${RUN_ID} crashed on or around line $1" >> "$STATUS_LOG"
}
trap 'failure_handler $LINENO' ERR

# ==============================================================================
# PIPELINE EXECUTION ENGINE (CONCURRENT RUNS)
# ==============================================================================
log_phase "START: Initializing array task ${REP_NUM}."

echo "=== Loading Global Cluster Python & GROMACS Environments ==="
module load Python/3.9.5-GCCcore-10.3.0 2>/dev/null || true
module load Gromacs/2023.2-Container 2>/dev/null || true
source ~/dyrk1a_venv/bin/activate

# --- PHASE 1: Trajectory Boundary Cleaning ---
echo "Starting Phase 1: Cleaning Trajectory..."
./01_clean_traj.sh "$RAW_TPR" "$RAW_XTC" "$LIGAND" "$RUN_ID"
log_phase "PHASE_1_COMPLETE: GROMACS trajectory boundary unwrapping finished."

# --- PHASE 1.5: GROMACS Topology Alignment ---
echo "Running Phase 1.5: Generating perfectly matched PDB from TPR..."
if ! command -v gmx &> /dev/null; then
    if [ -f "/home/lchill/miniconda3/envs/dyrk1a/bin.AVX2_256/gmx" ]; then
        export PATH="/home/lchill/miniconda3/envs/dyrk1a/bin.AVX2_256:$PATH"
    elif [ -f "/home/lchill/miniconda3/envs/dyrk1a/bin/gmx" ]; then
        export PATH="/home/lchill/miniconda3/envs/dyrk1a/bin:$PATH"
    fi
fi
echo -e "1 | r ${LIGAND}\nname 17 Protein_Ligand\nq" | gmx make_ndx -f "${PROJECT_ROOT}/input/${RAW_TPR}" -o "${LOGS_DIR}/${RUN_ID}_temp_match.ndx"
echo "Protein_Ligand" | gmx editconf -f "${PROJECT_ROOT}/input/${RAW_TPR}" -n "${LOGS_DIR}/${RUN_ID}_temp_match.ndx" -o "${PROJECT_ROOT}/input/${MATCHED_PDB}"
rm -f "${LOGS_DIR}/${RUN_ID}_temp_match.ndx"
log_phase "PHASE_1_5_COMPLETE: Replica-specific topology generated via GROMACS editconf."

# --- PHASE 2: Quality Validation ---
echo "Running Phase 2: Quality Validation..."
./02_quality_metrics.sh "$MATCHED_PDB" "$CLEAN_XTC" "$LIGAND" "$RUN_ID"
log_phase "PHASE_2_COMPLETE: Structural quality evaluations calculated."

# --- PARALLEL BLOCK (PHASE 3 & PHASE 4) ---
echo "🚀 Spawning parallel processing block..."
PID_PHASE3=""
PID_PHASE4=""

echo "-> Phase 3 (GetContacts Fingerprinting) running in background..."
./03_contacts.sh "$MATCHED_PDB" "$CLEAN_XTC" "$RUN_ID" "$LIGAND" > "${LOGS_DIR}/${RUN_ID}_phase3_stream.log" 2>&1 &
PID_PHASE3=$!

echo "-> Phase 4 (Per-Residue Energy Mapping) running in background..."
python -u 04_perres_energy.py \
    --run_id "$RUN_ID" \
    --pdb "$MATCHED_PDB" \
    --xtc "$CLEAN_XTC" \
    --ligand "$LIGAND" \
    --cutoff 12.0 \
    --cpus 16 > "${LOGS_DIR}/${RUN_ID}_phase4_stream.log" 2>&1 &
PID_PHASE4=$!

# Synchronize background tasks
echo "⏳ Syncing Phase 3 and Phase 4 threads..."
wait $PID_PHASE3 || { echo "❌ Phase 3 async task threw an error."; exit 1; }
log_phase "PHASE_3_COMPLETE: GetContacts geometric fingerprinting finished."
rm -f "${LOGS_DIR}/${RUN_ID}_phase3_stream.log"

wait $PID_PHASE4 || { echo "❌ Phase 4 async task threw an error."; exit 1; }
log_phase "PHASE_4_COMPLETE: Continuous interaction energy extraction complete."
rm -f "${LOGS_DIR}/${RUN_ID}_phase4_stream.log"

# --- PHASE 5: Bayesian Network Learning (BaNDyT) ---
echo "Running Phase 5: Discretized Bayesian Structure Learning..."
python 05_bandyt.py --run_id "$RUN_ID"
log_phase "PHASE_5_COMPLETE: BaNDyT causal network estimation finished."

# --- PHASE 6: Network Topology Analysis ---
echo "Running Phase 6: Headless Network Topology Analysis..."
python 06_network_analysis.py --run_id "$RUN_ID" --pdb "../input/$MATCHED_PDB"
log_phase "PHASE_6_COMPLETE: Topological metrics calculated."

# --- PHASE 7: PyMOL 3D Network Mapping ---
echo "Running Phase 7: PyMOL 3D Structural Network Automation..."
python 07_pymol_network.py \
    --run_id "$RUN_ID" \
    --pdb "$MATCHED_PDB" \
    --threshold 0.78
log_phase "PHASE_7_COMPLETE: 3D PyMOL PML visual script rendered and cached."

# --- PHASE 8: Active Site Local Network Isolation ---
echo "Running Phase 8: Active Site Network Isolation & PyMOL Generation..."
python 08_activecon_pymol.py \
    --run_id "$RUN_ID" \
    --pdb "$MATCHED_PDB" \
    --threshold 0.78 \
    --ligand "$LIGAND"
log_phase "PHASE_8_COMPLETE: Active site localized network maps built."

deactivate

log_phase "SUCCESS: Replica ${REP_NUM} completed cleanly."
