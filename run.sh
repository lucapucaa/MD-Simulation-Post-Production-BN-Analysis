#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=dyrk1a_bandyt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/slurm_%j.log

set -e

# 1. Absolute Path Anchor
PROJECT_ROOT="/home/lchill/work/pproduct_analysis"
cd "$PROJECT_ROOT/pipeline"

# ==========================================
# ⚙️ SYSTEM VARIABLES
# ==========================================
RAW_PDB="7o7k_ATP.pdb"      # Used for Phases 2, 3, 4, 6, 7, 8
RAW_TPR="WT_rep_1.tpr"         # Used ONLY for Phase 1 (GROMACS PBC cleaning)
RAW_XTC="WT_rep_1.xtc"
LIGAND="ATP"
RUN_ID="replica_1"

# Dynamically name cleaned output trajectory based on RUN_ID
CLEAN_XTC="${RUN_ID}_final.xtc"
# ==========================================

# 2. Setup Dedicated Logs Directory
LOGS_DIR="$PROJECT_ROOT/logs"
mkdir -p "$LOGS_DIR"

STATUS_LOG="${LOGS_DIR}/${RUN_ID}_pipeline_status.log"

# 💡 FORCED OVERRIDE: Always create a fresh log file
echo "==================================================" > "$STATUS_LOG"
echo "    PARALLEL PIPELINE LEDGER FOR: ${RUN_ID}" >> "$STATUS_LOG"
echo "    (MODE: FORCED OVERRIDE - ALL PHASES RERUN)" >> "$STATUS_LOG"
echo "==================================================" >> "$STATUS_LOG"

log_phase() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ✅ $1" >> "$STATUS_LOG"
}

failure_handler() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ❌ FATAL ERROR: Pipeline crashed on or around line $1" >> "$STATUS_LOG"
    echo "Check primary log ($PROJECT_ROOT/logs/slurm_${SLURM_JOB_ID}.log) or phase stream logs for details." >> "$STATUS_LOG"
}
trap 'failure_handler $LINENO' ERR

# ==============================================================================
# PIPELINE EXECUTION ENGINE
# ==============================================================================

log_phase "START: Initializing environments."

echo "=== Loading Global Cluster Python & GROMACS Environments ==="
if command -v module &> /dev/null; then
    module load Python/3.9.5-GCCcore-10.3.0 2>/dev/null || true
    module load Gromacs/2022.1-Container 2>/dev/null || true
fi

if [ -f ~/dyrk1a_venv/bin/activate ]; then
    source ~/dyrk1a_venv/bin/activate
fi


# --- PHASE 1: Trajectory Boundary Cleaning ---
echo "Starting Phase 1: Cleaning Trajectory..."
./01_clean_traj.sh "$RAW_TPR" "$RAW_XTC" "$LIGAND" "$RUN_ID"
log_phase "PHASE_1_COMPLETE: GROMACS trajectory boundary unwrapping finished."


# --- PHASE 2: Quality Metrics (RMSD, RMSF, Rg, DSSP) ---
echo "Running Phase 2: Quality Validation..."
./02_quality_metrics.sh "$RAW_PDB" "$CLEAN_XTC" "$LIGAND"
log_phase "PHASE_2_COMPLETE: Structural quality evaluations calculated."


# --- PARALLEL BLOCK (PHASE 3 & PHASE 4) ---
echo "🚀 Spawning parallel processing block..."

echo "-> Phase 3 (GetContacts Fingerprinting) running in background..."
./03_contacts.sh "$RAW_PDB" "$CLEAN_XTC" "$RUN_ID" "$LIGAND" > "${LOGS_DIR}/${RUN_ID}_phase3_stream.log" 2>&1 &
PID_PHASE3=$!

echo "-> Phase 4 (Per-Residue Energy Mapping) running in background..."
python -u 04_perres_energy.py \
    --run_id "$RUN_ID" \
    --pdb "$RAW_PDB" \
    --xtc "$CLEAN_XTC" \
    --ligand "$LIGAND" \
    --cutoff 12.0 \
    --cpus 16 > "${LOGS_DIR}/${RUN_ID}_phase4_stream.log" 2>&1 &
PID_PHASE4=$!

# Synchronize execution paths for active background processes
echo "⏳ Syncing Phase 3 and Phase 4 threads..."
wait $PID_PHASE3 || { echo "❌ Phase 3 (GetContacts) async task threw an error."; exit 1; }
log_phase "PHASE_3_COMPLETE: GetContacts geometric fingerprinting finished."
rm -f "${LOGS_DIR}/${RUN_ID}_phase3_stream.log"

wait $PID_PHASE4 || { echo "❌ Phase 4 (Energy Mapping) async task threw an error."; exit 1; }
log_phase "PHASE_4_COMPLETE: Continuous interaction energy extraction complete."
rm -f "${LOGS_DIR}/${RUN_ID}_phase4_stream.log"


# --- PHASE 5: Bayesian Network Learning (BaNDyT) ---
echo "Running Phase 5: Discretized Bayesian Structure Learning..."
python 05_bandyt.py --run_id "$RUN_ID"
log_phase "PHASE_5_COMPLETE: BaNDyT causal network estimation finished."


# --- PHASE 6: Network Topology Analysis ---
echo "Running Phase 6: Headless Network Topology Analysis..."
python 06_network_analysis.py --run_id "$RUN_ID" --pdb "../input/$RAW_PDB"
log_phase "PHASE_6_COMPLETE: Topological metrics calculated and final graph drafted."


# --- PHASE 7: PyMOL 3D Network Mapping ---
echo "Running Phase 7: PyMOL 3D Structural Network Automation..."
python 07_pymol_network.py \
    --run_id "$RUN_ID" \
    --pdb "$RAW_PDB" \
    --threshold 0.78
log_phase "PHASE_7_COMPLETE: 3D PyMOL PML visual script rendered and cached."


# --- PHASE 8: Active Site Local Network Isolation ---
echo "Running Phase 8: Active Site Network Isolation & PyMOL Generation..."
python 08_activecon_pymol.py \
    --run_id "$RUN_ID" \
    --pdb "$RAW_PDB" \
    --threshold 0.78 \
    --ligand "$LIGAND"
log_phase "PHASE_8_COMPLETE: Active site localized network maps built and PML written."


if command -v deactivate &> /dev/null; then
    deactivate
fi

log_phase "SUCCESS: Whole system pipeline execution concluded cleanly. Resources released."
echo "✅ All phases executed successfully."
