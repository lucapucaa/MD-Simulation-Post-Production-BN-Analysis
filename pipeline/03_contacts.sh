#!/bin/bash
set -e

# ==============================================================================
# PHASE 2.5: FRAME-BY-FRAME MOLECULAR CONTACT EXTRACTION (GetContacts)
# ==============================================================================

RAW_PDB="$1"
RAW_XTC="$2"
RUN_ID="${3:-replica_1}"
LIGAND="${4:-ATP}"

INPUT_DIR="../input"
DATA_DIR="../data"

mkdir -p "${DATA_DIR}"

RAW_PDB_PATH="${INPUT_DIR}/${RAW_PDB}"
RAW_XTC_PATH="${DATA_DIR}/${RAW_XTC}"

# Check input/ fallback for trajectory if not in data/
if [ ! -f "$RAW_XTC_PATH" ] && [ -f "${INPUT_DIR}/${RAW_XTC}" ]; then
    RAW_XTC_PATH="${INPUT_DIR}/${RAW_XTC}"
fi

# Sanity Checks
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "❌ Usage: $0 <pdb_file> <xtc_file> [run_id] [ligand_resname]"
    echo "   Example: $0 7o7k_ATP.pdb final.xtc rep1 ATP"
    exit 1
fi

if [ ! -f "$RAW_PDB_PATH" ]; then
    echo "❌ Error: Topology PDB file not found: ${RAW_PDB_PATH}"
    exit 1
fi

if [ ! -f "$RAW_XTC_PATH" ]; then
    echo "❌ Error: Trajectory XTC file not found: ${RAW_XTC_PATH}"
    exit 1
fi

# Convert to absolute paths
PDB_PATH=$(realpath "${RAW_PDB_PATH}")
TRAJ_PATH=$(realpath "${RAW_XTC_PATH}")
OUTPUT_TSV="$(realpath "${DATA_DIR}")/${RUN_ID}_contacts.tsv"

echo "=================================================="
echo "🔬 Phase 2.5: Molecular Contact Extraction"
echo " Target System : ${RUN_ID}"
echo " Ligand        : ${LIGAND}"
echo " PDB Path      : ${PDB_PATH}"
echo " Traj Path     : ${TRAJ_PATH}"
echo " Output TSV    : ${OUTPUT_TSV}"
echo "=================================================="

# 🔄 DYNAMIC ENVIRONMENT DETECTION
USING_VENV=false
USING_CONDA=false

if command -v module &> /dev/null; then
    echo "🖥️ Cluster environment detected! Loading cluster toolchains..."
    module load Python/3.9.5-GCCcore-10.3.0 2>/dev/null || true
    if [ -f ~/dyrk1a_venv/bin/activate ]; then
        source ~/dyrk1a_venv/bin/activate
        USING_VENV=true
    fi
elif command -v conda &> /dev/null; then
    echo "💻 Local/Miniconda environment detected!"
    eval "$(conda shell.bash hook)"
    if conda info --envs | grep -q "contacts_env"; then
        conda activate contacts_env
    else
        echo "ℹ️ 'contacts_env' not found, using active conda environment."
    fi
    USING_CONDA=true
fi

# GetContacts Directory Path
GETCONTACTS_DIR="/home/lchill/work/getcontacts"

if [ ! -d "$GETCONTACTS_DIR" ]; then
    echo "❌ Error: GetContacts installation not found at '${GETCONTACTS_DIR}'"
    exit 1
fi

# 🧠 CPU Thread Detection (Capped at 4/8 to prevent memory exhaustion)
DETECTED_CORES=$(nproc 2>/dev/null || echo 4)
if [ "$DETECTED_CORES" -gt 8 ]; then
    CPU_CORES=4
    echo "⚠️ System has ${DETECTED_CORES} cores. Capping at ${CPU_CORES} to prevent OOM crashes."
else
    CPU_CORES=$DETECTED_CORES
fi

echo "🚀 Spawning tracking engine using ${CPU_CORES} CPU cores..."

# Run frame-by-frame contact scan
python "${GETCONTACTS_DIR}/get_dynamic_contacts.py" \
    --topology "$PDB_PATH" \
    --trajectory "$TRAJ_PATH" \
    --output "$OUTPUT_TSV" \
    --cores "$CPU_CORES" \
    --sele "protein" \
    --sele2 "resname $LIGAND" \
    --ligand "resname $LIGAND" \
    --itypes all

echo "=================================================="
echo "🎉 SUCCESS: Contact profile saved to: ${OUTPUT_TSV}"
echo "=================================================="

# Cleanup environment state
if [ "$USING_VENV" = true ]; then
    deactivate 2>/dev/null || true
elif [ "$USING_CONDA" = true ]; then
    conda deactivate 2>/dev/null || true
fi
