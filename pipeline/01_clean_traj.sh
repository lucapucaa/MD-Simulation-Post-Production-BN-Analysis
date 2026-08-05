#!/bin/bash
set -e

# ==============================================================================
# ENVIRONMENT RESOLUTION (Find GROMACS)
# ==============================================================================
if ! command -v gmx &> /dev/null; then
    CONDA_GMX="/home/lchill/miniconda3/envs/dyrk1a/bin/gmx"
    CONDA_AVX_GMX="/home/lchill/miniconda3/envs/dyrk1a/bin.AVX2_256/gmx"

    if [ -f "$CONDA_AVX_GMX" ]; then
        export PATH="/home/lchill/miniconda3/envs/dyrk1a/bin.AVX2_256:$PATH"
    elif [ -f "$CONDA_GMX" ]; then
        export PATH="/home/lchill/miniconda3/envs/dyrk1a/bin:$PATH"
    elif command -v module &> /dev/null; then
        module load Gromacs 2>/dev/null || true
    fi
fi

if ! command -v gmx &> /dev/null; then
    echo "❌ Error: 'gmx' command not found. Please activate your environment."
    exit 1
fi

# ==============================================================================
# PHASE 1: TRAJECTORY CLEANING (TPR MODE)
# ==============================================================================

RAW_TPR="$1"            # e.g., "topol.tpr" (or "md.tpr")
RAW_XTC="$2"            # e.g., "rep_1.xtc"
LIGAND_NAME="${3:-ATP}"  # Default ligand residue name is ATP
RUN_ID="${4:-replica_1}" # Default run ID is replica_1

INPUT_DIR="../input"
DATA_DIR="../data"

mkdir -p "$DATA_DIR"

TPR_PATH="${INPUT_DIR}/${RAW_TPR}"
XTC_PATH="${DATA_DIR}/${RAW_XTC}"

if [ ! -f "$XTC_PATH" ] && [ -f "${INPUT_DIR}/${RAW_XTC}" ]; then
    XTC_PATH="${INPUT_DIR}/${RAW_XTC}"
fi

# Sanity Checks
if [ ! -f "$TPR_PATH" ]; then
    echo "❌ Error: Input TPR topology not found at ${TPR_PATH}"
    exit 1
fi

if [ ! -f "$XTC_PATH" ]; then
    echo "❌ Error: Input XTC trajectory not found at ${XTC_PATH}"
    exit 1
fi

# 🎯 RACE CONDITION FIX: Unique Temporary Files 
NDX_TEMP="${DATA_DIR}/${RUN_ID}_temp_clean.ndx"
TEMP_NOJUMP="${DATA_DIR}/${RUN_ID}_temp_nojump.xtc"
CLEAN_XTC="${DATA_DIR}/${RUN_ID}_final.xtc"

echo "=================================================="
echo "🧹 Phase 1: Trajectory Boundary Cleaning (TPR Mode)"
echo " TPR Topology : ${TPR_PATH}"
echo " Input Traj   : ${XTC_PATH}"
echo " Ligand Name  : ${LIGAND_NAME}"
echo " Run ID       : ${RUN_ID}"
echo " Output Traj  : ${CLEAN_XTC}"
echo "=================================================="

# 1. Create index group merging Group 1 (Protein) and residue name into 'Protein_Ligand'
echo -e "1 | r ${LIGAND_NAME}\nname 17 Protein_Ligand\nq" | gmx make_ndx -f "$TPR_PATH" -o "$NDX_TEMP"

# 2. Step 1: Remove coordinate jumps across periodic boundaries across the FULL system
echo -e "Protein\nSystem" | gmx trjconv \
    -s "$TPR_PATH" \
    -f "$XTC_PATH" \
    -n "$NDX_TEMP" \
    -o "$TEMP_NOJUMP" \
    -pbc nojump \
    -center

# 3. Step 2: Cluster ligand around protein and output ONLY Protein + Ligand
echo -e "Protein\nProtein\nProtein_Ligand" | gmx trjconv \
    -s "$TPR_PATH" \
    -f "$TEMP_NOJUMP" \
    -n "$NDX_TEMP" \
    -o "$CLEAN_XTC" \
    -pbc cluster \
    -center

# Housekeeping
rm -f "$NDX_TEMP" "$TEMP_NOJUMP"

echo "=================================================="
echo "✅ Phase 1 Complete!"
echo " Saved filtered trajectory to: ${CLEAN_XTC}"
echo "=================================================="
