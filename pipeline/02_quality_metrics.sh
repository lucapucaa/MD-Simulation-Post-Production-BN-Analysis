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
# PATH & ARGUMENT CONFIGURATION
# ==============================================================================
DATA_DIR="../data"
INPUT_DIR="../input"
OUTPUT_DIR="../results/figures"

PDB_PATH="${INPUT_DIR}/$1"
# Note: Since $2 is already RUN_ID_final.xtc, we look for it in DATA_DIR
TRAJ_PATH="${DATA_DIR}/$2"  
LIGAND_NAME="${3:-ATP}"
RUN_ID="${4:-replica_1}" # 🎯 Newly added to prevent race conditions

# 🎯 Define DH-Box Residue Range here
DH_RES="137-153"

# 🎯 RACE CONDITION FIX: Unique Temporary Files & Outputs
NDX_TEMP="${DATA_DIR}/${RUN_ID}_temp_quality.ndx"

mkdir -p "$OUTPUT_DIR"

echo "=================================================="
echo "📊 Phase 2: Structural Quality Controls & Analysis"
echo " Run ID       : ${RUN_ID}"
echo " Topology PDB : ${PDB_PATH}"
echo " Trajectory   : ${TRAJ_PATH}"
echo " Ligand Name  : ${LIGAND_NAME}"
echo " DH-Box Range : ${DH_RES}"
echo " Figures Dir  : ${OUTPUT_DIR}"
echo "=================================================="

# Create index groups: 
# 1. 'Complex' (Protein + Ligand)
# 2. 'CA_No_DH' (C-alpha atoms excluding the DH-box)
echo -e "1 | r ${LIGAND_NAME}\nname 14 Complex\n3 & ! r ${DH_RES}\nname 15 CA_No_DH\nq" | gmx make_ndx -f "$PDB_PATH" -o "$NDX_TEMP"

# 1a. Calculate Normal RMSD (Whole Protein C-alpha)
echo -e "C-alpha\nC-alpha" | gmx rms -s "$PDB_PATH" -f "$TRAJ_PATH" -n "$NDX_TEMP" -o "${DATA_DIR}/${RUN_ID}_rmsd.xvg"

# 1b. Calculate RMSD without DH-box (Fit & Compute on C-alpha excluding DH-box)
echo -e "CA_No_DH\nCA_No_DH" | gmx rms -s "$PDB_PATH" -f "$TRAJ_PATH" -n "$NDX_TEMP" -o "${DATA_DIR}/${RUN_ID}_rmsd_no_dh.xvg"

# 2. Calculate RMSF (Per-residue C-alpha atoms)
echo -e "C-alpha" | gmx rmsf -s "$PDB_PATH" -f "$TRAJ_PATH" -n "$NDX_TEMP" -o "${DATA_DIR}/${RUN_ID}_rmsf_per_res.xvg" -res

# 3. Calculate Radius of Gyration targeting 'Complex' index group
echo -e "Complex" | gmx gyrate -s "$PDB_PATH" -f "$TRAJ_PATH" -n "$NDX_TEMP" -o "${DATA_DIR}/${RUN_ID}_gyrate.xvg"

echo "=== Plotting Results (All Distance Metrics in Ångströms) ==="

# Inline Python plotting engine
python - <<EOF
import matplotlib
matplotlib.use('Agg')  # Headless mode
import matplotlib.pyplot as plt
import numpy as np
import mdtraj as md
import os
from matplotlib.colors import ListedColormap

# --- Part 1A: C-alpha RMSD Plotting (WITH DH-box) ---
try:
    rmsd_xvg = "${DATA_DIR}/${RUN_ID}_rmsd.xvg"
    if os.path.exists(rmsd_xvg):
        data_rmsd = np.loadtxt(rmsd_xvg, comments=["@", "#"])
        time_ns = data_rmsd[:, 0] / 1000.0 if data_rmsd[:, 0].max() > 1000 else data_rmsd[:, 0]
        rmsd_angstrom = data_rmsd[:, 1] * 10.0
        
        window = 10
        rmsd_avg = np.convolve(rmsd_angstrom, np.ones(window)/window, mode='same')

        plt.figure(figsize=(8, 4))
        plt.plot(time_ns, rmsd_angstrom, color="#D81B60", alpha=0.35, linewidth=1.0, label="Raw Trajectory")
        plt.plot(time_ns, rmsd_avg, color="#880E4F", alpha=1.0, linewidth=2.0, label="Running Average")

        plt.title("RMSD: Trajectory Stability and Robustness (with DH-box)", fontsize=12, fontweight="bold")
        plt.xlabel("Time (ns)", fontsize=10)
        plt.ylabel("RMSD (Å)", fontsize=10)
        plt.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="none")
        plt.grid(True, linestyle=":", alpha=0.6)
        plt.xlim(time_ns[0], time_ns[-1])
        plt.tight_layout()
        plt.savefig("${OUTPUT_DIR}/${RUN_ID}_rmsd_with_dh.png", dpi=300)
        plt.close()
        print("✅ Saved plot: ${OUTPUT_DIR}/${RUN_ID}_rmsd_with_dh.png")
except Exception as e:
    print(f"❌ RMSD (With DH-Box) plotting error: {e}")

# --- Part 1B: C-alpha RMSD Plotting (NO DH-box) ---
try:
    rmsd_no_dh_xvg = "${DATA_DIR}/${RUN_ID}_rmsd_no_dh.xvg"
    if os.path.exists(rmsd_no_dh_xvg):
        data_no_dh = np.loadtxt(rmsd_no_dh_xvg, comments=["@", "#"])
        time_ns = data_no_dh[:, 0] / 1000.0 if data_no_dh[:, 0].max() > 1000 else data_no_dh[:, 0]
        no_dh_angstrom = data_no_dh[:, 1] * 10.0
        
        window = 10
        no_dh_avg = np.convolve(no_dh_angstrom, np.ones(window)/window, mode='same')

        plt.figure(figsize=(8, 4))
        plt.plot(time_ns, no_dh_angstrom, color="#1565C0", alpha=0.35, linewidth=1.0, label="Raw Trajectory")
        plt.plot(time_ns, no_dh_avg, color="#0D47A1", alpha=1.0, linewidth=2.0, label="Running Average")

        plt.title("RMSD: Trajectory Stability and Robustness", fontsize=12, fontweight="bold")
        plt.xlabel("Time (ns)", fontsize=10)
        plt.ylabel("RMSD (Å)", fontsize=10)
        plt.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="none")
        plt.grid(True, linestyle=":", alpha=0.6)
        plt.xlim(time_ns[0], time_ns[-1])
        plt.tight_layout()
        plt.savefig("${OUTPUT_DIR}/${RUN_ID}_rmsd_no_dh.png", dpi=300)
        plt.close()
        print("✅ Saved plot: ${OUTPUT_DIR}/${RUN_ID}_rmsd_no_dh.png")
except Exception as e:
    print(f"❌ RMSD (No DH-Box) plotting error: {e}")


# --- Part 2: Per-Residue RMSF Plotting (Ångströms) ---
try:
    rmsf_xvg = "${DATA_DIR}/${RUN_ID}_rmsf_per_res.xvg"
    if os.path.exists(rmsf_xvg):
        data_rmsf = np.loadtxt(rmsf_xvg, comments=["@", "#"])
        res_num = data_rmsf[:, 0]
        rmsf_angstrom = data_rmsf[:, 1] * 10.0  # nm -> Å

        plt.figure(figsize=(8, 4))
        plt.plot(res_num, rmsf_angstrom, color="#D81B60", linewidth=1.5, label="C-alpha Fluctuation")

        plt.title("Per-Residue Structural Fluctuation (RMSF)", fontsize=12, fontweight="bold")
        plt.xlabel("Residue Number", fontsize=10)
        plt.ylabel("RMSF (Å)", fontsize=10)
        plt.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="none")
        plt.grid(True, linestyle=":", alpha=0.6)
        plt.xlim(res_num[0], res_num[-1])
        plt.tight_layout()
        plt.savefig("${OUTPUT_DIR}/${RUN_ID}_rmsf.png", dpi=300)
        plt.close()
        print("✅ Saved plot: ${OUTPUT_DIR}/${RUN_ID}_rmsf.png")
except Exception as e:
    print(f"❌ RMSF plotting error: {e}")


# --- Part 3: Radius of Gyration Plotting (Ångströms) ---
try:
    xvg_path = "${DATA_DIR}/${RUN_ID}_gyrate.xvg"
    if os.path.exists(xvg_path):
        data = np.loadtxt(xvg_path, comments=["@", "#"])
        time_ns = data[:, 0] / 1000.0 if data[:, 0].max() > 1000 else data[:, 0]

        # Extract only the Total Rg (column index 1)
        rg_total = data[:, 1] * 10.0  # nm -> Å

        plt.figure(figsize=(8, 4))
        plt.plot(time_ns, rg_total, label="Total Rg", color="black", linewidth=2)

        plt.title("Protein Compactness Profile (Radius of Gyration)", fontsize=12, fontweight="bold")
        plt.xlabel("Time (ns)", fontsize=10)
        plt.ylabel("Radius of Gyration (Å)", fontsize=10)
        plt.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="none")
        plt.grid(True, linestyle=":", alpha=0.6)
        plt.tight_layout()
        plt.savefig("${OUTPUT_DIR}/${RUN_ID}_radius_of_gyration.png", dpi=300)
        plt.close()
        print("✅ Saved plot: ${OUTPUT_DIR}/${RUN_ID}_radius_of_gyration.png")
except Exception as e:
    print(f"❌ Gyration plotting error: {e}")


# --- Part 4: MDTraj DSSP Secondary Structure Timeline ---
try:
    pdb_p = "${PDB_PATH}"
    xtc_p = "${TRAJ_PATH}"
    
    traj = md.load(xtc_p, top=pdb_p, stride=10)
    dssp = md.compute_dssp(traj, simplified=True)
    
    num_dssp = np.zeros(dssp.shape)
    num_dssp[dssp == 'H'] = 1  
    num_dssp[dssp == 'E'] = 2  
    num_dssp[dssp == 'C'] = 0  
    
    time_axis = traj.time / 1000.0 if traj.time.max() > 1000 else traj.time
    residues = [res.resSeq for res in traj.topology.residues if res.is_protein]
    num_dssp = num_dssp[:, :len(residues)]
    
    plt.figure(figsize=(10, 5))
    custom_cmap = ListedColormap(['#e0e0e0', '#d62728', '#1f77b4'])
    
    plt.imshow(num_dssp.T, aspect='auto', origin='lower', cmap=custom_cmap,
               extent=[time_axis[0], time_axis[-1], residues[0], residues[-1]])
    
    plt.title("Secondary Structure Evolutionary Timeline (DSSP)", fontsize=12, fontweight="bold")
    plt.xlabel("Time (ns)", fontsize=10)
    plt.ylabel("Residue ID", fontsize=10)
    
    cbar = plt.colorbar(ticks=[0.33, 1.0, 1.66])
    cbar.set_ticklabels(['Loop / Coil', 'Alpha-Helix', 'Beta-Sheet'])
    cbar.ax.tick_params(labelsize=9)
    
    plt.tight_layout()
    plt.savefig("${OUTPUT_DIR}/${RUN_ID}_secondary_structure_timeline.png", dpi=300)
    plt.close()
    print("✅ Saved plot: ${OUTPUT_DIR}/${RUN_ID}_secondary_structure_timeline.png")
    
except Exception as e:
    print(f"❌ DSSP Timeline processing error: {e}")
EOF

# Housekeeping
rm -f "$NDX_TEMP"
echo "=================================================="
echo "✅ Phase 2 Complete for ${RUN_ID}!"
echo "=================================================="
