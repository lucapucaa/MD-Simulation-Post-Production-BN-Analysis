import os
import argparse
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

def generate_dssp_heatmap():
    # 1. Setup paths matching your directory workflow
    data_dir = "../data"
    output_dir = "../results/figures"
    os.makedirs(output_dir, exist_ok=True)

    pdb_path = os.path.join(data_dir, "prot-lig.pdb")
    xtc_path = os.path.join(data_dir, "final.xtc")

    print("Loading trajectory into MDTraj...")
    # Load trajectory and strip out water/ions automatically to keep it light
    traj = md.load(xtc_path, top=pdb_path)
    protein_indices = traj.topology.select("protein")
    protein_traj = traj.atom_slice(protein_indices)

    print("Computing DSSP assignments (Simplified Scheme)...")
    # simplified=True yields 3 categories: 'H' (Helix), 'E' (Strand), 'C' (Coil)
    dssp_matrix = md.compute_dssp(protein_traj, simplified=True)

    # 2. Map character data to integers for heatmap rendering
    # Helix = 2, Strand = 1, Coil = 0
    char_to_int = {'H': 2, 'E': 1, 'C': 0, 'NA': 0}
    numeric_dssp = np.vectorize(char_to_int.get)(dssp_matrix)

    # Convert frames axis to nanoseconds (assuming 10ps or 100ps coordinate saving frequency)
    # Adjust the multiplier if your frame stride is different
    time_ns = protein_traj.time / 1000.0 
    
    # Extract structural residue numbering bounds
    residues = [r.resSeq for r in protein_traj.topology.residues]
    
    print("Generating publication-ready heatmap plot...")
    plt.figure(figsize=(11, 5))
    
    # Create a clean, categorical colormap: Grey for Coil, Black/Blue for Sheet, Red for Helix
    custom_cmap = ListedColormap(['#d3d3d3', '#1f77b4', '#d62728'])

    # Render matrix. Transpose (.T) so residues are on the Y-axis and Time is on the X-axis
    extent = [time_ns[0], time_ns[-1], residues[0], residues[-1]]
    plt.imshow(numeric_dssp.T, aspect='auto', origin='lower', cmap=custom_cmap, extent=extent)

    # 3. Labeling and styling matching your suite constraints
    plt.title("DYRK1A Secondary Structure Timeline", fontsize=14, fontweight='bold')
    plt.xlabel("Time (ns)", fontsize=12)
    plt.ylabel("Residue Number", fontsize=12)

    # Setup categorical colorbar legend
    cbar = plt.colorbar(ticks=[0.33, 1.0, 1.66])
    cbar.ax.set_yticklabels(['Coil/Turn', 'Beta Strand', 'Alpha Helix'], fontsize=10)
    cbar.set_label('Structural State', fontsize=11, fontweight='bold')

    plt.tight_layout()
    output_path = os.path.join(output_dir, "dssp_timeline.png")
    plt.savefig(output_path, dpi=300)
    print(f"Success! Structural timeline saved to: {output_path}")

if __name__ == "__main__":
    generate_dssp_heatmap()import os
import argparse
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

def generate_dssp_heatmap():
    # 1. Setup paths matching your directory workflow
    data_dir = "../data"
    output_dir = "../results/figures"
    os.makedirs(output_dir, exist_ok=True)

    pdb_path = os.path.join(data_dir, "prot-lig.pdb")
    xtc_path = os.path.join(data_dir, "final.xtc")

    print("Loading trajectory into MDTraj...")
    # Load trajectory and strip out water/ions automatically to keep it light
    traj = md.load(xtc_path, top=pdb_path)
    protein_indices = traj.topology.select("protein")
    protein_traj = traj.atom_slice(protein_indices)

    print("Computing DSSP assignments (Simplified Scheme)...")
    # simplified=True yields 3 categories: 'H' (Helix), 'E' (Strand), 'C' (Coil)
    dssp_matrix = md.compute_dssp(protein_traj, simplified=True)

    # 2. Map character data to integers for heatmap rendering
    # Helix = 2, Strand = 1, Coil = 0
    char_to_int = {'H': 2, 'E': 1, 'C': 0, 'NA': 0}
    numeric_dssp = np.vectorize(char_to_int.get)(dssp_matrix)

    # Convert frames axis to nanoseconds (assuming 10ps or 100ps coordinate saving frequency)
    # Adjust the multiplier if your frame stride is different
    time_ns = protein_traj.time / 1000.0 
    
    # Extract structural residue numbering bounds
    residues = [r.resSeq for r in protein_traj.topology.residues]
    
    print("Generating publication-ready heatmap plot...")
    plt.figure(figsize=(11, 5))
    
    # Create a clean, categorical colormap: Grey for Coil, Black/Blue for Sheet, Red for Helix
    custom_cmap = ListedColormap(['#d3d3d3', '#1f77b4', '#d62728'])

    # Render matrix. Transpose (.T) so residues are on the Y-axis and Time is on the X-axis
    extent = [time_ns[0], time_ns[-1], residues[0], residues[-1]]
    plt.imshow(numeric_dssp.T, aspect='auto', origin='lower', cmap=custom_cmap, extent=extent)

    # 3. Labeling and styling matching your suite constraints
    plt.title("DYRK1A Secondary Structure Timeline", fontsize=14, fontweight='bold')
    plt.xlabel("Time (ns)", fontsize=12)
    plt.ylabel("Residue Number", fontsize=12)

    # Setup categorical colorbar legend
    cbar = plt.colorbar(ticks=[0.33, 1.0, 1.66])
    cbar.ax.set_yticklabels(['Coil/Turn', 'Beta Strand', 'Alpha Helix'], fontsize=10)
    cbar.set_label('Structural State', fontsize=11, fontweight='bold')

    plt.tight_layout()
    output_path = os.path.join(output_dir, "dssp_timeline.png")
    plt.savefig(output_path, dpi=300)
    print(f"Success! Structural timeline saved to: {output_path}")

if __name__ == "__main__":
    generate_dssp_heatmap()
