import os
import sys
import argparse
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
import multiprocessing

def get_atom_vdw_parameters(atom):
    """Maps standard elements to universal Amber/GAFF-like Van der Waals parameters."""
    if hasattr(atom, 'element') and atom.element:
        element = atom.element.upper()
    else:
        name_upper = atom.name.upper()
        element = name_upper[:2] if name_upper[:2] in ['CL', 'BR'] else name_upper[0]
    
    vdw_table = {
        'C': (1.90, 0.086),
        'N': (1.82, 0.170),
        'O': (1.66, 0.210),
        'S': (2.00, 0.250),
        'H': (1.10, 0.015),
        'P': (2.10, 0.200),
        'F': (1.75, 0.061),
        'CL': (1.95, 0.227),
        'BR': (2.10, 0.320),
        'I': (2.35, 0.400)
    }
    return vdw_table.get(element, (1.80, 0.100))

def compute_all_pair_energies_worker(args):
    """
    Worker process evaluating all unique pairwise node interactions (Res-Res & Res-Lig)
    across an assigned chunk of trajectory frames.
    """
    pdb_path, xtc_path, cached_nodes, pair_indices, frame_indices, cutoff = args
    
    u_local = mda.Universe(pdb_path, xtc_path)
    node_atom_groups = [u_local.atoms[node['indices']] for node in cached_nodes]
    
    num_pairs = len(pair_indices)
    results = []
    
    for f_idx in frame_indices:
        u_local.trajectory[f_idx]
        frame_energies = np.zeros(num_pairs, dtype=np.float32)
        
        for p_idx, (i, j) in enumerate(pair_indices):
            group_i = node_atom_groups[i]
            group_j = node_atom_groups[j]
            
            # Distance cutoff filter using Center of Mass (COM)
            com_dist = np.linalg.norm(group_i.center_of_mass() - group_j.center_of_mass())
            if com_dist > cutoff:
                continue  # Skip distant pairs (> 12 Å)
                
            pos_i = group_i.positions
            pos_j = group_j.positions
            
            dists = distance_array(pos_i, pos_j)
            dists = np.maximum(dists, 0.1)  # Safeguard against zero division
            
            # Pairwise Lorentz-Berthelot parameters
            sig_i, eps_i = cached_nodes[i]['sigmas'], cached_nodes[i]['epsilons']
            sig_j, eps_j = cached_nodes[j]['sigmas'], cached_nodes[j]['epsilons']
            
            R_sig = (sig_i[:, np.newaxis] + sig_j[np.newaxis, :]) / 2.0
            R_eps = np.sqrt(eps_i[:, np.newaxis] * eps_j[np.newaxis, :])
            
            # Lennard-Jones 12-6 calculation
            inv_r6 = (R_sig / dists) ** 6
            inv_r12 = inv_r6 ** 2
            vdw_energy = 4 * R_eps * (inv_r12 - inv_r6)
            
            vdw_energy = np.clip(vdw_energy, -20.0, 50.0)
            frame_energies[p_idx] = np.sum(vdw_energy)
            
        results.append((f_idx, frame_energies))
        
    return results

def main():
    parser = argparse.ArgumentParser(description="All-Against-All Residue & Ligand Pairwise Interaction Energy Calculator")
    parser.add_argument("--run_id", required=True, help="Prefix name for the run")
    parser.add_argument("--pdb", required=True, help="Filename of system PDB topology")
    parser.add_argument("--xtc", required=True, help="Filename of trajectory XTC")
    parser.add_argument("--ligand", required=True, help="Residue name of ligand (e.g., ATP)")
    parser.add_argument("--cutoff", type=float, default=12.0, help="COM distance cutoff in Ångströms (default: 12.0 Å)")
    parser.add_argument("--cpus", type=int, default=16, help="Number of CPU cores")
    args = parser.parse_args()

    input_dir = "../input"
    data_dir = "../data"
    pdb_path = os.path.join(input_dir, args.pdb)
    xtc_path = os.path.join(data_dir, args.xtc)
    
    # Matches the default file name expected upstream by BaNDyT (Phase 5)
    output_csv = os.path.join(data_dir, f"{args.run_id}_interaction_energies.csv")

    if not os.path.exists(pdb_path) or not os.path.exists(xtc_path):
        print(f"❌ Error: Missing input files!\nPDB: {pdb_path}\nXTC: {xtc_path}")
        sys.exit(1)

    print(f"=== Phase 4: All-Against-All Interaction Energy Extraction for {args.run_id} ===")
    u = mda.Universe(pdb_path, xtc_path)
    
    # 1. Build Node Inventory (Protein Residues + Ligand)
    protein_atoms = u.select_atoms("protein")
    protein_residues = sorted(list(protein_atoms.residues), key=lambda r: r.resid)
    
    ligand_atoms = u.select_atoms(f"resname {args.ligand}")
    if len(ligand_atoms) == 0:
        print(f"❌ Error: Ligand residue '{args.ligand}' not found!")
        sys.exit(1)

    nodes = []
    for res in protein_residues:
        nodes.append({
            'label': f"A:{res.resname}:{res.resid}",
            'atoms': res.atoms
        })
    ligand_res = ligand_atoms.residues[0]
    nodes.append({
        'label': f"A:{args.ligand}:{ligand_res.resid}",
        'atoms': ligand_atoms
    })

    num_nodes = len(nodes)
    print(f"Identified {num_nodes} total interacting nodes ({len(protein_residues)} residues + 1 ligand).")

    # 2. Pre-calculate atomic VDW parameters
    cached_nodes = []
    for node in nodes:
        atoms = node['atoms']
        sigmas = np.array([get_atom_vdw_parameters(a)[0] for a in atoms])
        epsilons = np.array([get_atom_vdw_parameters(a)[1] for a in atoms])
        cached_nodes.append({
            'label': node['label'],
            'indices': atoms.indices,
            'sigmas': sigmas,
            'epsilons': epsilons
        })

    # 3. Generate all unique pairwise combinations
    pair_indices = []
    pair_labels = []
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            pair_indices.append((i, j))
            pair_labels.append(f"{cached_nodes[i]['label']}_{cached_nodes[j]['label']}")

    num_pairs = len(pair_indices)
    print(f"Total Unique Interacting Pairs: {num_pairs}")

    # 4. Multiprocessing Task Distribution
    num_frames = len(u.trajectory)
    num_workers = min(args.cpus, multiprocessing.cpu_count())
    frame_chunks = np.array_split(np.arange(num_frames), num_workers)
    
    pool_tasks = [
        (pdb_path, xtc_path, cached_nodes, pair_indices, chunk, args.cutoff)
        for chunk in frame_chunks if len(chunk) > 0
    ]

    print(f"🚀 Spawning {len(pool_tasks)} workers across {num_workers} CPUs (COM Cutoff: {args.cutoff} Å)...")
    
    with multiprocessing.Pool(processes=num_workers) as pool:
        chunked_results = pool.map(compute_all_pair_energies_worker, pool_tasks)

    # 5. Consolidate results into a single matrix
    energy_matrix = np.zeros((num_frames, num_pairs), dtype=np.float32)
    for worker_res in chunked_results:
        for f_idx, frame_energies in worker_res:
            energy_matrix[f_idx] = frame_energies

    df_energy = pd.DataFrame(energy_matrix, columns=pair_labels)
    
    # Filter out columns that remain zero across all frames
    non_zero_cols = df_energy.columns[(df_energy != 0).any(axis=0)]
    df_energy = df_energy[non_zero_cols]
    
    df_energy.to_csv(output_csv, index=False)

    print("==================================================")
    print("🎉 ALL-AGAINST-ALL ENERGY EXTRACTION COMPLETE!")
    print(f"Dimensions: {df_energy.shape[0]} Frames x {df_energy.shape[1]} Active Interacting Pairs")
    print(f"Saved to: {output_csv}")
    print("==================================================")

if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
