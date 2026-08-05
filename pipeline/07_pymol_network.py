import os
import sys
import argparse
import pandas as pd
import numpy as np

def generate_pymol_network(run_id, pdb_filename, threshold):
    data_dir = "../data"
    input_dir = "../input"
    output_dir = "../results/figures"
    
    os.makedirs(output_dir, exist_ok=True)
    
    edge_list_csv = os.path.join(data_dir, f"{run_id}_edge_list.csv")
    pdb_path = os.path.join(input_dir, pdb_filename)
    pdb_abs_path = os.path.abspath(pdb_path)  # Get absolute path for robust loading
    pml_output = os.path.join(output_dir, f"{run_id}_network_3d.pml")
    
    if not os.path.exists(edge_list_csv):
        print(f"❌ Error: Edge list file '{edge_list_csv}' not found!")
        sys.exit(1)
        
    print(f"=== Generating 3D PyMOL Network Visualization Script ===")
    print(f"System ID: {run_id}")
    print(f"PDB Source: {pdb_abs_path}")
    print(f"Correlation Threshold: >= {threshold}")
    
    # Load edges
    df = pd.read_csv(edge_list_csv)
    
    # Filter edges above absolute threshold
    df_filtered = df[df['Weight'].abs() >= threshold].copy()
    print(f"Found {len(df_filtered)} strong functional pathways above threshold.")
    
    if len(df_filtered) == 0:
        print("❌ Error: No edges found above the specified threshold!")
        print("Please lower the --threshold parameter and run again.")
        sys.exit(1)
        
    # Calculate Node Strength (Weighted Degree) based on filtered network
    node_strengths = {}
    for _, row in df_filtered.iterrows():
        src, tgt, w = row['Source'], row['Target'], abs(row['Weight'])
        node_strengths[src] = node_strengths.get(src, 0.0) + w
        node_strengths[tgt] = node_strengths.get(tgt, 0.0) + w
        
    max_strength = max(node_strengths.values()) if node_strengths else 1.0
    
    # Build the PML script contents
    pml_lines = []
    pml_lines.append("# ==========================================================================")
    pml_lines.append(f"# PyMOL 3D Network Visualization Script for {run_id}")
    pml_lines.append(f"# Generated automatically by 07_pymol_network.py")
    pml_lines.append("# ==========================================================================")
    pml_lines.append("reinitialize")
    pml_lines.append("bg_color white")
    pml_lines.append("")
    pml_lines.append("# Load and style the protein structure")
    pml_lines.append(f"load {pdb_abs_path}, system")
    pml_lines.append("show cartoon, system")
    pml_lines.append("color grey90, system")
    pml_lines.append("set cartoon_transparency, 0.0")  # Keep fully opaque ribbon
    pml_lines.append("hide lines, system")
    pml_lines.append("hide nonbonded, system")
    pml_lines.append("")
    pml_lines.append("# ==========================================================================")
    pml_lines.append("# 🔴 RENDER NODES (Sleek Micro-Beads scaled by Weighted Degree)")
    pml_lines.append("# ==========================================================================")
    
    rendered_nodes = set()
    
    def parse_node(node_str):
        parts = node_str.split(":")
        if len(parts) == 3:
            return parts[0], parts[1], parts[2]
        elif len(parts) == 2:
            return "A", parts[0], parts[1]
        else:
            return "A", "UNK", node_str

    for node, strength in node_strengths.items():
        chain, resname, resnum = parse_node(node)
        node_alias = f"node_{resname}_{resnum}"
        
        # Sleek sphere scale (0.25 to 0.75) prevents massive sphere collisions
        sphere_scale = 0.25 + (strength / max_strength) * 0.50
        
        pml_lines.append(f"select {node_alias}, system and resi {resnum} and name CA")
        pml_lines.append(f"show spheres, {node_alias}")
        pml_lines.append(f"set sphere_scale, {sphere_scale:.3f}, {node_alias}")
        pml_lines.append(f"color red, {node_alias}")
        rendered_nodes.add(node_alias)
        
    pml_lines.append("")
    pml_lines.append("# ==========================================================================")
    pml_lines.append("# 🔴 RENDER EDGES (Ultra-Thin Wireframe Cylinders)")
    pml_lines.append("# ==========================================================================")
    
    for idx, row in df_filtered.iterrows():
        src, tgt, w = row['Source'], row['Target'], abs(row['Weight'])
        
        src_chain, src_resname, src_resnum = parse_node(src)
        tgt_chain, tgt_resname, tgt_resnum = parse_node(tgt)
        
        edge_alias = f"edge_{src_resnum}_{tgt_resnum}"
        
        # Super-thin wireframe scaling (0.01 to 0.04)
        if 1.0 - threshold > 0:
            norm_w = (w - threshold) / (1.0 - threshold)
        else:
            norm_w = 1.0
        cylinder_radius = 0.01 + norm_w * 0.03
        
        src_sel = f"(system and resi {src_resnum} and name CA)"
        tgt_sel = f"(system and resi {tgt_resnum} and name CA)"
        
        pml_lines.append(f"distance {edge_alias}, {src_sel}, {tgt_sel}")
        pml_lines.append(f"set dash_gap, 0, {edge_alias}")
        pml_lines.append(f"set dash_radius, {cylinder_radius:.3f}, {edge_alias}")
        pml_lines.append(f"color red, {edge_alias}")
        
    pml_lines.append("")
    pml_lines.append("# Clean up layout and display settings")
    pml_lines.append("hide labels, edge_*")
    pml_lines.append("set dash_round_ends, 1")
    pml_lines.append("deselect")
    pml_lines.append("reset")
    pml_lines.append("zoom system")
    pml_lines.append("")
    pml_lines.append("python")
    pml_lines.append("print('\\n🎉 3D Allosteric Coupling Network Loaded Successfully!\\n')")
    pml_lines.append("python end")
    
    # Save the file
    with open(pml_output, "w") as f:
        f.write("\n".join(pml_lines))
        
    print(f"🎉 SUCCESS: PyMOL Network script saved to: {pml_output}")
    print("To open this in PyMOL, run:")
    print(f"  pymol {pml_output}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="3D Kinase Allosteric Network Mapper for PyMOL")
    parser.add_argument("--run_id", required=True, help="Prefix name of the system tracking folder")
    parser.add_argument("--pdb", required=True, help="Filename of the PDB topology inside ../input")
    parser.add_argument("--threshold", type=float, default=0.78, help="Absolute correlation threshold to display as 3D cylinders")
    
    args = parser.parse_args()
    generate_pymol_network(args.run_id, args.pdb, args.threshold)
