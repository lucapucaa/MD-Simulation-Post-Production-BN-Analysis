import os
import sys
import argparse
import re
import pandas as pd

def parse_node(node_str):
    """
    Robustly parses any node identifier (e.g., 'ASP_304', 'A:ASP:304', 'ASP304', '304')
    into a clean tuple of (chain, resname, resnum) where resnum is a string of digits.
    """
    node_str = str(node_str).strip()
    
    # Normalize common delimiters to spaces and split
    cleaned = re.sub(r'[-_:]', ' ', node_str).split()
    
    if len(cleaned) == 3:
        # Match pattern like: A ASP 304
        nums = [x for x in cleaned if x.isdigit()]
        letters = [x for x in cleaned if x.isalpha()]
        resnum = nums[0] if nums else "UNK"
        resname = [x for x in letters if len(x) == 3][0] if [x for x in letters if len(x) == 3] else "UNK"
        chain = [x for x in letters if len(x) == 1][0] if [x for x in letters if len(x) == 1] else "A"
        return chain, resname, resnum
    elif len(cleaned) == 2:
        # Match pattern like: ASP 304
        nums = [x for x in cleaned if x.isdigit()]
        letters = [x for x in cleaned if x.isalpha()]
        resnum = nums[0] if nums else "UNK"
        resname = letters[0] if letters else "UNK"
        return "A", resname, resnum
        
    # Regex fallback if the parts are smashed together (e.g. 'ASP304')
    resnum_match = re.search(r'\d+', node_str)
    resnum = resnum_match.group(0) if resnum_match else "UNK"
    
    resname_match = re.search(r'[A-Za-z]{3}', node_str)
    resname = resname_match.group(0) if resname_match else "UNK"
    
    chain_match = re.search(r'^[A-Za-z](?=:|_)', node_str)
    chain = chain_match.group(0) if chain_match else "A"
    
    return chain, resname, resnum

def generate_active_site_network(run_id, pdb_filename, threshold, ligand_name):
    data_dir = "../data"
    input_dir = "../input"
    output_dir = "../results/figures"
    
    os.makedirs(output_dir, exist_ok=True)
    
    contacts_tsv = os.path.join(data_dir, f"{run_id}_contacts.tsv")
    edge_list_csv = os.path.join(data_dir, f"{run_id}_edge_list.csv")
    pdb_path = os.path.join(input_dir, pdb_filename)
    pdb_abs_path = os.path.abspath(pdb_path)
    pml_output = os.path.join(output_dir, f"{run_id}_active_site_network.pml")
    
    # --- STEP 1: PARSE GETCONTACTS TO ISOLATE ACTIVE SITE RESIDUE NUMBERS ---
    if not os.path.exists(contacts_tsv):
        print(f"❌ Error: Contacts file '{contacts_tsv}' not found!")
        sys.exit(1)
        
    print(f"🔍 Scanning {contacts_tsv} for direct {ligand_name} interactions...")
    active_site_resnums = set()
    
    with open(contacts_tsv, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            tokens = line.strip().split()
            if len(tokens) >= 4:
                res1, res2 = tokens[2], tokens[3]
                # If either partner matches the target ligand name
                if ligand_name.upper() in res1.upper() or ligand_name.upper() in res2.upper():
                    for token in [res1, res2]:
                        if ligand_name.upper() not in token.upper():
                            # Extract residue number to handle all naming formats
                            match = re.search(r'(?<=:)\d+(?=:|$)|(?<=[A-Za-z])\d+', token)
                            if match:
                                active_site_resnums.add(match.group(0))
                            else:
                                for sub in token.split(':'):
                                    if sub.isdigit():
                                        active_site_resnums.add(sub)
                                        
    print(f"🎯 Isolated {len(active_site_resnums)} active site residues interacting with {ligand_name}: {sorted(list(active_site_resnums))}")
    
    if not active_site_resnums:
        print(f"⚠️ Warning: No direct contacts found for ligand pattern '{ligand_name}'. Check naming format.")
        sys.exit(1)

    # --- STEP 2: INGEST AND FILTER CAUSAL EDGES ---
    if not os.path.exists(edge_list_csv):
        print(f"❌ Error: Edge list file '{edge_list_csv}' not found!")
        sys.exit(1)
        
    df = pd.read_csv(edge_list_csv)
    df_filtered = df[df['Weight'].abs() >= threshold].copy()
    
    # Filter for pathways containing at least one active site residue number
    active_edges = []
    for idx, row in df_filtered.iterrows():
        src_chain, src_name, src_num = parse_node(row['Source'])
        tgt_chain, tgt_name, tgt_num = parse_node(row['Target'])
        
        if src_num in active_site_resnums or tgt_num in active_site_resnums:
            active_edges.append(row)
            
    df_active = pd.DataFrame(active_edges)
    print(f"📈 Filtered network down to {len(df_active)} allosteric pathways anchoring the active site pocket.")
    
    if df_active.empty:
        print("❌ Error: No network connections hit the active site residues at this threshold.")
        sys.exit(1)
        
    # --- STEP 3: MAP VISUAL COLOR NODES GROUPING ---
    unique_nodes = set(df_active['Source'].unique()) | set(df_active['Target'].unique())
    active_nodes = set()
    connecting_nodes = set()
    
    for node in unique_nodes:
        _, _, resnum = parse_node(node)
        if resnum in active_site_resnums:
            active_nodes.add(node)
        else:
            connecting_nodes.add(node)
            
    # --- STEP 4: GENERATE THE PYMOL SCRIPT (.PML) ---
    pml = []
    pml.append("# ==========================================================================")
    pml.append(f"# PyMOL Active Site Network Visualization for {run_id}")
    pml.append("# ==========================================================================")
    pml.append("reinitialize")
    pml.append("bg_color white")
    pml.append(f"load {pdb_abs_path}, system")
    pml.append("show cartoon, system")
    pml.append("color grey90, system")
    pml.append("set cartoon_transparency, 0.3")  # Transparent scaffolding highlight
    pml.append("hide lines, system")
    pml.append("hide nonbonded, system")
    pml.append("")
    
    # Show the ligand itself for positional reference
    pml.append(f"select ligand, system and resn {ligand_name}")
    pml.append("show sticks, ligand")
    pml.append("color magenta, ligand")
    pml.append("")
    
    pml.append("# 🔵 CRITICAL VISUAL RULE: MARK ACTIVE SITE RESIDUES FIRST (BLUE)")
    for node in active_nodes:
        chain, name, num = parse_node(node)
        alias = f"act_{name}_{num}"
        pml.append(f"select {alias}, system and resi {num} and name CA")
        pml.append(f"show spheres, {alias}")
        pml.append(f"set sphere_scale, 0.50, {alias}")
        pml.append(f"color marine, {alias}")  # Clean high-contrast blue
        
    pml.append("")
    pml.append("# 🔴 CRITICAL VISUAL RULE: MARK CONNECTING HOOK RESIDUES SECOND (RED)")
    for node in connecting_nodes:
        chain, name, num = parse_node(node)
        alias = f"conn_{name}_{num}"
        pml.append(f"select {alias}, system and resi {num} and name CA")
        pml.append(f"show spheres, {alias}")
        pml.append(f"set sphere_scale, 0.40, {alias}")
        pml.append(f"color red, {alias}")
        
    pml.append("")
    pml.append("# 🔴 CRITICAL VISUAL RULE: RENDER ALL ALLOSTERIC CONNECTIONS (RED LINES)")
    for idx, row in df_active.iterrows():
        src_chain, src_name, src_num = parse_node(row['Source'])
        tgt_chain, tgt_name, tgt_num = parse_node(row['Target'])
        edge_alias = f"act_edge_{src_num}_{tgt_num}"
        
        src_sel = f"(system and resi {src_num} and name CA)"
        tgt_sel = f"(system and resi {tgt_num} and name CA)"
        
        pml.append(f"distance {edge_alias}, {src_sel}, {tgt_sel}")
        pml.append(f"set dash_gap, 0, {edge_alias}")
        pml.append(f"set dash_radius, 0.03, {edge_alias}")
        pml.append(f"color red, {edge_alias}")
        
    pml.append("")
    pml.append("hide labels, act_edge_*")
    pml.append("set dash_round_ends, 1")
    pml.append("deselect")
    pml.append("reset")
    pml.append("zoom unique_nodes")
    
    with open(pml_output, "w") as f:
        f.write("\n".join(pml))
        
    print(f"\n🎉 SUCCESS: Active site pocket visual script generated at: {pml_output}")
    print(f"👉 Run 'pymol {pml_output}' on your workstation to inspect allosteric pocket loops.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Active Site Network Isolation Tool")
    parser.add_argument("--run_id", required=True)
    parser.add_argument("--pdb", required=True)
    parser.add_argument("--threshold", type=float, default=0.78)
    parser.add_argument("--ligand", default="ATP")
    
    args = parser.parse_args()
    generate_active_site_network(args.run_id, args.pdb, args.threshold, args.ligand)
