import os
import sys
import argparse
import pandas as pd
import numpy as np
import networkx as nx
import MDAnalysis as mda
import re

# Force Matplotlib to use a headless backend before importing pyplot
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def analyze_network_topology(run_id, pdb_file, mutations_csv=None):
    data_dir = "../data"
    output_dir = "../results/figures"
    
    os.makedirs(output_dir, exist_ok=True)
    
    edge_list_csv = os.path.join(data_dir, f"{run_id}_edge_list.csv")
    metrics_csv = os.path.join(data_dir, f"{run_id}_centrality_metrics.csv")
    
    if not os.path.exists(edge_list_csv):
        print(f"❌ Error: {edge_list_csv} missing!")
        sys.exit(1)

    # ==========================================================================
    # 🧬 ROBUST MUTATION PARSER & MATCHING ENGINE
    # ==========================================================================
    mutation_resnums = set()
    
    if not mutations_csv:
        mutations_csv = os.path.join(data_dir, "dyrk1a_mutations_list.csv")
        
    if os.path.exists(mutations_csv):
        print(f"📖 Loading disease mutation file: {mutations_csv}")
        
        try:
            df_mut = pd.read_csv(mutations_csv, sep=None, engine='python')
        except Exception as e:
            print(f"⚠️ Standard read failed, falling back to tab-separated: {e}")
            df_mut = pd.read_csv(mutations_csv, sep='\t')
            
        missense_rows = pd.Series(False, index=df_mut.index)
        for col in df_mut.columns:
            missense_rows |= df_mut[col].astype(str).str.lower().str.contains('missense', na=False)
            
        df_filtered = df_mut[missense_rows]
        
        if len(df_filtered) == 0:
            df_filtered = df_mut

        for _, row in df_filtered.iterrows():
            for val in row.values:
                val_str = str(val).strip()
                matches = re.findall(r'\d+', val_str)
                for m in matches:
                    if len(m) >= 2 and int(m) < 2000:
                        mutation_resnums.add(m)
                        
        print(f"🎯 Isolated target mutation resnums: {sorted(list(mutation_resnums), key=int)}")
    else:
        print(f"⚠️ [WARNING] Mutation CSV not found at '{mutations_csv}'. Proceeding without mutation highlights.")

    print(f"=== Running Global Topological Properties on BaNDyT Network: {run_id} ===")
    df_edges = pd.read_csv(edge_list_csv)
    
    G = nx.DiGraph() 
    u = mda.Universe(pdb_file)
    all_protein_residues = []
    
    # 🎯 TRANSLATOR DICTIONARY (Phase 5 pure numbers -> Phase 6 structural labels)
    node_translator = {}
    
    for res in u.residues:
        if res.resname not in ["SOL", "WAT", "HOH", "TIP3", "CL", "NA", "MG"]:
            chain = res.chainID if (hasattr(res, 'chainID') and res.chainID) else "A"
            node_name = f"{chain}:{res.resname}:{res.resnum}"
            all_protein_residues.append(node_name)
            
            # Accommodates both pure numbers ('207') and AA labels ('LEU207') from Phase 4
            node_translator[str(res.resnum)] = node_name
            node_translator[f"{res.resname}{res.resnum}"] = node_name
            node_translator[res.resname] = node_name
            
    G.add_nodes_from(all_protein_residues)
    
    for _, row in df_edges.iterrows():
        raw_source = str(row["Source"])
        raw_target = str(row["Target"])
        
        source_node = node_translator.get(raw_source, raw_source)
        target_node = node_translator.get(raw_target, raw_target)
        
        if source_node in G and target_node in G:
            G.add_edge(source_node, target_node, weight=float(row["Weight"]))
            
    print("Calculating native weighted degree metrics...")
    weighted_degree_dict = dict(G.degree(weight='weight'))

    print("Computing shortest-path centralities (Betweenness)...")
    betweenness_dict = nx.betweenness_centrality(G, weight='weight')
    
    metrics_df = pd.DataFrame({
        "Weighted_Degree": pd.Series(weighted_degree_dict),
        "Betweenness": pd.Series(betweenness_dict)
    }).sort_values(by="Weighted_Degree", ascending=False)
    metrics_df.to_csv(metrics_csv)

    # ==========================================================================
    # 🎨 VISUAL DESIGN CONFIGURATION
    # ==========================================================================
    STRONG_PATHWAY_THRESHOLD = 0.70  
    
    # Build a clean undirected layout graph containing ONLY connected edges
    G_layout = nx.Graph() 
    for u_node, v_node, d in G.edges(data=True):
        G_layout.add_edge(u_node, v_node, weight=d.get("weight", 0.0))
            
    # 🚀 FIX: Drop all isolated/orbit nodes completely. Only keep nodes with at least 1 connection.
    connected_nodes = [node for node, degree in G_layout.degree() if degree >= 1]
    G_layout = G_layout.subgraph(connected_nodes).copy()
    
    # 🚀 FIX: Lock thresholds exactly to the Top 10% (90th percentile)
    hub_threshold = np.percentile(list(weighted_degree_dict.values()), 90) if weighted_degree_dict else 0.0
    betweenness_threshold = np.percentile(list(betweenness_dict.values()), 90) if betweenness_dict else 0.0
    
    # Generate the organic circular physics layout
    pos = nx.spring_layout(G_layout, k=1.2, scale=2.5, iterations=200, seed=42, weight='weight')
        
    print("Applying global dynamic-radius relaxation to prevent node overlap...")
    for _ in range(150):      
        moved = False
        for node1 in G_layout.nodes():
            is_hub1 = weighted_degree_dict.get(node1, 0) >= hub_threshold
            parts1 = node1.split(":")
            resnum1 = parts1[2] if len(parts1) == 3 else ""
            is_mut1 = resnum1 in mutation_resnums
            r1 = 0.40 if (is_mut1 and is_hub1) else (0.25 if is_hub1 else 0.05)
            
            for node2 in G_layout.nodes():
                if node1 == node2:
                    continue
                is_hub2 = weighted_degree_dict.get(node2, 0) >= hub_threshold
                parts2 = node2.split(":")
                resnum2 = parts2[2] if len(parts2) == 3 else ""
                is_mut2 = resnum2 in mutation_resnums
                r2 = 0.40 if (is_mut2 and is_hub2) else (0.25 if is_hub2 else 0.05)
                
                dx = pos[node1][0] - pos[node2][0]
                dy = pos[node1][1] - pos[node2][1]
                dist = np.sqrt(dx**2 + dy**2)
                min_dist = r1 + r2 + 0.06
                
                if dist < min_dist:
                    moved = True
                    if dist == 0:
                        dx, dy = np.random.rand(), np.random.rand()
                        dist = np.sqrt(dx**2 + dy**2)
                    push = (min_dist - dist) * 0.18
                    pos[node1][0] += (dx / dist) * push
                    pos[node1][1] += (dy / dist) * push
                    pos[node2][0] -= (dx / dist) * push
                    pos[node2][1] -= (dy / dist) * push
        if not moved:
            break  
            
    node_colors, node_sizes, node_edge_colors, node_lw = [], [], [], []
    red_mutation_nodes = set()
    
    for node in G_layout.nodes():
        parts = node.split(":")
        resnum = parts[2] if len(parts) == 3 else ""
        is_mutation = resnum in mutation_resnums
        is_hub = weighted_degree_dict.get(node, 0) >= hub_threshold
        is_high_betweenness = betweenness_dict.get(node, 0) >= betweenness_threshold
        
        # Allosteric Bridge Borders (Top 10% Betweenness)
        if is_high_betweenness:
            node_edge_colors.append("#263238") 
            node_lw.append(1.5)
        else:
            node_edge_colors.append("#CFD8DC") 
            node_lw.append(0.4)

        # Node Fill Colors & Sizes
        if is_mutation:
            if is_hub:
                node_colors.append("#D32F2F") # Red Hub
                node_sizes.append(850)         
                red_mutation_nodes.add(node)
            else:
                node_colors.append("#FF9800") # Orange minor mutation
                node_sizes.append(60)          
        else:
            if is_hub:
                node_colors.append("#1565C0") # Blue Hub
                node_sizes.append(600)         
            else:
                node_colors.append("#90A4AE") # Standard node
                node_sizes.append(45)          

    clean_labels = {}
    for n in G_layout.nodes():
        parts = n.split(":")
        resnum = parts[2] if len(parts) == 3 else ""
        is_mutation = resnum in mutation_resnums
        is_hub = weighted_degree_dict.get(n, 0) >= hub_threshold
        
        if is_mutation and is_hub:
            clean_labels[n] = f"{parts[1]}-{parts[2]}" if len(parts) == 3 else n
        else:
            clean_labels[n] = ""

    background_edges = []
    standard_strong_edges = []
    allosteric_highways = [] 
    
    for u_node, v_node, d in G.edges(data=True):
        if u_node in G_layout and v_node in G_layout:
            w = d.get("weight", 0.0)
            if w >= STRONG_PATHWAY_THRESHOLD:
                if u_node in red_mutation_nodes or v_node in red_mutation_nodes:
                    allosteric_highways.append((u_node, v_node))
                else:
                    standard_strong_edges.append((u_node, v_node))
            else:
                background_edges.append((u_node, v_node))

    # ==========================================================================
    # 🎨 COMPACT LEGEND MATRIX 
    # ==========================================================================
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='Critical Hub + Disease Associated Mutation (Red)', markerfacecolor='#D32F2F', markeredgecolor='#D32F2F', markersize=10),
        Line2D([0], [0], marker='o', color='w', label='Critical Hub (Blue)', markerfacecolor='#1565C0', markeredgecolor='#1565C0', markersize=8),
        Line2D([0], [0], marker='o', color='w', label='Disease Associated Mutation (Orange)', markerfacecolor='#FF9800', markeredgecolor='#FF9800', markersize=5),
        Line2D([0], [0], marker='o', color='w', label='Allosteric Bridge Node (Charcoal Border)', markerfacecolor='none', markeredgecolor='#263238', markeredgewidth=1.5, markersize=8),
        
        # 🔗 EDGE / PATHWAY CLASSIFICATIONS
        Line2D([0], [0], lw=1.8, color='#E91E63', label='Allosteric Communication Highway'),
        Line2D([0], [0], lw=1.0, color='#37474F', label='Functional Pathway')
    ]

    def draw_base_graph():
        nx.draw_networkx_edges(G_layout, pos, edgelist=background_edges, width=0.6, edge_color="#90A4AE", alpha=0.75, arrows=True, arrowstyle="->", arrowsize=4, node_size=node_sizes)
        nx.draw_networkx_edges(G_layout, pos, edgelist=standard_strong_edges, width=0.6, edge_color="#37474F", alpha=0.75, arrows=True, arrowstyle="->", arrowsize=6, node_size=node_sizes)
        nx.draw_networkx_edges(G_layout, pos, edgelist=allosteric_highways, width=1.5, edge_color="#E91E63", alpha=0.95, arrows=True, arrowstyle="->", arrowsize=10, node_size=node_sizes)
        nx.draw_networkx_nodes(G_layout, pos, node_size=node_sizes, node_color=node_colors, edgecolors=node_edge_colors, linewidths=node_lw)
        
        plt.legend(
            handles=legend_elements, 
            loc='upper left', 
            fontsize=20,          
            labelspacing=0.35,     
            handletextpad=0.5,     
            handlelength=1.2,      
            borderpad=0.5,         
            frameon=True, 
            facecolor='white', 
            edgecolor='#CFD8DC'
        )

    # --------------------------------------------------------------------------
    # RENDER UNLABELED IMAGE
    # --------------------------------------------------------------------------
    plt.figure(figsize=(16, 16))
    draw_base_graph()
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{run_id}_network_topology_publication.png"), dpi=300)
    plt.close()

    # --------------------------------------------------------------------------
    # RENDER LABELED IMAGE
    # --------------------------------------------------------------------------
    plt.figure(figsize=(16, 16))
    draw_base_graph()
    
    for node, label in clean_labels.items():
        if label:
            x, y = pos[node]
            mag = np.sqrt(x**2 + y**2)
            dx, dy = (x / mag) * 0.22 if mag > 0 else 0.22, (y / mag) * 0.22 if mag > 0 else 0.22
            
            plt.text(
                x + dx, y + dy, label, 
                fontsize=16, 
                fontweight="bold", 
                color="#D32F2F", 
                ha='center', 
                va='center',
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#D32F2F", alpha=0.95, lw=1.6)
            )

    plt.axis("off")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{run_id}_network_topology_labeled.png"), dpi=300)
    plt.close()

    print(f"📸 Network topology plots rendered successfully to: {output_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True)
    parser.add_argument("--pdb", required=True)
    parser.add_argument("--mutations", required=False, default=None)
    args = parser.parse_args()
    analyze_network_topology(args.run_id, args.pdb, args.mutations)
