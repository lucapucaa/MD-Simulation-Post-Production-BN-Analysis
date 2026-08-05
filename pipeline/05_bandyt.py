import os
import sys
import pickle
import argparse
import pandas as pd
import numpy as np
import networkx as nx

def run_bandyt_analysis(run_id, restarts, bandyt_path="/home/lchill/work/dyrk1a-project/Simulation_data/bandyt/bandyt/"):
    # Anchor paths relative to pipeline execution folder
    data_dir = "../data"
    output_dir = "../results/figures"
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    
    # Handle custom package paths if specified
    if bandyt_path:
        bandyt_path = os.path.abspath(bandyt_path.rstrip('/'))
        print(f"Appending custom BaNDyT search path: {bandyt_path}")
        
        # Inject the direct folder AND its parent folder to cover all import styles
        if bandyt_path not in sys.path:
            sys.path.insert(0, bandyt_path)
        parent_dir = os.path.dirname(bandyt_path)
        if parent_dir not in sys.path:
            sys.path.insert(0, parent_dir)
        
    try:
        # Direct clean import
        import bandyt
        
        # If it loaded the outer folder namespace, step inside to the true module
        if hasattr(bandyt, 'bandyt'):
            print("ℹ️ Resolving nested package mapping...")
            bandyt = bandyt.bandyt
            
        # Final validation check for the target function
        if not hasattr(bandyt, 'read_input_file'):
            print("❌ Error: Loaded 'bandyt' module is missing core functions.")
            sys.exit(1)
            
    except ImportError as e:
        print(f"❌ Critical Import Error: {e}")
        print("Python tried searching in these directories:")
        for path in sys.path[:3]:
            print(f"  -> {path}")
        sys.exit(1)
        
    # Check for inputs matching standard Phase 4 outputs
    input_csv = os.path.join(data_dir, f"{run_id}_interaction_energies.csv")
    if not os.path.exists(input_csv):
        # Fallback to template feature label naming
        input_csv = os.path.join(data_dir, f"{run_id}_input_features.csv")
        
    model_pickle = os.path.join(data_dir, f"{run_id}_model.pkl")
    pickle_output = os.path.join(data_dir, f"{run_id}.pickle")
    properties_output = os.path.join(data_dir, f"{run_id}-properties")
    edge_list_csv = os.path.join(data_dir, f"{run_id}_edge_list.csv")
    
    if not os.path.exists(input_csv):
        print(f"❌ Error: Input dataset matrix '{input_csv}' missing!")
        sys.exit(1)
        
    print(f"=== Running BaNDyT Bayesian Network Structure Learning on: {run_id} ===")
    df = pd.read_csv(input_csv)
    
    # 📉 🚀 1000ns FIX: Set Stride to 100
    original_row_count = df.shape[0]
    df = df.iloc[::100].copy()
    print(f"📉 1000ns Trajectory downsampled using a stride of 100: {original_row_count} frames -> {df.shape[0]} frames.")

    # Drop the chronological frame index so it isn't treated as a physical variable
    if 'Frame' in df.columns:
        df = df.drop(columns=['Frame'])

    # 🧬 🚀 TOPOLOGY FIX: Convert Pairwise Contacts into Individual Residue Nodes
    print("🔄 Collapsing pairwise interaction energies into per-residue dynamic profiles...")
    original_pairs = df.shape[1]
    residue_df = pd.DataFrame(index=df.index)
    
    for col in df.columns:
        # Look for the dash that separates the two residues in your Phase 4 output
        if '-' in col:
            res1, res2 = col.split('-', 1)
            
            # Accumulate energy for Residue 1
            if res1 in residue_df:
                residue_df[res1] += df[col]
            else:
                residue_df[res1] = df[col].copy()
                
            # Accumulate energy for Residue 2
            if res2 in residue_df:
                residue_df[res2] += df[col]
            else:
                residue_df[res2] = df[col].copy()
        else:
            residue_df[col] = df[col]
            
    df = residue_df.copy()
    print(f"✅ Network restructured: {original_pairs} pairs collapsed into {df.shape[1]} unique individual residues.")
        
    # 🧹 RAM SAVER: Drop residues that are completely static (Variance < 0.2)
    print("🧹 Filtering out rigidly static residues to isolate functional pathways...")
    initial_cols = df.shape[1]
    variances = df.var(numeric_only=True)
    active_cols = variances[variances > 0.2].index
    df = df[active_cols].copy()
    
    print(f"Ingested dataset with {df.shape[1]} active residues (Dropped {initial_cols - df.shape[1]} static residues).")
    print(f"📊 Final Matrix Size going into BaNDyT: {df.shape[0]} frames x {df.shape[1]} residues")
    
    # Serialize temporary downsampled CSV because BaNDyT ingests data via paths
    temp_csv = os.path.join(data_dir, f"temp_downsampled_{run_id}.csv")
    df.to_csv(temp_csv, index=False)
    
    # Load input data and automatically discretize continuous energies into 8 bins
    print("Discretizing continuous energy spaces via maximum entropy algorithm...")
    dt = bandyt.read_input_file(temp_csv)
    
    # Initialize Bayesian network search loop and auto-detect C-scoring capabilities
    print("Initializing Bayesian topology architecture...")
    try:
        srch = bandyt.search(dt, ofunc=bandyt.cmu)
        print("⚡ Success: Compiled C-scoring acceleration (cmu) activated!")
    except (AttributeError, NameError):
        srch = bandyt.search(dt)
        print("ℹ️ Using default standard Python scoring engine (mu).")
        
    # Perform the recursive search of optimal network topologies
    print(f"Executing structure search across {restarts} randomized restarts for network convergence...")
    srch.restarts(nrestarts=restarts)
    
    # Force BaNDyT to save the file inside data_dir
    print("Exporting raw Graphviz network structure topology...")
    original_cwd = os.getcwd()
    try:
        os.chdir(os.path.abspath(data_dir))
        srch.dot(path=run_id)
    finally:
        os.chdir(original_cwd)
        
    # Dynamic identification of the Dot artifact path
    dot_file_path = None
    if os.path.exists(os.path.join(data_dir, f"{run_id}.dot")):
        dot_file_path = os.path.join(data_dir, f"{run_id}.dot")
    elif os.path.exists(os.path.join(data_dir, run_id)):
        dot_file_path = os.path.join(data_dir, run_id)
        
    if dot_file_path:
        print(f"💾 Verified dot file asset path: {dot_file_path}")
        print("Compiling network topology into a high-resolution vector PDF format...")
        try:
            import pydot
            pdf_output_path = os.path.join(output_dir, f"{run_id}_bayesian_network.pdf")
            (graph,) = pydot.graph_from_dot_file(dot_file_path)
            graph.write_pdf(pdf_output_path)
            print(f"📸 SUCCESS: Vector network map rendered to: {pdf_output_path}")
        except Exception as e:
            print(f"⚠️ Automated PDF rendering skipped: {e}")
    else:
        print("⚠️ Warning: Could not locate exported dot file structure.")
        
    # Convert Bayesian network to an igraph instance
    print(f"Converting architecture to igraph representation: {pickle_output}")
    bandyt.convert_bn_to_igraph(srch, fout=pickle_output, format="pickle")
    
    # Generate weighted degrees and graph property logs
    print(f"Exporting global network metrics: {properties_output}")
    bandyt.getGraphProp(pickle_output, properties_output)
    
    # Clean up temporary file artifact
    if os.path.exists(temp_csv):
        os.remove(temp_csv)
        
    # Translate igraph metadata properties into your downstream edge lists
    print("Translating graph topologies into clean dataframes...")
    try:
        with open(pickle_output, "rb") as f:
            g = pickle.load(f)
            
        # Version-agnostic attribute lookup mapping
        vertex_attrs = g.vs.attribute_names() if hasattr(g.vs, 'attribute_names') else []
        edge_attrs = g.es.attribute_names() if hasattr(g.es, 'attribute_names') else []
        
        name_attr = "name" if "name" in vertex_attrs else ("label" if "label" in vertex_attrs else None)
        
        edges_data = []
        for edge in g.es:
            source_node = g.vs[edge.source][name_attr] if name_attr else str(edge.source)
            target_node = g.vs[edge.target][name_attr] if name_attr else str(edge.target)
            
            weight = 1.0
            for attr in ["weight", "influence", "probability", "score"]:
                if attr in edge_attrs:
                    weight = edge[attr]
                    break
            edges_data.append({"Source": source_node, "Target": target_node, "Weight": weight})
            
        edges_df = pd.DataFrame(edges_data)
    except Exception as e:
        print(f"⚠️ Direct pickle ingestion skipped ({e}). Extracting from edge property logs...")
        prop_edges_csv = f"{properties_output}_edges.csv"
        if os.path.exists(prop_edges_csv):
            edges_df = pd.read_csv(prop_edges_csv)
            edges_df.rename(columns={"source": "Source", "target": "Target", "weight": "Weight"}, inplace=True)
        else:
            edges_df = pd.DataFrame(columns=["Source", "Target", "Weight"])
            
    # Clean and save the edge lists without an arbitrary weight threshold
    if not edges_df.empty:
        # Eliminate structural self-interactions loops (Source == Target)
        edges_df = edges_df[edges_df["Source"] != edges_df["Target"]].copy()
        
        # Save complete parsed continuous edge list for visual scripts (Phase 6 & 7)
        edges_df.to_csv(edge_list_csv, index=False)
        print(f"📊 Flat edge data list written to: {edge_list_csv}")
        
        # Build and pickle a NetworkX graph structure for backward compatibility
        G_nx = nx.from_pandas_edgelist(
            edges_df, 
            source="Source", 
            target="Target", 
            edge_attr="Weight", 
            create_using=nx.DiGraph()
        )
        with open(model_pickle, "wb") as f:
            pickle.dump(G_nx, f)
            
        print(f"🎉 SUCCESS: Causal network analysis complete for {run_id} ({len(edges_df)} relationships saved).")
    else:
        print("❌ Error: Network parsing returned an empty edge array structure.")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Integrated BaNDyT Structure Logic and Downsampler")
    parser.add_argument("--run_id", required=True, help="Prefix name of the system tracking folder")
    parser.add_argument("--restarts", type=int, default=50, help="Number of search restarts for network convergence")
    
    args = parser.parse_args()
    run_bandyt_analysis(args.run_id, args.restarts)
