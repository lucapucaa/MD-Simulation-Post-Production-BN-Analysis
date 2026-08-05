import os
import re
import sys
import pandas as pd

def parse_missense_dams(mutations_csv):
    """
    Parses the mutation CSV:
    - Column 8 (index 7): Filters for 'missense'
    - Column 9 (index 8): Extracts residue numbers
    """
    mutation_resnums = set()
    if not os.path.exists(mutations_csv):
        print(f"❌ Error: Mutation CSV missing at {mutations_csv}")
        sys.exit(1)

    print(f"📖 Reading mutation file: {mutations_csv}")
    try:
        df_mut = pd.read_csv(mutations_csv, sep='\t')
        if df_mut.shape[1] < 9:
            df_mut = pd.read_csv(mutations_csv) # Fallback to comma delimiter
    except Exception as e:
        print(f"❌ Error reading {mutations_csv}: {e}")
        sys.exit(1)

    try:
        missense_mask = df_mut.iloc[:, 7].astype(str).str.lower().str.contains('missense', na=False)
        df_filtered = df_mut[missense_mask]
        
        raw_resnums = df_filtered.iloc[:, 8].dropna().astype(str).str.strip().tolist()
        for val in raw_resnums:
            match = re.search(r'\d+', val)
            if match:
                mutation_resnums.add(match.group())
                
        print(f"🎯 Isolated {len(mutation_resnums)} unique missense DAM residue numbers.")
        return mutation_resnums
    except IndexError:
        print(f"❌ Error: Mutation CSV must contain at least 9 columns!")
        sys.exit(1)


def parse_contacts_tsv(filepath):
    """
    Extracts protein residue numbers from getContact log TSVs.
    Filters out ATP, solvent, and non-protein residues.
    """
    contact_resnums = set()
    if not os.path.exists(filepath):
        print(f"⚠️ Warning: Contact file not found at {filepath}")
        return contact_resnums

    contact_pattern = re.compile(r'(?:([A-Za-z0-9]+):)?([A-Z0-9]{3,4}):([0-9\-]+)')
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            for token in line.strip().split():
                match = contact_pattern.search(token)
                if match:
                    _, resname, resnum = match.groups()
                    if resname not in ["ATP", "SOL", "WAT", "HOH"] and resnum.isdigit() and int(resnum) >= 100:
                        contact_resnums.add(resnum)
    return contact_resnums


def main():
    data_dir = "../data"
    output_dir = "../results/pymol"
    os.makedirs(output_dir, exist_ok=True)

    mutations_csv = os.path.join(data_dir, "dyrk1a_mutations_list.csv")
    
    # Locate Substrate contact file
    sub_tsv_candidates = ["2w06_sub_contacts.tsv", "2wo6_sub_contacts.tsv"]
    sub_tsv_path = next((os.path.join(data_dir, f) for f in sub_tsv_candidates if os.path.exists(os.path.join(data_dir, f))), os.path.join(data_dir, "2w06_sub_contacts.tsv"))
    
    # Locate ATP contact file
    atp_tsv_candidates = ["dyrk1a_ATP_run_contacts.tsv", "2w06_ATP_contacts.tsv", "2wo6_ATP_contacts.tsv"]
    atp_tsv_path = next((os.path.join(data_dir, f) for f in atp_tsv_candidates if os.path.exists(os.path.join(data_dir, f))), os.path.join(data_dir, "dyrk1a_ATP_run_contacts.tsv"))

    # 1. Parse Missense DAMs
    all_dams = parse_missense_dams(mutations_csv)

    # 2. Parse Contacts
    print(f"📖 Parsing substrate contacts from: {sub_tsv_path}")
    sub_contacts = parse_contacts_tsv(sub_tsv_path)
    
    print(f"📖 Parsing ATP contacts from: {atp_tsv_path}")
    atp_contacts = parse_contacts_tsv(atp_tsv_path)

    # 3. Categorize DAMs according to overlapping pocket logic
    both_dams = all_dams & sub_contacts & atp_contacts
    substrate_only_dams = (all_dams & sub_contacts) - both_dams
    atp_only_dams = (all_dams & atp_contacts) - both_dams
    non_active_dams = all_dams - (sub_contacts | atp_contacts)

    print("\n📊 Mutation Classification Summary:")
    print(f"  🔴 Non-Active Site DAMs: {len(non_active_dams)}")
    print(f"  🔵 Substrate Only DAMs:   {len(substrate_only_dams)}")
    print(f"  🟡 ATP Only DAMs:         {len(atp_only_dams)}")
    print(f"  🟢 Both ATP & Substrate: {len(both_dams)}")

    # 4. Generate Master Combined PML
    master_pml_path = os.path.join(output_dir, "all_categorized_dams.pml")
    
    with open(master_pml_path, "w") as pml:
        pml.write("# ======================================================================\n")
        pml.write("# PYMOL 3D MAP: CATEGORIZED DISEASE-ASSOCIATED MUTATIONS (DAMS)\n")
        pml.write("# ======================================================================\n\n")
        
        pml.write("bg_color white                  # Canvas background\n")
        pml.write("hide everything\n")
        pml.write("show cartoon, all               # Background cartoon backbone\n")
        pml.write("color gray85, all\n\n")

        pml.write("# 🎨 PALETTE DEFINITIONS\n")
        pml.write("set_color non_active_red, [0.827, 0.184, 0.184]   # Non-active site\n")
        pml.write("set_color substrate_blue, [0.082, 0.396, 0.753]   # Substrate only\n")
        pml.write("set_color atp_yellow,     [0.984, 0.753, 0.176]   # ATP only\n")
        pml.write("set_color both_green,     [0.220, 0.557, 0.235]   # Both ATP and Substrate\n\n")

        # helper function to write selection block
        def write_selection(name, res_set, color, label_text):
            res_str = "+".join(sorted(list(res_set), key=int)) if res_set else "none"
            pml.write(f"# {label_text}\n")
            if res_str != "none":
                pml.write(f"select {name}, polymer.protein and resi {res_str} and name CA\n")
                pml.write(f"show spheres, {name}\n")
                pml.write(f"color {color}, {name}\n")
                pml.write(f"set sphere_scale, 0.70, {name}\n\n")
            else:
                pml.write(f"# No residues found for {name}\n\n")

        write_selection("dams_non_active", non_active_dams, "non_active_red", "🔴 NON-ACTIVE SITE DAMS (RED)")
        write_selection("dams_substrate_only", substrate_only_dams, "substrate_blue", "🔵 SUBSTRATE ONLY DAMS (BLUE)")
        write_selection("dams_atp_only", atp_only_dams, "atp_yellow", "🟡 ATP ONLY DAMS (YELLOW)")
        write_selection("dams_both_active", both_dams, "both_green", "🟢 BOTH ATP & SUBSTRATE DAMS (GREEN)")

        pml.write("deselect\n")
        pml.write("set ray_shadows, 0\n")

    print(f"\n✨ Successfully exported master PyMOL script: {master_pml_path}")

if __name__ == "__main__":
    main()
