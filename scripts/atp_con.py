import sys
import os

# ==============================================================================
# CONFIGURATION
# ==============================================================================
CONTACTS_TSV = "../data/ogxtc/dyrk1a_ATP_run_contacts.tsv" 
OUTPUT_TXT = "pymol_magenta_commands.txt"
LIGAND_NAME = "ATP"

# ==============================================================================
# PARSING ENGINE
# ==============================================================================
contact_residues = set()

if not os.path.exists(CONTACTS_TSV):
    print(f"❌ Error: Contacts file '{CONTACTS_TSV}' not found!")
    sys.exit(1)

print(f"📖 Parsing contacts from: {CONTACTS_TSV}")

with open(CONTACTS_TSV, "r") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        
        tokens = line.strip().split()
        for token in tokens:
            if ":" in token:
                parts = token.split(":")
                if len(parts) >= 3:
                    chain = parts[0]
                    resname = parts[1]
                    resnum = parts[2]
                    
                    if resname != LIGAND_NAME:
                        contact_residues.add((chain, resnum, resname))

if not contact_residues:
    print("⚠️ Warning: No contact residues were found!")
    sys.exit(0)

print(f"🎯 Successfully extracted {len(contact_residues)} unique contact residues.")

chain_groups = {}
sorted_contacts = sorted(contact_residues, key=lambda x: int(x[1]) if x[1].isdigit() else x[1])

for chain, resi, resn in sorted_contacts:
    if chain not in chain_groups:
        chain_groups[chain] = []
    chain_groups[chain].append(resi)

# ==============================================================================
# WRITE PYMOL COMMANDS
# ==============================================================================
with open(OUTPUT_TXT, "w") as out:
    out.write("# ==================================================\n")
    out.write("# PYMOL COMMANDS TO ISOLATE & COLOR CONTACT RESIDUES\n")
    out.write(f"# Compiled from: {CONTACTS_TSV}\n")
    out.write("# ==================================================\n\n")
    
    selection_names = []
    
    for chain, resis in chain_groups.items():
        unique_resis = sorted(list(set(resis)), key=lambda x: int(x) if x.isdigit() else x)
        resi_str = "+".join(unique_resis)
        sel_name = f"contacts_ch{chain}"
        
        out.write(f"select {sel_name}, (chain {chain} and resi {resi_str})\n")
        selection_names.append(sel_name)
    
    out.write("\n# Combine selection groups, apply magenta coloring, and show as sticks\n")
    all_selections = " or ".join(selection_names)
    out.write(f"select active_site_contacts, {all_selections}\n")
    out.write("color magenta, active_site_contacts\n")
    out.write("show sticks, active_site_contacts\n")
    out.write("deselect\n")

print(f"🎉 SUCCESS: PyMOL command script written to '{OUTPUT_TXT}'!")
