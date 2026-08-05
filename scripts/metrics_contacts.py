import sys
import re

def parse_contacts(file_path):
    unique_atom_pairs = set()
    unique_residue_pairs = set()

    # Regex to recognize a standard chain:resname:resnum:atom identifier (e.g., B:GLU:291:OE1)
    # This also safely captures multi-character chain IDs or non-standard residues (e.g., PTR:321)
    contact_pattern = re.compile(r'([A-Za-z0-9]+:[A-Za-z0-9]+:[0-9\-]+):[A-Za-z0-9\-\'\*]+')

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip empty lines and comments/headers
            if not line or line.startswith('#'):
                continue
            
            # Split line into tokens
            tokens = line.split()
            
            # Find all tokens matching the atom identifier pattern
            matches = []
            res_matches = []
            for t in tokens:
                match = contact_pattern.search(t)
                if match:
                    matches.append(match.group(0))      # Full atom specifier
                    res_matches.append(match.group(1))  # Residue part only (chain:resname:resnum)
            
            # Ensure we found exactly one pair of interacting partners
            if len(matches) >= 2:
                atom1, atom2 = matches[0], matches[1]
                res1, res2 = res_matches[0], res_matches[1]
                
                # Sort pairs alphabetically to avoid duplicate counting of directionality (A->B vs B->A)
                atom_pair = tuple(sorted([atom1, atom2]))
                res_pair = tuple(sorted([res1, res2]))
                
                unique_atom_pairs.add(atom_pair)
                unique_residue_pairs.add(res_pair)

    print(f"--- Analysis Summary for {file_path} ---")
    print(f"Total Unique Atom-to-Atom Contacts:    {len(unique_atom_pairs)}")
    print(f"Total Unique Residue-to-Residue Contacts: {len(unique_residue_pairs)}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <path_to_getcontacts_output>")
        sys.exit(1)
        
    parse_contacts(sys.argv[1])

