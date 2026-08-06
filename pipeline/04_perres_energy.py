from energyextract import Simul
import MDAnalysis as mda
import numpy as np
import argparse

print("start")

parser = argparse.ArgumentParser(description="Extraction of Interaction Energy Streamlined for BaNDyT")
parser.add_argument("--pdb", required=True, help="Filename of system PDB topology")
parser.add_argument("--xtc", required=True, help="Filename of trajectory XTC")
parser.add_argument("--run_id", required=True, help="Prefix name for the run")
parser.add_argument("--top", required=True, help="Filename of topology file")

args = parser.parse_args()

pdb=args.pdb
xtc=args.xtc
run_id=args.run_id
top=args.top
input_dir="../input/"

system=Simul(pdb_file=pdb,xtc_file=xtc,verbose=True)

system.IncludeTopology(path_top=top,include_dir=inp_dir)
system.get_exclusions()
print("hi")
energy=system.get_Energy_perres(start=0,stop=-1,selection="protein")

lj = energy[0]
col = energy[1]

u = mda.Universe(pdb)
protein_residues = u.select_atoms("protein").residues

strings = [f"{res.resname}{res.resnum}" for res in protein_residues]
header_line = ",".join(strings)

np.savetxt(f"../data/{run_id}_interaction_energies.csv", lj+col, delimiter=",", header=header_line, comments="")
