import os
import MDAnalysis as mda
from MDAnalysis.transformations import wrap

PROJECT_ROOT = "/home/lchill/work/pproduct_analysis"
topology_file = os.path.join(PROJECT_ROOT, "input/7o7k_ATP.pdb")
output_pdb = os.path.join(PROJECT_ROOT, "data/extracted_frame.pdb")

TARGET_FRAME = 100 

# Candidate trajectory locations
possible_trajectories = [
    os.path.join(PROJECT_ROOT, "data/replica_1.xtc"),
    os.path.join(PROJECT_ROOT, "data/final.xtc"),
    os.path.join(PROJECT_ROOT, "data/replica_data_1/final.xtc"),
    os.path.join(PROJECT_ROOT, "input/results/combined_1.xtc")
]

trajectory_file = None
for candidate in possible_trajectories:
    if os.path.exists(candidate):
        trajectory_file = candidate
        break

if not trajectory_file:
    print("❌ Could not find trajectory file in any of these locations:")
    for path in possible_trajectories:
        print(f"   - {path}")
    exit(1)

print(f"📖 Loading Topology: {topology_file}")
print(f"🎬 Loading Trajectory: {trajectory_file}")

u = mda.Universe(topology_file, trajectory_file)

if TARGET_FRAME >= len(u.trajectory):
    print(f"❌ Target frame {TARGET_FRAME} exceeds total available frames ({len(u.trajectory)}).")
    exit(1)

# Jump to frame
u.trajectory[TARGET_FRAME]

# Apply PBC wrap transformation to prevent PDB coordinate overflow errors
transform = wrap(u.atoms)
u.trajectory.add_transformations(transform)

atoms = u.select_atoms("all")
atoms.write(output_pdb)

print(f"✅ Successfully saved Frame {TARGET_FRAME} to: {output_pdb}")
