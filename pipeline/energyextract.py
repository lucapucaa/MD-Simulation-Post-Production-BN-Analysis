"""This module contains function to compute LJ and Coulomb potential
        from an xtc trajectory
"""
from openmm.app import GromacsTopFile
from openmm.app import PDBFile, PME, HBonds
from openmm import unit, NonbondedForce, CustomNonbondedForce, LangevinMiddleIntegrator
from scipy.special import erfcinv
from openmm.app import ForceField, Simulation
from scipy.special import erfc
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
import MDAnalysis as mda
from MDAnalysis.lib.nsgrid import FastNS
import numpy as np
import pandas as pd
from numba import njit
from math import erfc as matherfc
from math import erf as matherf
from math import sqrt
import time

@njit
def get_EperResidue_numba( 
                    positions, 
                    resids,
                    n_res, 
                    nbindices, 
                    Acoef, 
                    Bcoef, 
                    charges,
                    neigh_res,
                    cutoff,
                    beta,                    
                    exc_begin,
                    exc_i,
                    exc_j,
                    exc_qprod,
                    exc_aij,
                    exc_bij):
    """Function to compute per-residue Lennard-Jones and Coulomb energies using
    Numba for optimization. This function uses a neighbor list for atom pairs, 
    and then computes the energies for all atom pairs within those residues.

    Parameters
    ----------
    positions : np.ndarray(n_at, 3)
        Array of atom positions. Must be the positions of all atoms in the system
    resids : np.ndarray(n_at)
        Array of residue indices for each atom.
    n_res :  int
        Total number of residues.
    nbindices : np.ndarray
        Array of lj nonbonded indices for each atom.
    Acoef : np.ndarray
        A coefficients for Lennard-Jones potential.
    Bcoef : np.ndarray
            B coefficients for Lennard-Jones potential.
    charges : np.ndarray(n_at)
        Array of charges for each atom.
    neigh_res : np.ndarray
        Array of neighboring atom pairs.
    cutoff : float
        Cutoff distance for nonbonded interactions.
    beta : float
        Parameter for the error function in Coulomb potential.

    Returns
    -------
    _type_
        _description_
    """
    cutoff_2 = cutoff*cutoff
    cutoff_6 = cutoff_2*cutoff_2*cutoff_2
    inv_cutt6 = 1.0/cutoff_6 
    LJ_mat = np.zeros((n_res, n_res))
    Coul_mat = np.zeros((n_res, n_res))
    for k in range(neigh_res.shape[0]):
        i = neigh_res[k,0]
        j = neigh_res[k,1]
        if resids[i] != resids[j]:
            dx = positions[i,0] - positions[j,0]
            dy = positions[i,1] - positions[j,1]
            dz = positions[i,2] - positions[j,2]
            r2 = dx*dx + dy*dy + dz*dz
            if r2 < cutoff_2: #May delete
                r = sqrt(r2)
                inv_r = 1.0/r
                inv_r6 = inv_r**6
                aij = Acoef[nbindices[i], nbindices[j]]
                bij = Bcoef[nbindices[i], nbindices[j]]
                qprod = charges[i] * charges[j]

                Coul = 138.935456 * qprod * (inv_r) * matherfc(beta * r) #-1/cutoff)  # Coulomb's constant in kJ·nm/(mol·e²)
                if i > j:
                    ii = j
                    jj = i
                else:  
                    ii = i
                    jj = j
                # Lookup in the encoded list
                start_idx = exc_begin[ii]
                final_idx = exc_begin[ii+1]
                # Lookup in the elements of the encoded list
                for s in range(start_idx,final_idx):
                    if exc_j[s] == jj and exc_i[s] == ii:
                        aij = 0#exc_aij[s] # Gromacs excludes all the interactions 1-2, 1-3, 1-4
                        bij = 0#exc_bij[s] # Gromacs excludes all the interactions 1-2, 1-3, 1-4 for LJ potentials
                        #Coul = 138.935456 * qprod * (inv_r) * matherfc(beta * r)
                        Coul = -138.935456 * qprod * (inv_r) * matherf(beta * r)  # Coulomb correction for PME (May remove eventually) 
                        break


                LJ_val = ((inv_r6 * aij)**2 - inv_r6 * bij) - ((aij*inv_cutt6)**2 - (bij*inv_cutt6))  # Lennard-Jones potential with cutoff
                LJ_mat[resids[i]-1, resids[j]-1] += LJ_val
                LJ_mat[resids[j]-1, resids[i]-1] += LJ_val

                Coul_mat[resids[i]-1, resids[j]-1] += Coul
                Coul_mat[resids[j]-1, resids[i]-1] += Coul
    return LJ_mat, Coul_mat



# Running with numba and residue optimization
# Warning, this potentally only work is the residues and numbered from 1-n_res without any jump
@njit
def get_EperResidue_numba_res( 
                    positions, 
                    resids,
                    n_res,
                    res_limits, 
                    nbindices, 
                    Acoef, 
                    Bcoef, 
                    charges,
                    neigh_res,
                    cutoff,
                    beta,
                    exc_begin,
                    exc_i,
                    exc_j,
                    exc_qprod,
                    exc_aij,
                    exc_bij,): # Notice that exc arrays have not been used because I realized that when computing energy with 
                                # gromacs, those corrections are not added. The corresponding exceptions are just deleted 
    LJ_mat = np.zeros((n_res, n_res))
    Coul_mat = np.zeros((n_res, n_res))
    """Function to compute per-residue Lennard-Jones and Coulomb energies using 
    Numba for optimization. This function uses a neighbor list only 
    for residues, and then computes the energies for all atom pairs within those residues.
    In this way, we avoid computing energies for atom pairs that are not in neighboring residues,
    which can significantly reduce the number of computations for large systems.

    Parameters
    ----------
    positions : np.ndarray(n_at, 3)
        Array of atom positions. Must be the positions of all atoms in the system
    resids : np.ndarray(n_at)
        Array of residue indices for each atom.
    n_res : int
        Total number of residues.
    res_limits : np.ndarray(n_res)
        Array of indices that mark the end of each residue in the positions array.
    nbindices : np.ndarray
        Array of nonbonded indices for each atom.
    Acoef : np.ndarray
        Array of A coefficients for Lennard-Jones potential.
    Bcoef : np.ndarray  
         Array of B coefficients for Lennard-Jones potential.
    charges : np.ndarray(n_at)
        Array of charges for each atom.
    neigh_res : np.ndarray
        Array of neighboring residue pairs.
    cutoff : float
        Cutoff distance for nonbonded interactions.
    beta : float
        Parameter for the error function in Coulomb potential.
    exc_begin : np.ndarray
        Array of indices that mark the beginning of the exceptions for each atom.
    exc_i : np.ndarray
        Array of atom indices for the first atom in each exception.
    exc_j : np.ndarray
        Array of atom indices for the second atom in each exception.
    exc_qprod : np.ndarray
        Array of charge products for each exception.
    exc_aij : np.ndarray
        Array of A^2 coefficients for each exception.
    exc_bij : np.ndarray
        Array of B coefficients for each exception.




    Returns
    -------
    LJ_mat : np.ndarray(n_res, n_res)
        Matrix of Lennard-Jones energies between residues.
    Coul_mat : np.ndarray(n_res, n_res)
        Matrix of Coulomb energies between residues.
    """
    
    # May consired to ask this quantities as input
    cutoff_2 = cutoff*cutoff
    cutoff_6 = cutoff_2*cutoff_2*cutoff_2
    inv_cutt6 = 1.0/cutoff_6



    for k in range(neigh_res.shape[0]):
        i = neigh_res[k,0] # These will be residue indices, not atom indices
        j = neigh_res[k,1] # These will be residue indices, not atom indices
        
        # Get the ids of the atoms in the residues
        begin_i = res_limits[i-1] if i > 0 else 0
        begin_j = res_limits[j-1] if j > 0 else 0


        if resids[begin_i] != resids[begin_j]:  # Check if the residues are different
            for l in range(begin_i, res_limits[i]):
                for m in range(begin_j, res_limits[j]):
                    dx = positions[l,0] - positions[m,0]
                    dy = positions[l,1] - positions[m,1]
                    dz = positions[l,2] - positions[m,2]
                    r2 = (dx*dx + dy*dy + dz*dz)
                    if r2 < cutoff_2:
                        r = sqrt(r2)
                        #print(r)
                        inv_r = 1.0/r
                        inv_r6 = inv_r**6
                        aij = Acoef[nbindices[l], nbindices[m]]
                        bij = Bcoef[nbindices[l], nbindices[m]]
                        qprod = charges[l] * charges[m]

                        Coul = 138.935456 * qprod * (inv_r) * matherfc(beta * r)#-1/cutoff)  # Coulomb's constant in kJ·nm/(mol·e²)

                        # Only add bonded exceptions if the residues are adjacent
                        if resids[l] == resids[m] + 1  or resids[l] == resids[m]-1:
                            if l > m:
                                ii = m
                                jj = l
                            else:  
                                ii = l
                                jj = m
                            # Lookup in the encoded list
                            start_idx = exc_begin[ii]
                            final_idx = exc_begin[ii+1]
                            # Lookup in the elements of the encoded list
                            for s in range(start_idx,final_idx):
                                if exc_j[s] == jj and exc_i[s] == ii:
                                    aij = 0#exc_aij[s] # Gromacs excludes all the interactions 1-2, 1-3, 1-4
                                    bij = 0#exc_bij[s] # Gromacs excludes all the interactions 1-2, 1-3, 1-4 for LJ potentials
                                    #Coul = 138.935456 * qprod * (inv_r) * matherfc(beta * r)
                                    Coul = -138.935456 * qprod * (inv_r) * matherf(beta * r)  # Coulomb correction for PME (May remove eventually) 
                                    break
                        
                        

                        LJ_val = ((inv_r6 * aij)**2 - inv_r6 * bij) - ((aij*inv_cutt6)**2 - (bij*inv_cutt6))  # Lennard-Jones potential with cutoff


                        # New store approach:

                        LJ_mat[resids[l], resids[m]] += LJ_val
                        LJ_mat[resids[m], resids[l]] += LJ_val
                        Coul_mat[resids[l], resids[m]] += Coul
                        Coul_mat[resids[m], resids[l]] += Coul

    return LJ_mat, Coul_mat











class Simul():
    def __init__(self, 
                 pdb_file, 
                 xtc_file, 
                 tpr_file = None,
                 verbose = False):

        self.xtc_file = xtc_file
        self.tpr_file = tpr_file if tpr_file is not None else None
        self.pdb_file = pdb_file
        self.cutoff = 1.2
        self.emtol = 1e-5
        self.beta = erfcinv(self.emtol)/self.cutoff
        self.universe = mda.Universe(self.pdb_file, self.xtc_file)

        # Store numbre of atoms
        self.n_atoms = len(self.universe.atoms)

        # Store the resindices (Unique and starting from 0)        
        self.resids = self.universe.atoms.resindices

        # Store the number of residues
        self.n_res = len(self.universe.residues.resindices)

        # Turns on/off printing
        self.verbose = verbose

        # Get the limits of the residues in the atom list, to be used in the energy calculation
        resindices = self.universe.atoms.resids
        res_limits = np.flatnonzero(np.diff(resindices)) + 1
        res_limits = np.append(res_limits, len(resindices))

        self.resid_limits = np.array(res_limits)

        # Create an inverse map of the original residues with the resindices
        # Notice that this map is not a function since original resids can have identical resids
        self.original_resids = self.universe.residues.resids 
        self.map_resids = {int(new_resid): int(resid_or)  for resid_or, new_resid in zip(self.original_resids, self.universe.residues.resindices)}

        

        

    def IncludeTopology(self, 
                        path_top = None, 
                        include_dir = None, 
                        openmm_param = False,
                        simulation = False):
        """Extract topology information from the GROMACS files or parametrize them
        with OpenMM ForceField under the CHARMM36 force field. Then, include the 
        topology information into the MDAnalysis Universe.

        Parameters
        ----------
        path_top : str, optional
            Path to the topology file (.top), by default None
        include_dir : str, optional
            Path to the attachment in the .top file, by default None
        openmm_param : bool, optional
            If true add parameters from charmm36 forcefield using openmm, by default False
        """
        
        if openmm_param: # Parametrizes the topology using openmm
            pdb = PDBFile(self.pdb_file)
            forcefield = ForceField("charmm36.xml")
            self.system = forcefield.createSystem(pdb.topology, 
                                                  nonbondedMethod=PME, 
                                                  constraints=HBonds, 
                                                  nonbondedCutoff=1.2*unit.nanometer)
        else:  # uses gromacs files to build the topology/parametrization
            pdb = PDBFile(self.pdb_file)
            top = GromacsTopFile(path_top, 
                                 periodicBoxVectors=pdb.topology.getPeriodicBoxVectors(), 
                                 includeDir = include_dir)
            self.system = top.createSystem(nonbondedMethod=PME, 
                                       constraints=HBonds, 
                                       nonbondedCutoff=1.2*unit.nanometer)
        self.simulation_obj = None
        if simulation: # Create a simulation object needed to compute things wiht openmm
            integrator = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.004*unit.picoseconds)
            self.simulation_obj = Simulation(pdb.topology, self.system, integrator)
        
        # Variable not used yet, still thinking if worth it or not
        self.optimization_methods = ["residue-optimized", "atom-optimized", "openmm", "python"]

        # Obtain topology information
        self.forces = self.system.getForces()
        self.custom_nb = None
        self.simulation = simulation if simulation else None

        # When using gromcas topology, openmm build the LJ parameters and exceptions in custom non bonded force
        for force in self.forces:
            if isinstance(force, NonbondedForce):
                self.nb = force
            if isinstance(force, CustomNonbondedForce):
                self.custom_nb = force

        # Get the information from the forcefield
        n = self.nb.getNumParticles()
        attributes = {"charges" : np.empty(n),
                      "radii" : np.empty(n),
                      "epsilons" : np.empty(n),
                      "nbindex" : np.empty(n, dtype=int)}
        
        # MAtrices to store the LJ parameters
        self.Acoef = None
        self.Bcoef = None
        if self.custom_nb is not None:
            acoef = self.custom_nb.getTabulatedFunction(0)
            bcoef = self.custom_nb.getTabulatedFunction(1)

            xsize, ysize, val_a = acoef.getFunctionParameters()
            xsize, ysize, val_b = bcoef.getFunctionParameters()
            self.Acoef = np.array(val_a).reshape((xsize, ysize))
            self.Bcoef = np.array(val_b).reshape((xsize, ysize))
        for i in range(n):
            q, s, e = self.nb.getParticleParameters(i)
            attributes["charges"][i] = q.value_in_unit(unit.elementary_charge)
            attributes["radii"][i] = s.value_in_unit(unit.nanometers)
            attributes["epsilons"][i] = e.value_in_unit(unit.kilojoule_per_mole)
            attributes["nbindex"][i] = self.custom_nb.getParticleParameters(i)[0]
        
        self.charges = attributes["charges"]
        self.radii = attributes["radii"]
        self.epsilons = attributes["epsilons"]
        self.nbindices = attributes["nbindex"]
        self.resindices = self.universe.atoms.resindices

        for attr in attributes:
            self.universe.add_TopologyAttr(attr, attributes[attr])


    # Compute exclusions for the nonbonded interactions
    # Still not implemented in energy calculation
    def get_exclusions(self):

        exclusions_residues = set()
        n_exceptions = self.nb.getNumExceptions()

        # Create arrays to hold the exception parameters
        # Machinery for adding exceptions
        id_i = np.empty(n_exceptions, dtype=np.int32)
        id_j = np.empty(n_exceptions, dtype=np.int32)
        sigmas = np.empty(n_exceptions)
        epsilons = np.empty(n_exceptions)
        qprod = np.empty(n_exceptions)

        exc = []


        for k in range(self.nb.getNumExceptions()):
            p1, p2, chargeprod, sigma, epsilon = self.nb.getExceptionParameters(k)
            id_i[k] = p1
            id_j[k] = p2
            exc.append([p1, p2])

            sigmas[k] = sigma.value_in_unit(unit.nanometer)
            epsilons[k] = epsilon.value_in_unit(unit.kilojoule_per_mole)
            qprod[k] = chargeprod.value_in_unit(unit.elementary_charge**2)

        
        self.exc = np.array(exc)
        order = np.lexsort((id_j, id_i))

        self.exc_i = id_i[order]
        self.exc_j = id_j[order]
        sigmas = sigmas[order]
        self.exc_sigma6 = sigmas*sigmas*sigmas*sigmas*sigmas*sigmas*sigmas

        self.exc_epsilons = epsilons[order]
        self.exc_qprod = qprod[order]
        self.exc_aij = np.sqrt(4*self.exc_epsilons*self.exc_sigma6*self.exc_sigma6)
        self.exc_bij = 4*self.exc_epsilons*self.exc_sigma6

        counts = np.bincount(id_i, minlength=self.n_atoms)
        begin = np.zeros(len(counts)+1, dtype=np.int32)
        begin[1:] = np.cumsum(counts)
        self.exc_begin = begin


        return exclusions_residues
    
    def get_Energy_perres(self, start = 0, stop = -1, step = 1, selection = "protein"):
        """Compute the energy of the system from the trajectory

        Parameters
        ----------
        start : int, optional
            Start frame, by default 0
        stop : int, optional
            Stop frame, by default -1 (last frame)
        step : int, optional
            Step size, by default 1
        selection: str, optional
            Atom selection string for MDAnalysis, by default "protein"

        Returns
        -------
        result_lj : np.ndarray
            Array of Lennard-Jones energies for each frame
        result_coul : np.ndarray
            Array of Coulomb energies for each frame
        """
        
        print("started_perres")
        sel_string = selection if isinstance(selection, str) else " or ".join([f"({sel})" for sel in selection])
        sel_indices = self.universe.select_atoms(sel_string).residues.resindices


        select_indices = False
        if isinstance(selection, list):
            select_indices = [] 
            for select in selection:
                select_indices.append(self.universe.select_atoms(select).residues.resindices)
        


        #print(sel_string, "sel_string")
        #print(sel_indices, "resindices")

        all_atoms = self.universe.select_atoms("protein") 
        ca_atoms = all_atoms.select_atoms(sel_string) # Only used for residue-optimized energy calculation
        ca_atoms = all_atoms.select_atoms(f"({sel_string}) and name CA") # Only used for residue-optimized energy calculation
        ca_atoms_ids = ca_atoms.indices

        #print(ca_atoms_ids, "ca_atoms_ids")
        Lj_data = []
        Coul_data = []
        for ts in self.universe.trajectory[start:stop:step]:
            if self.verbose:
                print(f"Frame {ts.frame}")

            # Get positions of all atoms in the system
            positions = all_atoms.positions/10 # Convert from Angstroms (MDAnalysis default units) to nanometers

            t1 = time.perf_counter()
            energies_per_res = self.call_Energy_resbased(positions=positions, 
                                                         ca_atoms_ids=ca_atoms_ids, 
                                                         selections = select_indices, 
                                                         selected_residues = sel_indices)
            t2 = time.perf_counter()
            print(f"Frame {ts.frame}: Energy calculation took {t2-t1:.4f} seconds####")

            if isinstance(selection, str):
                lj = energies_per_res[0][sel_indices][:, sel_indices]
                coul = energies_per_res[1][sel_indices][:, sel_indices]
            elif len(select_indices) == 2:
                lj = energies_per_res[0][select_indices[0]][:, select_indices[1]]
                coul = energies_per_res[1][select_indices[0]][:, select_indices[1]]
            else:
                print("The code does not support interaction between three atom groups")

            Lj_data.append(np.sum(lj, axis=1)) # Sum over all residues to get total energy for the frame
            Coul_data.append(np.sum(coul, axis = 1)) # Sum over all residues to get total energy for the frame

        result_lj = np.asarray(Lj_data, dtype = np.float32)
        
        result_coul = np.asarray(Coul_data, dtype = np.float32)
        
        return result_lj, result_coul
    

    def call_Energy_atombased(self,positions):

        neigh_res = FastNS(1.3, positions, box = self.universe.dimensions)

        # Get neighboring residue pairs
        neigh_res = neigh_res.self_search()
        neigh_res = neigh_res.get_pairs()

        energies = get_EperResidue_numba( 
                        positions, 
                        self.resids,
                        self.n_res, 
                    self.nbindices, 
                    self.Acoef, 
                    self.Bcoef, 
                    self.charges,
                    neigh_res,
                    self.cutoff,
                    self.beta,                    
                    self.exc_begin,
                    self.exc_i,
                    self.exc_j,
                    self.exc_qprod,
                    self.exc_aij,
                    self.exc_bij)
        
        return energies


    
    def call_Energy_resbased(self,positions, ca_atoms_ids, selections = False, selected_residues = None):

        neigh_res = FastNS(3, positions[ca_atoms_ids], box = self.universe.dimensions)
        
        # Get neighboring residue pairs
        neigh_res = neigh_res.self_search()
        neigh_res = neigh_res.get_pairs()
        if selections:
            # Filter intra residue interactions and interactions between residue groups
            group_mask = np.zeros(len(selected_residues), dtype=np.uint8)
            idx1 = np.searchsorted(selected_residues, selections[0])
            idx2 = np.searchsorted(selected_residues, selections[1])
            group_mask[idx1] |= 1
            group_mask[idx2] |= 2
            g0 = group_mask[neigh_res[:,0]]
            g1 = group_mask[neigh_res[:,1]]
            mask = (
                    (((g0 & 1) != 0) & ((g1 & 2) != 0)) |
                    (((g0 & 2) != 0) & ((g1 & 1) != 0))
                    )
            
            neigh_res = neigh_res[mask]

        neigh_res[:,0] = selected_residues[neigh_res[:,0]]
        neigh_res[:,1] = selected_residues[neigh_res[:,1]]


        energies_per_res = get_EperResidue_numba_res(positions,
                                        self.resids,
                                        self.n_res, 
                                        self.resid_limits, 
                                        self.nbindices, 
                                        self.Acoef, 
                                        self.Bcoef, 
                                        self.charges,
                                        neigh_res,
                                        self.cutoff,
                                        self.beta,
                self.exc_begin,
                self.exc_i,
                self.exc_j,
                self.exc_qprod,
                self.exc_aij,
                self.exc_bij,)        

        return energies_per_res
    

    def get_EperResidue(self, frame, positions, resids, nbindices):
        """First attemp to compute per-residue Lennard-Jones and Coulomb energies using MDAnalysis and FastNS for neighbor search.


        Parameters
        ----------
        frame : int
            Frame number
        positions : np.ndarray
            Array of atom positions. Must be the positions of all atoms in the system
        resids : np.ndarray
            Array of residue indices for each atom.
        nbindices : np.ndarray
            Array of lj values fro each atom pair.

        Returns
        -------
        LJ_mat : np.ndarray
            Matrix of Lennard-Jones energies between residues.
        Coul_mat : np.ndarray
            Matrix of Coulomb energies between residues.
        """
        #print(self.cutoff)
        neigh = FastNS(self.cutoff, positions/10, box = self.universe.dimensions)
        neigh_res = neigh.self_search()
        #print(type(neigh_res.get_pairs()))

        n_res = len(set(resids))

        LJ_mat = np.zeros((n_res, n_res))
        Coul_mat = np.zeros((n_res, n_res))


        for i,j in neigh_res.get_pairs():
            #print(i,j,resids[i], resids[j])
            if resids[i] == resids[j]:
                continue
            else:
                r = np.linalg.norm(positions[i] - positions[j])/10
                inv_r = 1.0/r
                inv_r6 = inv_r**6

                aij = self.Acoef[nbindices[i], nbindices[j]]
                bij = self.Bcoef[nbindices[i], nbindices[j]]

                qprod = self.universe.atoms.charges[i] * self.universe.atoms.charges[j]
                LJ_val = ((inv_r6 * aij)**2 - inv_r6 * bij) - ((aij/self.cutoff**6)**2 - (bij/self.cutoff**6))  # Lennard-Jones potential with cutoff
                LJ_mat[resids[i]-1, resids[j]-1] += LJ_val
                LJ_mat[resids[j]-1, resids[i]-1] += LJ_val

                Coul = 138.935456 * qprod * (inv_r) * erfc(self.beta * r)#-1/cutoff)  # Coulomb's constant in kJ·nm/(mol·e²)
                Coul_mat[resids[i]-1, resids[j]-1] += Coul
                Coul_mat[resids[j]-1, resids[i]-1] += Coul

        return LJ_mat, Coul_mat
        





    def get_Energy_opmm(self, start = 0, stop = -1, step = 1):
        """Compute the energy of the system from the trajectory using OpenMM

        Parameters
        ----------
        start : int, optional
            Start frame, by default 0
        stop : int, optional
            Stop frame, by default -1 (last frame)
        step : int, optional
            Step size, by default 1

        Returns
        -------
        dict
            Dictionary with the energies for each frame
        """
        energies = {"frame" : [], "LJ" : [], "Coulomb" : []}

        
        

        for i, force in  enumerate(self.forces):
            print(force)
            try:
                print(force.getEnergyFunction())
            except:
                continue
            if isinstance(force, CustomNonbondedForce):
                print(i,force)
                cforce = force
        print(cforce.getEnergyFunction())        
        
        cforce.addInteractionGroup({100},{107})


        
        #for ts in self.universe.trajectory[start:stop:step]:
        #    self.simulation_obj.context.setPositions(self.universe.atoms.positions/10)



        #    state = self.simulation_obj.context.getState(getEnergy=True)
        #    energies["frame"].append(ts.frame)
        #    energies["LJ"].append(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
        #    energies["Coulomb"].append(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
        #    if ts.frame == 10:
        #        print(state.getPotentialEnergy().count)
        #return energies









"""
import matplotlib.pyplot as plt


fig,ax = plt.subplots(1,2)
frames = energies[0][:,0]
lj = energies[0][:,6]
coul = energies[1][:,6]

print(system.map_resids)

#lj1 = energies2[0][:,5]
#coul1 = energies2[1][:,5]
print(frames, lj, coul, energies[0].shape, energies[1].shape)
cb2 = ax[0].plot(frames, coul)
ax[0].set_title("Coulomb Energy")
cb1 = ax[1].plot(frames, lj)
ax[1].set_title("Lennard-Jones Energy")
#ax[2].plot(frames, lj)
#ax[3].plot(frames, coul)


 
 
plt.show()


import pandas as pd
fig, ax = plt.subplots(3,1)
data = pd.read_csv("../../testing_env/job_1/energy_213.xvg", 
                   skiprows=14,
                   comment='@',
                   sep=r"\s+", 
                   header=None,
                   names=["time", "coul", "lj"])



ax[0].plot(energies["frame"], energies["LJ"], label = "LJ-own")
ax[1].plot(energies["frame"], energies["Coulomb"], label = "Coulomb-own")
ax[0].plot(data["time"]/3000, data["lj"], label = "gromcas")
ax[1].plot(data["time"]/3000, data["coul"], label = "gromacs")
ax[1].legend()

ax[0].legend()

plt.show()

print(np.max(system.universe.atoms.nbindices))
"""                 
