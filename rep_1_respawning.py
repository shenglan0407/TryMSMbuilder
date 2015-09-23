#---------------------------------------------------------------------
# Author: Shenglan Qiao <shenglan@stanford.edu>
# Contributors:
# Copyright (c) 2015, Stanford University
# All rights reserved.
# 
# This script takes simulation data in replicate 1 of the b2AR-alprenolol binding 
# simulation and demonstrates how to determine the starting conformations if the next
# epoch of simulations. 

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import numpy as np
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import utilities as util

import mdtraj as md
from mdtraj.geometry import distance as md_dist

from msmbuilder.featurizer import RawPositionsFeaturizer
from msmbuilder.cluster import KCenters
from msmbuilder.cluster import KMedoids
from msmbuilder.msm import MarkovStateModel
#-----------------------------------------------------------------------------
# matplotlib settings
#-----------------------------------------------------------------------------
mpl.rcParams['font.size'] = 24.0
#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------

sim_path = '/scratch/PI/rondror/MD_simulations/amber/b2AR_ligand_binding/alprenolol/ten_ligands/production/ten_ligands/1/'
sim_file = 'b2AR_ALP_Prod1to9_skip1_reimaged.nc'
topology = '/scratch/PI/rondror/MD_simulations/amber/b2AR_ligand_binding/alprenolol/ten_ligands/system.psf'

# load simulation as mdtraj object
traj = md.load(sim_path+sim_file,top = topology, stride = 100)

# track atoms in ligands
inds_N =[] #indices of N atom
atoms_to_track =['N']
ligands = ['ALP%d'%(ii+1) for ii in range(10)]

for ligand in ligands:
    iis = [atom.index for atom in traj.topology.atoms if atom.element.symbol in atoms_to_track
          and str(atom.residue)==ligand]
    inds_N.append(iis)
    
# track ASP 113
res_to_track = ['ASP113']
target_res = []
for residue in traj.topology.residues:
    if str(residue) in res_to_track:
        if str(residue) in str(target_res):
            pass
        else:
            target_res.append(residue)

ind_res = []
for this_res in target_res:
    ii = [atom.index for atom in this_res.atoms if atom.name in ['OD1']]
    ind_res.append(ii)
atom_pairs = []
for this_lig in inds_N:
    for lig_atom,rec_atom in zip(this_lig,ind_res[0]):
        atom_pairs.append([lig_atom,rec_atom])
atom_pairs = np.array(atom_pairs)
print 'We are tracking distances between the following pairs of atoms:'     
print atom_pairs