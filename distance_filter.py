#---------------------------------------------------------------------
# Author: Shenglan Qiao <shenglan@stanford.edu>
# Contributors:
# Copyright (c) 2015, Stanford University
# All rights reserved.

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import numpy as np
import re
import pickle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import utilities as util

from msmbuilder.cluster import KCenters
from msmbuilder.cluster import KMedoids
from msmbuilder.msm import MarkovStateModel

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------

#global parameters
use_COM = False
parent_dir = '/scratch/PI/rondror/DesRes-Simulations/B2AR_PNAS_2011a_Dror_Desres/DESRES-Trajectory_pnas2011a-B-no-water-no-lipid.1'
number_of_sims = 1

run_dirs = []
for ii in range(number_of_sims):
    run_dirs.append('/DESRES-Trajectory_pnas2011a-B-'+str(ii)+'-no-water-no-lipid/pnas2011a-B-'+str(ii)+'-no-water-no-lipid')
dir_top = '/home/shenglan/topologies'
times_path = parent_dir+run_dirs[0]+'/'+run_dirs[0].split('/')[-1]+'_times.csv'

LOAD_STRIDE = 1000


#load list of mdtraj objects
simulations = []
for this_run_dir in run_dirs:
    simulations.append(util.load_trajs(this_run_dir,parent_dir,dir_top,load_stride = LOAD_STRIDE))

#track atoms in ligands
inds_N =[] #indices of N atom
atoms_to_track =['N']
ligands = [residue for residue in simulations[0][0].topology.chain(0).residues][0:10]

for ligand in ligands:
    iis = [atom.index for atom in simulations[0][0].topology.chain(0).atoms if atom.element.symbol in atoms_to_track
          and atom.residue == ligand]
    inds_N.append(iis)

#track residue Asp113
first_res = str(simulations[0][0].topology.chain(0).residue(0))
APSs = [residue.index for residue in simulations[0][0].topology.chain(0).residues if residue.name == 'ASP']

match = re.match(r"([a-z]+)([0-9]+)", first_res, re.I)
if match:
    items = match.groups()
first_res_number = int(items[-1])
target_res = simulations[0][0].topology.chain(0).residue(113-18)
res_ind = [[atom.index for atom in target_res.atoms if atom.name in ['OD1']]]

atom_pairs = []
for this_lig in inds_N:
    print this_lig
    print res_ind
    for lig_atom,rec_atom in zip(this_lig,res_ind[0]):
        atom_pairs.append([lig_atom,rec_atom])
print 'We are tracking distances between the following pairs of atoms:'     
print atom_pairs
    