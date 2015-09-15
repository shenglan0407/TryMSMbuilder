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

from mdtraj.geometry import distance as md_dist

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
    for lig_atom,rec_atom in zip(this_lig,res_ind[0]):
        atom_pairs.append([lig_atom,rec_atom])
atom_pairs = np.array(atom_pairs)
print 'We are tracking distances between the following pairs of atoms:'     
print atom_pairs

distances = []
for this_sim in simulations:
    this_dist = []
    for this_traj in this_sim:
        
        this_dist.extend(md_dist.compute_distances(this_traj,atom_pairs))
    this_dist = np.array(this_dist)
    distances.append(this_dist)

sequences_all = []
for this_sim in simulations:
    if use_COM:
        this_seq = util.featurize_RawPos(inds_all,this_sim,average = True)
    else:
        this_seq = util.featurize_RawPos(inds_N,this_sim)
    sequences_all.extend(this_seq)
    
# convert to pdb so I can view in vmd
# will use chian ID to store if ligand is in bulk or near receptor

BULK_CUTOFF = 2.0 #unit is nm

assignments = []
this_assign = []
for this_sim in distances:
    natoms = this_sim.shape[1]
    for ii in range(natoms):
        this_assign = []
        for this_frame in this_sim[:,ii]:
            if this_frame > BULK_CUTOFF:
                #B for bulk
                this_assign.append('B')
            else:
                #R for receptor
                this_assign.append('R')
        assignments.append(this_assign)
    