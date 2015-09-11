#---------------------------------------------------------------------
# Author: Shenglan Qiao <shenglan@stanford.edu>
# Contributors:
# Copyright (c) 2015, Stanford University
# All rights reserved.

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import numpy as np
import pickle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import utilities as util

from msmbuilder.cluster import KCenters
from msmbuilder.cluster import KMedoids
from msmbuilder.cluster import KMeans
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

LOAD_STRIDE = None
N_CLUSTER = 10

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
print inds_N
#sequences of coordinates of ligands
sequences_all = []
for this_sim in simulations:
    if use_COM:
        this_seq = util.featurize_RawPos(inds_all,this_sim,average = True)
    else:
        this_seq = util.featurize_RawPos(inds_N,this_sim)
    sequences_all.extend(this_seq)
    
seqfile = '/home/shenglan/TryMSMbuilder/output/sequences'+'_s'+str(LOAD_STRIDE)+'.out'
pickle.dump(sequences_all, open(seqfile,'wb'))
    
KC_clustering = KCenters(n_clusters = N_CLUSTER)
KC_assignments = KC_clustering.fit_predict(sequences_all)
KC_centers = KC_clustering.cluster_centers_

KM_clustering = KCenters(n_clusters = N_CLUSTER)
KM_assignments = KM_clustering.fit_predict(sequences_all)
KM_centers = KM_clustering.cluster_centers_

KC_output_file = '/home/shenglan/TryMSMbuilder/output/KC_centers_c'+str(N_CLUSTER)+'_s'+str(LOAD_STRIDE)+'.out'
KM_output_file = '/home/shenglan/TryMSMbuilder/output/KM_centers_c'+str(N_CLUSTER)+'_s'+str(LOAD_STRIDE)+'.out'
np.savetxt(KC_output_file,KC_centers,fmt = '%10.4g')
np.savetxt(KM_output_file,KM_centers,fmt = '%10.4g')

KC_assign_file = '/home/shenglan/TryMSMbuilder/output/KC_assign_'+str(N_CLUSTER)+'_s'+str(LOAD_STRIDE)+'.out'
KM_assign_file = '/home/shenglan/TryMSMbuilder/output/KM_assign_'+str(N_CLUSTER)+'_s'+str(LOAD_STRIDE)+'.out'
pickle.dump(KC_assignments,open(KC_assign_file,'wb'))
pickle.dump(KM_assignments,open(KM_assign_file,'wb'))

KC_pdb_path = KC_output_file.split('.')[0]+'.pdb' 
KM_pdb_path = KM_output_file.split('.')[0]+'.pdb' 

util.convert_to_pdb(KC_output_file,KC_pdb_path)
util.convert_to_pdb(KM_output_file,KM_pdb_path)


KC_pdb_file = '/home/shenglan/TryMSMbuilder/output/KC_c'+str(N_CLUSTER)+'.pdb'
util.convert_sequences_to_pdb(seqfile,KC_assign_file,KC_pdb_file)

KM_pdb_file = '/home/shenglan/TryMSMbuilder/output/KM_c'+str(N_CLUSTER)+'.pdb'
util.convert_sequences_to_pdb(seqfile,KM_assign_file,KM_pdb_file)

