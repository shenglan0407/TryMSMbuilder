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
import copy
import os

import mdtraj as md
from mdtraj.geometry import distance as md_dist

from msmbuilder.cluster import KCenters
from msmbuilder.cluster import KMedoids
from msmbuilder.msm import MarkovStateModel
from msmbuilder.tpt import mfpt
#-----------------------------------------------------------------------------
# matplotlib settings
#-----------------------------------------------------------------------------
mpl.rcParams['font.size'] = 24.0
#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------


N_CLUSTER = 60
LOAD_STRIDE = None
BULK_CUTOFF = 4.0 #unit is nm
N_LIGAND = 10*9

sim_path = '/scratch/PI/rondror/MD_simulations/amber/b2AR_ligand_binding/alprenolol/ten_ligands/production/ten_ligands/reimaged_trajs2/'
topology = '/scratch/PI/rondror/MD_simulations/amber/b2AR_ligand_binding/alprenolol/ten_ligands/system.psf'

# load simulation as mdtraj object
simulations = []
for this_file in os.listdir(sim_path):
    if this_file.endswith('.nc'):
        this_traj = md.load(sim_path+this_file,top = topology, stride = LOAD_STRIDE)
        simulations.append(this_traj)
print len(simulations)
print simulations[-1]

# gather some properties of the trajectories
print ("There are %d replicates in the simulation." % len(simulations))
n_frames = []

for this_sim in simulations:
    n_frames.append(this_sim.n_frames)
n_frames = np.array(n_frames)
time_step = simulations[0].timestep * 1e-3 # unit is ns
# time_step2 = simulations[4].timestep * 1e-3
# print time_step
# print time_step2
total_sim_time = np.sum(n_frames)*time_step

print ("The time step between frames is %.2f." % time_step)
print ("Total simulation time is %.2f ns." % total_sim_time)
print ("The longest simulation is %.2f ns and the shortest %.5f ns." % (np.max(n_frames)*time_step,np.min(n_frames)*time_step))
print ("Total number of frames is %d." % np.sum(n_frames))

print ("We are building an MSM with %d clusters." % N_CLUSTER)



# track atoms in ligands
inds_N =[] #indices of N atom
atoms_to_track =['N']
ligands = ['ALP%d'%(ii+1) for ii in range(10)]

for ligand in ligands:
    iis = [atom.index for atom in simulations[0].topology.atoms if atom.element.symbol in atoms_to_track
          and str(atom.residue)==ligand]
    inds_N.append(iis)
    
# track ASP 113
res_to_track = ['ASP113']
target_res = []
for residue in simulations[0].topology.residues:
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

# comput and record distances between ligand and binding pocket ASP113
distances = []
for this_sim in simulations:
    for this_atom_pair in atom_pairs:
        this_lig = []
        for this_traj in [this_sim]: 
            this_lig.extend(md_dist.compute_distances(this_traj,[this_atom_pair]))
        distances.append(this_lig)

dist_path = '/home/shenglan/TryMSMbuilder/output/ten_ligands/dist_to_binding'\
+'_s'+str(LOAD_STRIDE)+'.out'
pickle.dump(distances,open(dist_path,'wb'))

# get N positions
sequences_all = []
for this_sim in simulations:
    this_seq = util.featurize_RawPos(inds_N,[this_sim])
    sequences_all.extend(this_seq)
seq_path = '/home/shenglan/TryMSMbuilder/output/ten_ligands/sequences'+'_s'+str(LOAD_STRIDE)+'.out'
pickle.dump(sequences_all,open(seq_path,'wb'))

clustering = KCenters(n_clusters = N_CLUSTER)
geo_assign = clustering.fit_predict(sequences_all)
centers = clustering.cluster_centers_

geo_assign_path = '/home/shenglan/TryMSMbuilder/output/ten_ligands/KC_geoassign_c' \
+str(N_CLUSTER)+'_s'+str(LOAD_STRIDE)+'.out'
pickle.dump(geo_assign,open(geo_assign_path,'wb'))

micro_msm = MarkovStateModel(lag_time=1, reversible_type = 'transpose', 
ergodic_cutoff = 'off'
,verbose=True).fit(geo_assign)

msm_path = '/home/shenglan/TryMSMbuilder/output/ten_ligands/KC_msm_c'+str(N_CLUSTER)+ \
'_s'+str(LOAD_STRIDE)+'.out'
pickle.dump(micro_msm,open(msm_path,'wb'))

# map assignments
print('There are %d microstates in msm' % micro_msm.n_states_)

raw_clusters = []
for this_assign in geo_assign:
    raw_clusters.extend(np.unique(this_assign))
raw_clusters = np.unique(np.array(raw_clusters))
print('There are %d clusters in the original geometric clustering.'%len(raw_clusters))

#generate a mapping of geometric clustering to microstates
raw_to_micro_mapping = copy.copy(micro_msm.mapping_)
N_MICRO = len(raw_to_micro_mapping)
dummy_assign = N_MICRO
for ii in range(len(raw_clusters)):
    if ii in raw_to_micro_mapping.keys():
        pass
    else:
        raw_to_micro_mapping[ii] = dummy_assign
        #dummy_assign = dummy_assign+1
        
micro_assign = []

for this_assign in geo_assign:
    this_list = np.zeros(len(this_assign),dtype = int)
    micro_assign.append(this_list)
    
#map assignments for microstates    
for nn in range(len(micro_assign)):
    for ii in range(len(micro_assign[nn])):
        micro_assign[nn][ii] = raw_to_micro_mapping[geo_assign[nn][ii]]


#check for mismatch
for this_assign, that_assign in zip(micro_assign,geo_assign):
    for this_item,that_item in zip(this_assign,that_assign):
        if this_item == raw_to_micro_mapping[that_item]:
            pass
        else:
            print ('mismatch %d, %d' % (this_item, that_item))

micro_assign_path = '/home/shenglan/TryMSMbuilder/output/ten_ligands/KC_microassign_c'+str(N_MICRO)+'_s'+str(LOAD_STRIDE)+'.out'
pickle.dump(micro_assign,open(micro_assign_path,'wb'))

save_micro_pdb_path = '/home/shenglan/TryMSMbuilder/output/ten_ligands/KC_microassign_c'+str(N_MICRO)+'_s'+str(LOAD_STRIDE)+'.pdb'
util.convert_sequences_to_pdb(seq_path,micro_assign_path,save_micro_pdb_path)


# dictionary; key: microstate id; 
# item: list of distances to binding pocket for points in microstate
clusters_dist = {}

for ii in range(len(distances)):
    if len(distances) == len(micro_assign) and len(distances[ii]) == len(micro_assign[ii]):
        for nn in range(len(distances[ii])):
            if micro_assign[ii][nn] in clusters_dist.keys():
                clusters_dist[micro_assign[ii][nn]].append(distances[ii][nn])
            else:
                clusters_dist.update({micro_assign[ii][nn] : [distances[ii][nn]]})
                
    else:
        print ('something is wrong. distances and micro_assign not matching in dimensions')
        print len(distances[ii]) 
        print len(micro_assign[ii])

# Identify 'sinks' for mfpt calculations
bulk_clusters = []
for ii in clusters_dist.keys():
    if np.average(clusters_dist[ii]) > BULK_CUTOFF:
        bulk_clusters.append(ii)
print "These microstates are designated as bulk states:"
print bulk_clusters
    
mfpt_mat = mfpt.mfpts(micro_msm,sinks = bulk_clusters,lag_time = 1.0)
# print mfpt_mat.shape
# print mfpt_mat
print mfpt_mat[mfpt_mat > 0.0]
bound_clusters = np.where(mfpt_mat> 0.0)[0]
print 'The states to not in bulk are:'
print np.where(mfpt_mat> 0.0)


# calculate respawning probability
log_t = np.log(mfpt_mat[mfpt_mat > 0.0])
print log_t
norm = np.sum(np.log(mfpt_mat[mfpt_mat > 0.0]))
respawn_pr = log_t/np.sum(log_t)

print respawn_pr
print np.sum(respawn_pr)
respawn_pop = np.array(respawn_pr*N_LIGAND,dtype = int)
print respawn_pop

# map microstate to respawning populations
respawning = {}
for ii in range(len(bound_clusters)):
    respawning.update({bound_clusters[ii]:respawn_pop[ii]})
respawning_path = '/home/shenglan/TryMSMbuilder/output/ten_ligands/respawning_c'+str(N_MICRO)+'_s'+str(LOAD_STRIDE)+'.out'
pickle.dump(respawning,open(respawning_path,'wb'))

print respawning