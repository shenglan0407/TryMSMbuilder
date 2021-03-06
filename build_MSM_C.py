#---------------------------------------------------------------------
# Author: Shenglan Qiao <shenglan@stanford.edu>
# Contributors:
# Copyright (c) 2015, Stanford University
# All rights reserved.
# 
# This script take all the simulation data in condition C to build MSM. The MSM saved here 
# does not have an ergodic cutoff, i.e. all geometric clusters are used to make the 
# transition matrix even if some are not strongly connected to the network. There is also
# a choice of adding a prior.

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import numpy as np
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import utilities as util

from mdtraj.geometry import distance as md_dist

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

#global parameters
use_COM = False
parent_dir = '/scratch/PI/rondror/DesRes-Simulations/B2AR_PNAS_2011a_Dror_Desres/DESRES-Trajectory_pnas2011a-C-no-water-no-lipid.1'
number_of_sims = 10

run_dirs = []
for ii in range(number_of_sims):
    run_dirs.append('/DESRES-Trajectory_pnas2011a-C-'+str(ii)+'-no-water-no-lipid/pnas2011a-C-'+str(ii)+'-no-water-no-lipid')
dir_top = '/home/shenglan/topologies'
times_path = parent_dir+run_dirs[0]+'/'+run_dirs[0].split('/')[-1]+'_times.csv'

LOAD_STRIDE = 10
N_CLUSTER = 1300
PRIOR = 0.0#1.0/N_CLUSTER

#load list of mdtraj objects
simulations = []
for this_run_dir in run_dirs:
    simulations.append(util.load_trajs(this_run_dir,parent_dir,dir_top,load_stride = LOAD_STRIDE))

#track atoms in ligands
inds_N =[] #indices of N atom
atoms_to_track =['N']
ligands = ['ALP%d'%(ii+1) for ii in range(10)]


for ligand in ligands:
    iis = [atom.index for atom in simulations[0][0].topology.chain(0).atoms if atom.element.symbol in atoms_to_track
          and str(atom.residue)==ligand]
    inds_N.append(iis)

inds_all = [] #indices of all atoms, separated by ligands
for ligand in ligands:
    iis = [atom.index for atom in simulations[0][0].topology.chain(0).atoms if str(atom.residue) == ligand]
    inds_all.append(iis)

#track ASP 113
res_to_track = ['ASP113']
target_res = []
for residue in simulations[0][0].topology.chain(0).residues:
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

distances = []
for this_sim in simulations:
    for this_atom_pair in atom_pairs:
        this_lig = []
        for this_traj in this_sim: 
            this_lig.extend(md_dist.compute_distances(this_traj,[this_atom_pair]))
        distances.append(this_lig)

    


#sequences of coordinates of ligands
total_frames = 0
sequences_all = []
for this_sim in simulations:
    if use_COM:
        this_seq = util.featurize_RawPos(inds_all,this_sim,average = True)
    else:
        this_seq = util.featurize_RawPos(inds_N,this_sim)
    total_frames = this_seq[0].shape[0]*10+total_frames
    sequences_all.extend(this_seq)
print('the total number of conformations we are clustering is %d.' % total_frames)



time_step = util.calc_time_step(times_path,stride = LOAD_STRIDE)
 
clustering = KCenters(n_clusters = N_CLUSTER)
assignments = clustering.fit_predict(sequences_all)
centers = clustering.cluster_centers_
# 
# print len(assignments)
# print assignments[1].shape
# 
msm = MarkovStateModel(lag_time=80, reversible_type = 'transpose', 
ergodic_cutoff = 'off'
,verbose=True).fit(assignments)
countsmat = msm.countsmat_
transmat = msm.transmat_
print countsmat.shape

#try different lag_times
msmts0 = {}
msmts1 = {}
msmts2 = {}
lag_times = [1,5,10,20,30,40,50,60,70,80,90,100]#,150,200,250,300]
n_states = [N_CLUSTER]

for n in n_states:
    msmts0[n] = []
    msmts1[n] = []
    msmts2[n] = []

    for lag_time in lag_times:
        this_msm = MarkovStateModel(lag_time=lag_time, reversible_type = 'transpose', 
        ergodic_cutoff = 'off'
        , verbose=False).fit(assignments)
        timescales = this_msm.timescales_
        msmts0[n].append(timescales[0])
        msmts1[n].append(timescales[1])
        msmts2[n].append(timescales[2])
        print('n_states=%d\tlag_time=%.1f\ttimescales=%s (ns)' % (n, lag_time*time_step, np.array(timescales[0:3])*time_step))
    print('-------------------')


#----------------------------------------------------------------------------------
# Save data
#----------------------------------------------------------------------------------
msm_path = '/home/shenglan/TryMSMbuilder/output/C/all_clusters/KC_msm_c'+str(N_CLUSTER)+'_s'+str(LOAD_STRIDE)+'.out'
pickle.dump(msm,open(msm_path,'wb'))

seq_path = '/home/shenglan/TryMSMbuilder/output/C/all_clusters/sequences'+'_s'+str(LOAD_STRIDE)+'.out'
pickle.dump(sequences_all,open(seq_path,'wb'))

dist_path = '/home/shenglan/TryMSMbuilder/output/C/all_clusters/dist_to_binding'+'_s'+str(LOAD_STRIDE)+'.out'
pickle.dump(distances,open(dist_path,'wb'))

assign_path = '/home/shenglan/TryMSMbuilder/output/C/all_clusters/KC_assign_c'+str(N_CLUSTER)+'_s'+str(LOAD_STRIDE)+'.out'
pickle.dump(assignments,open(assign_path,'wb'))

countsmat_path = '/home/shenglan/TryMSMbuilder/output/C/all_clusters/KC_countsmat_c'+str(N_CLUSTER)+'_s'+str(LOAD_STRIDE)+'.out'
np.savetxt(countsmat_path,countsmat,fmt = '%8.4g')

transmat_path = '/home/shenglan/TryMSMbuilder/output/C/all_clusters/KC_transmat_c'+str(N_CLUSTER)+'_s'+str(LOAD_STRIDE)+'.out'
np.savetxt(transmat_path,transmat,fmt = '%10.4g')

centers_path = '/home/shenglan/TryMSMbuilder/output/C/all_clusters/KC_centers_c'+str(N_CLUSTER)+'_s'+str(LOAD_STRIDE)+'.out'
np.savetxt(centers_path,centers,fmt = '%10.4g')

#----------------------------------------------------------------------------------
# Convert to pdb
#----------------------------------------------------------------------------------
pdb_path = '/home/shenglan/TryMSMbuilder/output/C/all_clusters/KC_c'+str(N_CLUSTER)+'_s'+str(LOAD_STRIDE)+'.pdb'
util.convert_sequences_to_pdb(seq_path,assign_path,pdb_path)

#----------------------------------------------------------------------------------
# Visualization of data
#----------------------------------------------------------------------------------
# fig5 = plt.figure(figsize=(14,3))
# 
# plt.plot(assignments[0],'r.',alpha = 0.5)
# plt.plot(assignments[1],'b.',alpha = 0.5)
# plt.plot(assignments[2],'g.',alpha = 0.5)
# plt.xlabel('Frame')
# plt.ylabel('Cluster number')
# plt.title('cluster numbers, condition B')
# plt.savefig('/home/shenglan/TryMSMbuilder/output/fig5.png')
# plt.close(fig5)
# 
#----------------------------------------------------------------------------------
fig4 = plt.figure(figsize=(18,8))
plt.title('lag time vs relaxation time condition B')

for i, n in enumerate(n_states):
    plt.subplot(1,len(n_states),1+i)
    plt.plot(np.array(lag_times)*time_step, np.array(msmts0[n])*time_step)
    if i == 0:
        plt.ylabel('Relaxation Timescale 1(ns)')
    plt.xlabel('Lag Time (ns)')
    plt.title('%d states' % n)

plt.savefig('/home/shenglan/TryMSMbuilder/output/C/all_clusters/lagtime_ts1.png')
plt.close(fig4)

#----------------------------------------------------------------------------------
fig6 = plt.figure(figsize=(18,8))
plt.title('lag time vs relaxation time condition B')

for i, n in enumerate(n_states):
    plt.subplot(1,len(n_states),1+i)
    plt.plot(np.array(lag_times)*time_step, np.array(msmts1[n])*time_step)
    if i == 0:
        plt.ylabel('Relaxation Timescale 2(ns)')
    plt.xlabel('Lag Time (ns)')
    plt.title('%d states' % n)

plt.savefig('/home/shenglan/TryMSMbuilder/output/C/all_clusters/lagtime_ts2.png')
plt.close(fig6)

#----------------------------------------------------------------------------------
fig7 = plt.figure(figsize=(18,8))
plt.title('lag time vs relaxation time condition B')

for i, n in enumerate(n_states):
    plt.subplot(1,len(n_states),1+i)
    plt.plot(np.array(lag_times)*time_step, np.array(msmts2[n])*time_step)
    if i == 0:
        plt.ylabel('Relaxation Timescale 3(ns)')
    plt.xlabel('Lag Time (ns)')
    plt.title('%d states' % n)

plt.savefig('/home/shenglan/TryMSMbuilder/output/C/all_clusters/lagtime_ts3.png')
plt.close(fig7)

#----------------------------------------------------------------------------------
# fig1 = plt.figure(figsize = (100,100))
# 
# for this_seq in sequences_all:
#     for frame in this_seq:
#         plt.scatter(frame[0],frame[1],color = 'Blue', marker = '.')
#         plt.xlabel('x')
#         plt.ylabel('y')
# 
# for point in res_pos_ave:
#         plt.scatter(point[0],point[1],color = 'Green',marker = '*')
# plt.savefig('/home/shenglan/TryMSMbuilder/output/fig1.png')
# plt.close(fig1)


#----------------------------------------------------------------------------------
# 
# # fig2 = plt.figure()
# # for seq in zip(sequences_A_1,assignments):
# #     seq = zip(seq[0],seq[1])
# #     for point in seq:
# #         plt.scatter(point[0][0],point[0][2],color = colors[point[1]])
# #         plt.xlabel('x')
# #         plt.ylabel('z')
# # plt.scatter(res_pos_ave[0],res_pos_ave[2],color = 'Yellow', marker = '*',s= 500)
# # plt.savefig('/home/shenglan/TryMSMbuilder/output/fig2.png')
# # plt.close(fig2)
# # 
# # #----------------------------------------------------------------------------------
# # fig3 = plt.figure()
# # ax = fig3.add_subplot(111,projection = '3d')
# # ax.view_init(elev=30., azim=100)
# # 
# # for seq in zip(sequences_A_1,assignments):
# #     seq = zip(seq[0],seq[1])
# #     for point in seq:
# #         ax.scatter(point[0][0],point[0][1],point[0][2],color = colors[point[1]])
# #         ax.set_xlabel('x')
# #         ax.set_ylabel('y')
# #         ax.set_zlabel('z')
# # ax.scatter(res_pos_ave[0],res_pos_ave[1],res_pos_ave[2],color = 'Yellow', marker = '*',s=500)
# # plt.savefig('/home/shenglan/TryMSMbuilder/output/fig3.png')
# plt.close(fig3)

# colors = np.array([x for x in 'bgrcmykbgrcmykbgrcmykbgrcmyk'])
# colors = np.hstack([colors] * 20)

# #track residue Aps113
# first_res = str(simulations[0][0].topology.chain(0).residue(0))
# import re
# match = re.match(r"([a-z]+)([0-9]+)", first_res, re.I)
# if match:
#     items = match.groups()
# first_res_number = int(items[-1])
# target_res = simulations[0][0].topology.chain(0).residue(113-first_res_number)
# res_ind = [[atom.index for atom in target_res.atoms if atom.name in ['OD1']]]
