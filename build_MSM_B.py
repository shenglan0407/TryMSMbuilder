#---------------------------------------------------------------------

import utilities as util
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

from msmbuilder.cluster import KCenters
from msmbuilder.msm import MarkovStateModel

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------

#global parameters
use_COM = False
parent_dir = '/scratch/PI/rondror/DesRes-Simulations/B2AR_PNAS_2011a_Dror_Desres/DESRES-Trajectory_pnas2011a-B-no-water-no-lipid.1'
number_of_sims = 20

run_dirs = []
for ii in range(number_of_sims):
    run_dirs.append('/DESRES-Trajectory_pnas2011a-B-'+str(ii)+'-no-water-no-lipid/pnas2011a-B-'+str(ii)+'-no-water-no-lipid')
dir_top = '/home/shenglan/topologies'
times_path = parent_dir+run_dirs[0]+'/'+run_dirs[0].split('/')[-1]+'_times.csv'

LOAD_STRIDE = None

colors = np.array([x for x in 'bgrcmykbgrcmykbgrcmykbgrcmyk'])
colors = np.hstack([colors] * 20)

#load list of mdtraj objects
simulations = []
for this_run_dir in run_dirs:
    simulations.append(util.load_trajs(this_run_dir,parent_dir,dir_top,load_stride = LOAD_STRIDE))

#track atoms in ligands
inds_N =[] #indices of N atom
atoms_to_track =['N']
ligands = [residue for residue in simulations[0][0].topology.chain(1).residues][:-1]

for ligand in ligands:
    iis = [atom.index for atom in simulations[0][0].topology.chain(1).atoms if atom.element.symbol in atoms_to_track
          and atom.residue == ligand]
    inds_N.append(iis)
    
inds_all = [] #indices of all atoms, separated by ligands
for ligand in ligands:
    iis = [atom.index for atom in simulations[0][0].topology.chain(1).atoms if atom.residue == ligand]
    inds_all.append(iis)


#track residue Aps113
first_res = str(simulations[0][0].topology.chain(0).residue(0))
import re
match = re.match(r"([a-z]+)([0-9]+)", first_res, re.I)
if match:
    items = match.groups()
first_res_number = int(items[-1])
target_res = simulations[0][0].topology.chain(0).residue(113-first_res_number)
res_ind = [[atom.index for atom in target_res.atoms if atom.name in ['OD1']]]

#sequences of coordinates of ligand aromatic ring and Aps113
sequences_all = []
for this_sim in simulations:
    if use_COM:
        this_seq = util.featurize_RawPos(inds_all,this_sim,average = True)
    else:
        this_seq = util.featurize_RawPos(inds_N,this_sim)
    sequences_all.extend(this_seq)

#print len(sequences_all)
#print sequences_all[-1].shape


res_pos = []
res_pos_ave =[]
total_frames = 0
for this_sim in simulations:
    this_seq = util.featurize_RawPos(res_ind,this_sim)
    res_pos.extend(this_seq)
    total_frames = total_frames + len(this_seq[0])
    res_pos_ave.append(np.mean(this_seq[0],axis=0))

print total_frames


time_step = util.calc_time_step(times_path,stride = LOAD_STRIDE)
 
clustering = KCenters(n_clusters = 2000)
assignments = clustering.fit_predict(sequences_all)
centers = clustering.cluster_centers_

#print len(assignments)
#print assignments[1].shape

msm = MarkovStateModel(lag_time=180, verbose=True).fit(assignments)
countsmat = msm.countsmat_
transmat = msm.transmat_
#print np.sum(countsmat)

#np.savetxt('/home/shenglan/TryMSMbuilder/output/assignments.out',assignments, fmt = '%3.0f')
np.savetxt('/home/shenglan/TryMSMbuilder/output/countsmat.out',countsmat,fmt = '%8.4g')
np.savetxt('/home/shenglan/TryMSMbuilder/output/transmat.out',transmat,fmt = '%10.4g')
np.savetxt('/home/shenglan/TryMSMbuilder/output/cluster_centers.out',centers,fmt = '%10.4g')
np.savetxt('/home/shenglan/TryMSMbuilder/output/res_pos_ave.out',res_pos_ave,fmt = '%10.4g')


#try different lag_times
msmts0 = {}
lag_times = [20,60,100,150,180]
n_states = [2000,2200,2400,2600]

for n in n_states:
    msmts0[n] = []
    for lag_time in lag_times:
        assignments = KCenters(n_clusters=n).fit_predict(sequences_all)
        msm = MarkovStateModel(lag_time=lag_time, verbose=False).fit(assignments)
        timescales = msm.timescales_
        msmts0[n].append(timescales[0])
        print('n_states=%d\tlag_time=%d\ttimescales=%s' % (n, lag_time, timescales[0:2]))
    print('-------------------')


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

#----------------------------------------------------------------------------------
fig4 = plt.figure(figsize=(18,5))
plt.title('lag time vs relaxation time condition B')

for i, n in enumerate(n_states):
    plt.subplot(1,len(n_states),1+i)
    plt.plot(np.array(lag_times)*time_step, np.array(msmts0[n])*time_step)
    if i == 0:
        plt.ylabel('Relaxation Timescale (ns)')
    plt.xlabel('Lag Time (ns)')
    plt.title('%d states' % n)

plt.savefig('/home/shenglan/TryMSMbuilder/output/fig4.png')
plt.close(fig4)

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