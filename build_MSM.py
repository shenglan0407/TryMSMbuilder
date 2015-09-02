import utilities as util
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from msmbuilder.cluster import KCenters


#global parameters
parent_dir = '/scratch/PI/rondror/DesRes-Simulations/B2AR_PNAS_2011a_Dror_Desres/DESRES-Trajectory_pnas2011a-A-no-water-no-lipid.1'
A_run_dir = '/DESRES-Trajectory_pnas2011a-A-1-no-water-no-lipid/pnas2011a-A-1-no-water-no-lipid'
dir_top = '/home/shenglan/topologies'
times_path = parent_dir+A_run_dir+'/'+A_run_dir.split('/')[-1]+'_times.csv'

LOAD_STRIDE = 1000

colors = np.array([x for x in 'bgrcmykbgrcmykbgrcmykbgrcmyk'])
colors = np.hstack([colors] * 20)

#load list of mdtraj objects

trajs_A_1 = util.load_trajs(A_run_dir,parent_dir,dir_top,load_stride = LOAD_STRIDE)

#track atoms in ligands
inds =[] #indices of Cs in the aromatic ring
atoms_in_ring = ['C8','C9','C10','C11','C12','C13']
ligands = [residue for residue in trajs_A_1[0].topology.chain(1).residues][:-1]
for ligand in ligands:
    iis = [atom.index for atom in trajs_A_1[0].topology.chain(1).atoms if atom.name in atoms_in_ring
          and atom.residue == ligand]
    inds.append(iis)

#track residue D113
first_res = str(trajs_A_1[0].topology.chain(0).residue(0))
import re
match = re.match(r"([a-z]+)([0-9]+)", first_res, re.I)
if match:
    items = match.groups()
first_res_number = int(items[-1])
target_res = trajs_A_1[0].topology.chain(0).residue(113-first_res_number)
res_ind = [atom.index for atom in target_res.atoms if atom.name in ['OD1']]

sequences_A_1 = util.featurize_RawPos(inds,trajs_A_1)
res_pos_A_1 = util.featurize_RawPos([res_ind],trajs_A_1)

res_pos_ave = np.mean(res_pos_A_1[0],axis = 0)
print res_pos_ave

fig1 = plt.figure()

for this_seq in sequences_A_1:
    for frame in this_seq:
        plt.scatter(frame[0],frame[1])
        plt.xlabel('x')
        plt.ylabel('y')
for this_seq in res_pos_A_1:
    for frame in this_seq:
        plt.scatter(frame[0],frame[1],color = 'Green')
plt.savefig('/home/shenglan/TryMSMbuilder/output/fig1.png')
plt.close(fig1)


time_step = util.calc_time_step(times_path,stride = LOAD_STRIDE)

print time_step

clustering = KCenters(n_clusters = 4)
assignments = clustering.fit_predict(sequences_A_1)
centers = clustering.cluster_centers_

fig2 = plt.figure()
for seq in zip(sequences_A_1,assignments):
    seq = zip(seq[0],seq[1])
    for point in seq:
        plt.scatter(point[0][0],point[0][2],color = colors[point[1]])
        plt.xlabel('x')
        plt.ylabel('z')
plt.scatter(res_pos_ave[0],res_pos_ave[2],color = 'Yellow', marker = '*',s= 500)
plt.savefig('/home/shenglan/TryMSMbuilder/output/fig2.png')
plt.close(fig2)

fig3 = plt.figure()
ax = fig3.add_subplot(111,projection = '3d')
ax.view_init(elev=30., azim=100)

for seq in zip(sequences_A_1,assignments):
    seq = zip(seq[0],seq[1])
    for point in seq:
        ax.scatter(point[0][0],point[0][1],point[0][2],color = colors[point[1]])
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
ax.scatter(res_pos_ave[0],res_pos_ave[1],res_pos_ave[2],color = 'Yellow', marker = '*',s=500)
plt.savefig('/home/shenglan/TryMSMbuilder/output/fig3.png')
plt.close(fig3)