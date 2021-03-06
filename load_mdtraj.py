import utilities as util
import numpy as np
import os

import mdtraj as md

use_COM = False
parent_dir = '/scratch/PI/rondror/DesRes-Simulations/B2AR_PNAS_2011a_Dror_Desres/DESRES-Trajectory_pnas2011a-C-no-water-no-lipid.1'
number_of_sims = 1


run_dirs = []
for ii in range(number_of_sims):
    run_dirs.append('/DESRES-Trajectory_pnas2011a-C-'+str(ii)+'-no-water-no-lipid/pnas2011a-C-'+str(ii)+'-no-water-no-lipid')
dir_top = '/home/shenglan/topologies'
times_path = parent_dir+run_dirs[0]+'/'+run_dirs[0].split('/')[-1]+'_times.csv'

files = []
for this_file in os.listdir(parent_dir+run_dirs[0]):
    
    if this_file.endswith('.dcd'):
        files.append(this_file)
files = sorted(files)

LOAD_STRIDE = 10000

test_run = run_dirs[0]
test_traj = util.load_trajs(test_run,parent_dir,dir_top,load_stride = LOAD_STRIDE)

inds_N =[] #indices of N atom
atoms_to_track =['N']
ligands = ['ALP%d'%(ii+1) for ii in range(10)]
print ligands

print test_traj[0].topology.chain(0).atom(4667).residue
for ligand in ligands:
    iis = [atom.index for atom in test_traj[0].topology.chain(0).atoms if atom.element.symbol in atoms_to_track
          and str(atom.residue)==ligand]
    inds_N.append(iis)

# inds_N = [[atom.index] for atom in test_traj[0].topology.chain(0).atoms if atom.element.symbol in atoms_to_track
#            and atom.residue.name== 'ALP']
print inds_N

test_seq = util.featurize_RawPos(inds_N,test_traj)
print test_seq[0][0]
print test_seq[1][0]

for atom_idx in inds_N: 
    print(tuple(test_traj[0].xyz[0, atom_idx,:]))

print(tuple(test_traj[0].xyz[0, 4979,:]))

print test_traj[0]
# 
# test_traj2 = md.load_dcd(parent_dir + run_dirs[0] + '/' + files[0]
#                                     , top = '/home/shenglan/topologies/test.pdb'
#                                     , stride = LOAD_STRIDE
#                                     )
# 
# test_seq2 = util.featurize_RawPos(inds_N,[test_traj2])
# print test_seq2[0][0]
# print test_seq2[1][0]
# 
# inds_N2 = []
# for ligand in ligands:
#     iis = [atom.index for atom in test_traj2.topology.chain(1).atoms if atom.element.symbol in atoms_to_track
#           and atom.residue == ligand]
#     inds_N2.append(iis)
# 
# print inds_N2
# 
# 
# for atom_idx in inds_N2: 
#     print(tuple(test_traj2.xyz[0, atom_idx,:]))
# 
# print(tuple(test_traj2.xyz[0, 223,:]))