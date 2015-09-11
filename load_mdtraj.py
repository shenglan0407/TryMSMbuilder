import utilities as util
import numpy as np
import os

use_COM = False
parent_dir = '/scratch/PI/rondror/DesRes-Simulations/B2AR_PNAS_2011a_Dror_Desres/DESRES-Trajectory_pnas2011a-B-no-water-no-lipid.1'
number_of_sims = 20



run_dirs = []
for ii in range(number_of_sims):
    run_dirs.append('/DESRES-Trajectory_pnas2011a-B-'+str(ii)+'-no-water-no-lipid/pnas2011a-B-'+str(ii)+'-no-water-no-lipid')
dir_top = '/home/shenglan/topologies'
times_path = parent_dir+run_dirs[0]+'/'+run_dirs[0].split('/')[-1]+'_times.csv'

# files = []
# for this_file in os.listdir(parent_dir+run_dirs[2]):
#     
#     if this_file.endswith('.dcd'):
#         files.append(this_file)
#         
# for this_file in sorted(files):
#         print this_file

LOAD_STRIDE = None

test_run = run_dirs[0]
test_traj = util.load_trajs(test_run,parent_dir,dir_top,load_stride = LOAD_STRIDE)

inds_N =[] #indices of N atom
atoms_to_track =['N']
ligands = [residue for residue in test_traj[0].topology.chain(1).residues][:-1]

for ligand in ligands:
    iis = [atom.index for atom in test_traj[0].topology.chain(1).atoms if atom.element.symbol in atoms_to_track
          and atom.residue == ligand]
    inds_N.append(iis)

print inds_N

frame_ind = 10000
test_seq = util.featurize_RawPos(inds_N,test_traj)
print test_seq[0][10000]
print test_seq[1][10000]

for atom_idx in inds_N: 
    print(tuple(test_traj[1].xyz[0, atom_idx,:]))