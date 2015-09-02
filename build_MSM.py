import utilities as util
import numpy as np
import matplotlib.pyplot as plt


parent_dir = '/scratch/PI/rondror/DesRes-Simulations/B2AR_PNAS_2011a_Dror_Desres/DESRES-Trajectory_pnas2011a-A-no-water-no-lipid.1'
A_run_dir = '/DESRES-Trajectory_pnas2011a-A-1-no-water-no-lipid/pnas2011a-A-1-no-water-no-lipid'
dir_top = '/home/shenglan/topologies'
times_path = parent_dir+A_run_dir+'/'+A_run_dir.split('/')[-1]+'_times.csv'

LOAD_STRIDE = 2000

trajs_A_1 = util.load_trajs(A_run_dir,parent_dir,dir_top,load_stride = LOAD_STRIDE)

inds =[] #indices of Cs in the aromatic ring
atoms_in_ring = ['C8','C9','C10','C11','C12','C13']
ligands = [residue for residue in trajs_A_1[0].topology.chain(1).residues][:-1]
for ligand in ligands:
    iis = [atom.index for atom in trajs_A_1[0].topology.chain(1).atoms if atom.name in atoms_in_ring
          and atom.residue == ligand]
    inds.append(iis)

print inds

sequences_A_1 = util.featurize_RawPos(inds,trajs_A_1)

fig1 = plt.figure()

for this_seq in sequences_A_1:
    for frame in this_seq:
        plt.scatter(frame[0],frame[1])
plt.savefig('/home/shenglan/TryMSMbuilder/output/fig1.png')
plt.close(fig1)


time_step = util.calc_time_step(times_path,stride = LOAD_STRIDE)

print time_step