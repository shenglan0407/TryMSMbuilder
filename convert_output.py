import utilities as util
import numpy as np
import pickle

# crd_file = '/home/shenglan/TryMSMbuilder/output/cluster_centers.out'
# pdb_path = '/home/shenglan/TryMSMbuilder/output/cluster_centers.pdb'
# 
# util.convert_to_pdb(crd_file,pdb_path)
# 
# crd_file = '/home/shenglan/TryMSMbuilder/output/KM_centers_50.out'
# pdb_path = '/home/shenglan/TryMSMbuilder/output/KM_centers_50.pdb'
# 
# util.convert_to_pdb(crd_file,pdb_path)
# 
# crd_file = '/home/shenglan/TryMSMbuilder/output/KC_centers_50.out'
# pdb_path = '/home/shenglan/TryMSMbuilder/output/KC_centers_50.pdb'
# 
# util.convert_to_pdb(crd_file,pdb_path)


seqfile = '/home/shenglan/TryMSMbuilder/output/sequences.out'
assignfile = '/home/shenglan/TryMSMbuilder/output/KC_assign_5.out'
savefile = '/home/shenglan/TryMSMbuilder/output/KC_rec_crds.pdb'
util.convert_sequences_to_pdb(seqfile,assignfile,savefile)