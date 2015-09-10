import utilities as util
import numpy as np

crd_file = '/home/shenglan/TryMSMbuilder/output/cluster_centers.out'
pdb_path = '/home/shenglan/TryMSMbuilder/output/cluster_centers.pdb'

util.convert_to_pdb(crd_file,pdb_path)