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

sequences = pickle.load(open('/home/shenglan/TryMSMbuilder/output/sequences.out','rb'))
print len(sequences)
print sequences[0].shape

KC_assign = pickle.load(open('/home/shenglan/TryMSMbuilder/output/KC_assign_50.out','rb'))
outfile = open('/home/shenglan/TryMSMbuilder/output/rec_crds.pdb','wb')
outfile.write('%s     %s \n' % ('TITLE', 'cluster centers'))
outfile.write('%s        %d \n' % ('MODEL' , 1))
atom_count = 0
for this_seq,this_assign in zip(sequences,KC_assign):
    atom_num = this_seq.shape[0]
    
    for ii in range(atom_num):
        atom_count = atom_count + 1
        outfile.write("%6s%5d %4s %3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n" % 
        ('HETATM', ii+1,'C','CLU','I',this_assign[ii]+1,'D'
        ,this_seq[ii,0]*10,this_seq[ii,1]*10,this_seq[ii,2]*10, 1.0, 0.0,'C')) 
        #the factor of 10 multiplied to the coordinate is a conversion from nm to angstrom, to match the .mae file provided by DesRes
        
outfile.write('%s \n' % 'ENDMDL')
outfile.write('%s    \n' % 'END')
    
print('done converting!')