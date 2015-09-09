import numpy as np

crds = np.loadtxt('/home/shenglan/TryMSMbuilder/output/cluster_centers.out')
outfile = open('/home/shenglan/TryMSMbuilder/output/cluster_centers.pdb','wb')
atom_num = crds.shape[0]

outfile.write('%s     %s \n' % ('TITLE', 'cluster centers'))
outfile.write('%s        %d \n' % ('MODEL' , 1))
for ii in range(atom_num):
    outfile.write("%6s%5.2d %4s %s %s  %s   %7.4f %7.4f %7.4f %s %5s %12s  \n" % ('HETATM', ii,'C','CLU','I',str(ii)+'D'
    ,crds[ii,0]*10,crds[ii,1]*10,crds[ii,2]*10,'1.00','0.00','C'))
    
outfile.write('%s \n' % 'ENDMDL')
outfile.write('%s    \n' % 'END')