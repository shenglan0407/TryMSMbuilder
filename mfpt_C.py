from msmbuilder.tpt import mfpt
from msmbuilder.cluster import KCenters

import pickle
import numpy as np
import copy


LOAD_STRIDE = 10
N_MICRO = 1300
BULK_CUTOFF = 6.75 #unit is nm

model = 'c'+str(N_MICRO)+'_s'+str(LOAD_STRIDE)
micro_msm_path = '/home/shenglan/TryMSMbuilder/output/C/all_clusters/KC_msm_'+model+'.out'
micro_msm = pickle.load(open(micro_msm_path,'rb'))

micro_assign_path = '/home/shenglan/TryMSMbuilder/output/C/all_clusters/KC_microassign_c'+str(N_MICRO)+'_s'+str(LOAD_STRIDE)+'.out'
micro_assign = pickle.load(open(micro_assign_path,'rb'))

#distances
dist_path = '/home/shenglan/TryMSMbuilder/output/C/all_clusters/dist_to_binding'+'_s'+str(LOAD_STRIDE)+'.out'
distances = pickle.load(open(dist_path,'rb'))

# dictionary; key: microstate id; 
# item: list of distances to binding pocket for points in microstate
clusters_dist = {}

for ii in range(len(distances)):
    if len(distances) == len(micro_assign) and len(distances[ii]) == len(micro_assign[ii]):
        for nn in range(len(distances[ii])):
            if micro_assign[ii][nn] in clusters_dist.keys():
                clusters_dist[micro_assign[ii][nn]].append(distances[ii][nn])
            else:
                clusters_dist.update({micro_assign[ii][nn] : [distances[ii][nn]]})
                
    else:
        print ('something is wrong. distances and micro_assign not matching in dimensions')
        print len(distances[ii]) 
        print len(micro_assign[ii])

# Identify 'sinks' for mfpt_off calculations
bulk_clusters = []
bound_clusters = [288,1020]
for ii in clusters_dist.keys():
    if np.average(clusters_dist[ii]) > BULK_CUTOFF:
        bulk_clusters.append(ii)
   #  else:
#         bound_clusters.append(ii)
print len(bulk_clusters)
    
mfpt_off = mfpt.mfpts(micro_msm,sinks = bulk_clusters,lag_time = 80*1.8)
mfpt_on = mfpt.mfpts(micro_msm,sinks = bound_clusters,lag_time = 80*1.8)
#print mfpt_mat.shape
print "Here are so mfpt_off times (ns) and corresponding cluster ids:"
print mfpt_off[mfpt_off > 7.354e6]
print np.where(mfpt_off> 7.354e6)
print ("Maximum mfpt_off time is %.3g ns in cluster %d." % (np.max(mfpt_off),np.argmax(mfpt_off)))
print ("and the mean mfpt_off is %.3g ns"%np.mean(mfpt_off))

print ("Maximum mfpt_on time is %.g ns in cluster %d." % (np.max(mfpt_on),np.argmax(mfpt_on)))
print ("and the mean mfpt_on is %.3g ns"%np.mean(mfpt_on))

mfpt_on_path = '/home/shenglan/TryMSMbuilder/output/C/all_clusters/KC_mfpton_c'+str(N_MICRO)+'_s'+str(LOAD_STRIDE)+'.out'
mfpt_off_path = '/home/shenglan/TryMSMbuilder/output/C/all_clusters/KC_mfptoff_c'+str(N_MICRO)+'_s'+str(LOAD_STRIDE)+'.out'

pickle.dump(mfpt_on,open(mfpt_on_path,'wb'))
pickle.dump(mfpt_off,open(mfpt_off_path,'wb'))

volume = ()