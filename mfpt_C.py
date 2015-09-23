from msmbuilder.tpt import mfpt
from msmbuilder.cluster import KCenters

import pickle
import numpy as np
import copy


LOAD_STRIDE = 10
N_MICRO = 1300
BULK_CUTOFF = 2.0 #unit is nm

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

# Identify 'sinks' for mfpt calculations
bulk_clusters = []
for ii in clusters_dist.keys():
    if np.average(clusters_dist[ii]) > BULK_CUTOFF:
        bulk_clusters.append(ii)
print len(bulk_clusters)
    
mfpt_mat = mfpt.mfpts(micro_msm,sinks = bulk_clusters,lag_time = 80*1.8)
print mfpt_mat.shape
print mfpt_mat[mfpt_mat > 1000.0]
print np.where(mfpt_mat> 1000.0)
