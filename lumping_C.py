#---------------------------------------------------------------------
# Author: Shenglan Qiao <shenglan@stanford.edu>
# Contributors:
# Copyright (c) 2015, Stanford University
# All rights reserved.

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import numpy as np
import pickle
import matplotlib.pyplot as plt

from msmbuilder.lumping import PCCA, PCCAPlus
from msmbuilder.msm import MarkovStateModel

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------
model = 'c2000_s10.out'
micro_msm_path = '/home/shenglan/TryMSMbuilder/output/C/KC_msm_'+model
micro_msm = pickle.load(open(micro_msm_path,'rb'))

geo_assign_path = '/home/shenglan/TryMSMbuilder/output/C/KC_assign_'+model
geo_assign = pickle.load(open(geo_assign_path,'rb'))
#print geo_assign.shape
print geo_assign[0][0]
print('There are %d microstates in msm' % micro_msm.n_states_)

raw_clusters = []
for this_assign in geo_assign:
    raw_clusters.extend(np.unique(this_assign))
raw_clusters = np.unique(np.array(raw_clusters))
print('There are %d clusters in the original geometric clustering.'%len(raw_clusters))


raw_to_micro_mapping = micro_msm.mapping_

dummy_assign = 1308
for ii in range(len(raw_clusters)):
    if ii in raw_to_micro_mapping.keys():
        pass
    else:
        raw_to_micro_mapping[ii] = dummy_assign
        dummy_assign = dummy_assign+1
print raw_to_micro_mapping[geo_assign[0][0]]
print geo_assign[0][0]

micro_assign = []

for this_assign in geo_assign:
    this_list = np.zeros(len(this_assign),dtype = int)
    micro_assign.append(this_list)
    
    
for nn in range(len(micro_assign)):
    for ii in range(len(micro_assign[nn])):
        micro_assign[nn][ii] = raw_to_micro_mapping[geo_assign[nn][ii]]

#check for mismatch
for this_assign, that_assign in zip(micro_assign,geo_assign):
    for this_item,that_item in zip(this_assign,that_assign):
        if this_item == raw_to_micro_mapping[that_item]:
            pass
        else:
            print ('mismatch %d, %d' % (this_item, that_item))
        
# pcca = PCCA.from_msm(micro_msm,2)
# macro_assign = pcca.fit_transform(geo_assign)[0]
# 
# print len(macro_assign)
