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
import copy

import utilities as util

from msmbuilder.lumping import PCCA, PCCAPlus
from msmbuilder.msm import MarkovStateModel

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------
LOAD_STRIDE = 10
model = 'c1300_s'+str(LOAD_STRIDE)
micro_msm_path = '/home/shenglan/TryMSMbuilder/output/C/all_clusters/KC_msm_'+model+'.out'
micro_msm = pickle.load(open(micro_msm_path,'rb'))

geo_assign_path = '/home/shenglan/TryMSMbuilder/output/C/all_clusters/KC_assign_'+model+'.out'
geo_assign = pickle.load(open(geo_assign_path,'rb'))
#print geo_assign.shape
print('There are %d microstates in msm' % micro_msm.n_states_)

raw_clusters = []
for this_assign in geo_assign:
    raw_clusters.extend(np.unique(this_assign))
raw_clusters = np.unique(np.array(raw_clusters))
print('There are %d clusters in the original geometric clustering.'%len(raw_clusters))

#generate a mapping of geometric clustering to microstates
raw_to_micro_mapping = copy.copy(micro_msm.mapping_)
N_MICRO = len(raw_to_micro_mapping)
dummy_assign = N_MICRO
for ii in range(len(raw_clusters)):
    if ii in raw_to_micro_mapping.keys():
        pass
    else:
        raw_to_micro_mapping[ii] = dummy_assign
        #dummy_assign = dummy_assign+1
        
micro_assign = []

for this_assign in geo_assign:
    this_list = np.zeros(len(this_assign),dtype = int)
    micro_assign.append(this_list)
    
#map assignments for microstates    
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

partial_raw_to_micro_mapping = copy.copy(micro_msm.mapping_)
partial_micro_assign = []
for this_assign in geo_assign:
    this_list = []
    for this_item in this_assign:
        try:
            this_list.append(partial_raw_to_micro_mapping[this_item])
        except KeyError:
            pass
    partial_micro_assign.append(this_list)
    
unique_assign = []
for this_assign in partial_micro_assign:
    unique_assign.extend(np.unique(this_assign))
unique_assign = np.unique(np.array(unique_assign))

N_MACRO = 5
pcca = PCCA.from_msm(micro_msm,N_MACRO)
micro_to_macro_mapping = {}
for ii in range(len(pcca.microstate_mapping_)):
    micro_to_macro_mapping[ii] = pcca.microstate_mapping_[ii]

# for n macrostates, any frame that does not belong to any microstate and thus macrostate is labeled to be in the nth macrostate
for ii in range(len(raw_clusters)):
    if ii in micro_to_macro_mapping.keys():
        pass
    else:
        micro_to_macro_mapping[ii] = N_MACRO
    
# assignments to macro states
macro_assign = []

for this_assign in geo_assign:
    this_list = np.zeros(len(this_assign),dtype = int)
    macro_assign.append(this_list)

for nn in range(len(macro_assign)):
    for ii in range(len(macro_assign[nn])):
        macro_assign[nn][ii] = micro_to_macro_mapping[micro_assign[nn][ii]]

macro_assign_path = '/home/shenglan/TryMSMbuilder/output/C/all_clusters/KC_macroassign_c'+str(N_MACRO)+'_s'+str(LOAD_STRIDE)+'.out'
pickle.dump(macro_assign,open(macro_assign_path,'wb'))

micro_assign_path = '/home/shenglan/TryMSMbuilder/output/C/all_clusters/KC_microassign_c'+str(N_MICRO)+'_s'+str(LOAD_STRIDE)+'.out'
pickle.dump(micro_assign,open(micro_assign_path,'wb'))

seq_path = '/home/shenglan/TryMSMbuilder/output/C/all_clusters/sequences_s'+str(LOAD_STRIDE)+'.out'

save_macro_pdb_path = '/home/shenglan/TryMSMbuilder/output/C/all_clusters/KC_macroassign_c'+str(N_MACRO)+'_s'+str(LOAD_STRIDE)+'.pdb'
save_micro_pdb_path = '/home/shenglan/TryMSMbuilder/output/C/all_clusters/KC_microassign_c'+str(N_MICRO)+'_s'+str(LOAD_STRIDE)+'.pdb'

util.convert_sequences_to_pdb(seq_path,macro_assign_path,save_macro_pdb_path)
util.convert_sequences_to_pdb(seq_path,micro_assign_path,save_micro_pdb_path)


# macro_assign = pcca.fit(partial_micro_assign)
# 
# print macro_assign
