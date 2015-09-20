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

micro_assign_path = '/home/shenglan/TryMSMbuilder/output/C/KC_assign_'+model
micro_assign = pickle.load(open(micro_assign_path,'rb'))
#print micro_assign.shape

print('There are %d microstates in msm' % micro_msm.n_states_)

raw_clusters = []
for this_assign in micro_assign:
    raw_clusters.extend(np.unique(this_assign))
raw_clusters = np.unique(np.array(raw_clusters))
print('There are %d clusters in the original geometric clustering.'%len(raw_clusters))
print np.min(raw_clusters)
print np.max(raw_clusters)

raw_to_micro_mapping = micro_msm.mapping_
# 
# pcca = PCCA.from_msm(micro_msm,2)
# macro_assign = pcca.fit_transform(micro_assign)[0]
# 
# print len(macro_assign)
