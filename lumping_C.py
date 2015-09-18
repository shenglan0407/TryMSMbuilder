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
print micro_assign.shape

print('there are %d microstates in msm' % micro_msm.n_states_)

#pcca = PCCA.from_msm(micro_msm,50)
#macro_assign = pcca.fit_transform(micro_assign)[0]

print micro_assign.shape
