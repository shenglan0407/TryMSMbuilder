import pickle
import numpy as np
import matplotlib.pyplot as plt

model = 'c2000_s10.out'
micro_msm_path = '/home/shenglan/TryMSMbuilder/output/C/KC_msm_'+model
micro_msm = pickle.load(open(micro_msm_path,'rb'))

timescales = micro_msm.timescales_
# print timescales.shape
# print timescales[1290:]

eigenvalues = micro_msm.eigenvalues_
right_eigenvectors = micro_msm.right_eigenvectors_

test_vector = right_eigenvectors[:,-1]
print test_vector[0]
# print right_eigenvectors.shape
# print eigenvalues[1290:]
# fig1 = plt.figure()
# plt.plot(timescales[30:],'o')
# plt.title('1000th to 1307th timescales')
# plt.savefig('/home/shenglan/TryMSMbuilder/output/timescales.png')
# plt.close(fig1)