import pickle
import numpy as np
import matplotlib.pyplot as plt

model = 'c1000_sNone.out'
micro_msm_path = '/home/shenglan/TryMSMbuilder/output/C/all_clusters/KC_msm_'+model
micro_msm = pickle.load(open(micro_msm_path,'rb'))

countsmat = micro_msm.countsmat_
# for ii in range(countsmat.shape[0]):
#     if np.sum(countsmat[ii][:]) == 0.0:
#         print ii
#         print np.sum(countsmat[ii][:])
rev_counts = (countsmat+countsmat.T)*0.5

for ii in range(rev_counts.shape[0]):
    if np.sum(rev_counts[ii][:]) == 0.0:
        print ii
        print np.sum(rev_counts[ii][:])

transmat = micro_msm.transmat_

print countsmat[(countsmat < 1.0) & (countsmat > 0.0)]
if True in np.isnan(transmat):
    print('there is nan')
# timescales = micro_msm.timescales_
# # print timescales.shape
# # print timescales[1290:]
# 
# eigenvalues = micro_msm.eigenvalues_
# right_eigenvectors = micro_msm.right_eigenvectors_
# 
# test_vector = right_eigenvectors[:,-1]
# print test_vector[0]
# print right_eigenvectors.shape
# print eigenvalues[1290:]
# fig1 = plt.figure()
# plt.plot(timescales[30:],'o')
# plt.title('1000th to 1307th timescales')
# plt.savefig('/home/shenglan/TryMSMbuilder/output/timescales.png')
# plt.close(fig1)