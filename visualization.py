import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmap
import pickle

# transmat = np.loadtxt('/home/shenglan/TryMSMbuilder/output/transmat.out')
# countsmat = np.loadtxt('/home/shenglan/TryMSMbuilder/output/countsmat.out')
# 
# fig1 = plt.figure()
# plt.imshow(transmat,cmap = cmap.Blues)
# plt.title('transition matrix')
# plt.colorbar()
# plt.savefig('/home/shenglan/TryMSMbuilder/output/transmat.png')
# plt.close(fig1)
# 
# fig2 = plt.figure()
# plt.imshow(countsmat,cmap = cmap.BuGn)
# plt.title('count matrix')
# plt.colorbar()
# plt.savefig('/home/shenglan/TryMSMbuilder/output/countsmat.png')
# plt.close(fig2)
# 
# assign = pickle.load(open('/home/shenglan/TryMSMbuilder/output/KC_assign_5_sNone.out','rb'))
# fig3 = plt.figure()
# print len(assign[0])
# plt.title("one ligand's cluster assignment")
# plt.plot(assign[0],'b.',alpha = 0.5)
# plt.plot(assign[1],'g.',alpha = 0.5)
# plt.plot(assign[3],'y.',alpha = 0.5)
# plt.ylim(-0.5,4)
# plt.savefig('/home/shenglan/TryMSMbuilder/output/assignment.png')
# plt.close(fig3)

transC_path = '/home/shenglan/TryMSMbuilder/output/C/KC_transmat_c500_s10.out'
countsC_path = '/home/shenglan/TryMSMbuilder/output/C/KC_countsmat_c500_s10.out'
transmat_C = np.loadtxt(transC_path)
countsmat_C = np.loadtxt(countsC_path)

fig4 = plt.figure()
plt.imshow(transmat_C,cmap = cmap.Blues)
plt.title('transition matrix')
plt.colorbar()
plt.savefig(transC_path.split('.')[0]+'.png')
plt.close(fig4)

fig5 = plt.figure()
plt.imshow(countsmat_C,cmap = cmap.BuGn)
plt.title('count matrix')
plt.colorbar()
plt.savefig(countsC_path.split('.')[0]+'.png')
plt.close(fig5)

