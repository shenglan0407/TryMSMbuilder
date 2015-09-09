import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmap

transmat = np.loadtxt('/home/shenglan/TryMSMbuilder/output/transmat.out')
countsmat = np.loadtxt('/home/shenglan/TryMSMbuilder/output/countsmat.out')

fig1 = plt.figure()
plt.imshow(transmat,cmap = cmap.Blues)
plt.title('transition matrix')
plt.colorbar()
plt.savefig('/home/shenglan/TryMSMbuilder/output/transmat.png')
plt.close(fig1)

fig2 = plt.figure()
plt.imshow(countsmat,cmap = cmap.BuGn)
plt.title('count matrix')
plt.colorbar()
plt.savefig('/home/shenglan/TryMSMbuilder/output/countsmat.png')
plt.close(fig2)
