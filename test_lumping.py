from msmbuilder.tests import test_lumping as tlp
from msmbuilder.lumping import PCCA, PCCAPlus
import matplotlib.pyplot as plt
import numpy as np

from msmbuilder.msm import MarkovStateModel

assign, ref_macro_assign = tlp._metastable_system()

fig1 = plt.figure()
plt.plot(assign,'o')
plt.title('microstate assignments')
plt.savefig('/home/shenglan/TryMSMbuilder/output/test_microassign.png')
plt.close(fig1)

print assign.shape
print ref_macro_assign.shape

print assign[0]

test = PCCA(2).fit(assign)
macro_assign = test.fit_transform(assign)[0]
print test.n_states_
print test.mapping_

print macro_assign.shape
print np.sum(macro_assign)

print len(np.unique(macro_assign))==2

fig2 = plt.figure()
plt.plot(macro_assign,'o')
plt.title('macro assignments')
plt.savefig('/home/shenglan/TryMSMbuilder/output/test_macroassign.png')
plt.close(fig2)


macro_msm = MarkovStateModel().fit(macro_assign)
print macro_msm.n_states_
print macro_msm.mapping_
