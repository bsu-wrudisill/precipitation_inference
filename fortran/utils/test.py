import numpy as np
import sys
from myfastmodule import fast_foward_model as ffm
from matplotlib import pyplot as plt
import matplotlib as mpl
mpl.use("TkAgg")

daily_temp = np.load("../data/daily_temp.npy")
daily_precip = np.load("../data/daily_precip.npy")
daily_q_observed = np.load("../data/daily_q_observed.npy")
day_len_hrs = np.load("../data/day_len_hrs.npy")
dz = np.load("../data/dz_reduced.npy")

PdVec = daily_precip
TdVec = daily_temp
LVec = day_len_hrs
dz = dz

# for i in range(1000):
frtdir = .09 #+ .02*np.random.randn(1)
frtgw  = .06 #+ .02*np.random.randn(1)
smcap = 200 #+ 50.0*np.random.randn(1)      # Parameter
etpar = max(.001, 0)#.009+ .01*np.random.randn(1))    # Parameter
tmelt = .01 #+ .1*np.random.randn(1)      # Parameter
t_snow = 0.0# + .1*np.random.randn(1)     # Parameter
t_melt = 1.0# + .1*np.random.randn(1)      # Parameter
t_base = 1.2# + .1*np.random.randn(1)     # Parameter
t_power = .5#  + .01*np.random.randn(1) # Parame
opg  = .005 #+ .001*np.random.randn(1)
bias = 0 #+ 5.*np.random.randn(1)


Q, E, PET, Snow, PdVecOut = ffm.fwd(3, PdVec, TdVec, LVec, dz, frtdir, frtgw, smcap, etpar, t_snow, t_melt, t_base, t_power, bias, opg)
print(PdVecOut.sum(0).sum())


Q, E, PET, Snow, PdVecOut = ffm.fwd(2, PdVec, TdVec, LVec, dz, frtdir, frtgw, smcap, etpar, t_snow, t_melt, t_base, t_power, bias, opg)
print(PdVecOut.sum(0).sum())


print(PdVec.sum())

# plt.plot(Q)
# plt.show()
# plt.plot(PdVec, label='snotel')
# plt.legend()
# plt.show()
# ettot = E.sum()
# ptot = PdVecOut.sum()
# pettot = PET.sum()

# x = pettot/ptot
# y = ettot/ptot
# plt.scatter(x, y)
# plt.plot([0,1], [0,1])
# plt.plot([1,4], [1,1])

# # # analytical budyko
# # def bdk(pet_p_ratio, n=1.5):
# # 	return pet_p_ratio/((pet_p_ratio)**n + 1.0)**(1/n)
# np.save('FakeQ', Q)

# plt.plot(model, color='blue', alpha=.1)
# plt.plot(E)
# #plt.plot(daily_q_observed, label='obs', color='black')
# plt.legend()
# plt.show()