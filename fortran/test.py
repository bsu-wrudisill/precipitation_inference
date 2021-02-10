import numpy as np
import sys
from myfastmodule import fast_foward_model as ffm
from matplotlib import pyplot as plt

daily_temp = np.load("../data/daily_temp.npy")
daily_precip = np.load("../data/daily_precip.npy")
daily_q_observed = np.load("../data/daily_q_observed.npy")
day_len_hrs = np.load("../data/day_len_hrs.npy")
dz = np.load("../data/dz_reduced.npy")

PdVec = daily_precip
TdVec = daily_temp
LVec = day_len_hrs
dz = dz


for i in range(100):
	frtdir = .09 + .02*np.random.randn(1)
	frtgw  = .06 + .02*np.random.randn(1)
	smcap = 200 + 50.0*np.random.randn(1)      # Parameter
	etpar = .09 + .1*np.random.randn(1)    # Parameter
	tmelt = .01 + .1*np.random.randn(1)      # Parameter
	t_snow = 0.0 + .1*np.random.randn(1)     # Parameter
	t_melt = 1.0 + .1*np.random.randn(1)      # Parameter
	t_base = 1.2 + .1*np.random.randn(1)     # Parameter
	t_power = .5  + .01*np.random.randn(1) # Parame
	opg  = .002 + .001*np.random.randn(1)
	bias = 0 + 5.*np.random.randn(1)


	model = ffm.fwd(3, PdVec, TdVec, LVec, dz, frtdir, frtgw, smcap, etpar, t_snow, t_melt, t_base, t_power, bias, opg)
	plt.plot(model, color='blue', alpha=.1)

plt.plot(daily_q_observed, label='obs', color='black')
plt.legend()
plt.show()