from myfastmodule import fast_foward_model as ffm
import numpy as np
import time
import matplotlib.pyplot as plt

frtdir = .0001,     # Parameter
frtgw = .06,      # Parameter
smcap = 100,      # Parameter
etpar = .005,     # Parameter
tmelt = .01,      # Parameter
t_snow = 0.0,     # Parameter
t_melt = 1.0,      # Parameter
t_base = 0,      # Parameter
t_power = .5   # Parameter


N = 365

TdVec = np.load("../daily_temp.npy")
PdVec = np.load("../daily_precip.npy")
daily_q_observed = np.load("../daily_q_observed.npy")
LVec = np.load("../day_len_hrs.npy")


#ffm.fwd

# YOU DON'T NEED TO PROVIDE THE DIMENSION ARGUMENTS!!!
start = time.time()
a = ffm.fwd(PdVec, TdVec, LVec, frtdir, frtgw, smcap, etpar, t_snow, t_melt, t_base, t_power)
print('It took', time.time()-start, 'seconds.')


plt.plot(a, label='model')
plt.plot(daily_q_observed, label='obs')
plt.xlabel("DOY")
plt.ylabel("Runoff (mm)")
plt.legend()
plt.show()