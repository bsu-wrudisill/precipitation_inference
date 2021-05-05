import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import metpy.calc as mpcalc
from metpy.units import units
import sys
from ForcingModel import *


opg = .01
lapse_rate = -.0065
mult = 1.1


daily_temp = np.load("daily_temp.npy")
daily_precip = np.load("daily_precip.npy")
dz = np.load("dz_reduced.npy")

idz = len(dz)

newtd = np.zeros((365))
newpd = np.zeros((365))
newpd2 = np.zeros((365))


for t in range(365):
    dTd = 0
    dPd = 0
    for z in range(idz):
        dTd = dTd + dz[z] * -.0065
        dPd = dPd + daily_precip[t] * mult * dz[z] * opg

    dTd = dTd/idz
    dPd = dPd/idz

    newtd[t] = daily_temp[t] + dTd
    newpd[t] = daily_precip[t] + dPd


    ### Second method -- vecotr math. same as for loop, confirmed
    # dpvec =  np.array([daily_precip[t]]*idz).reshape(idz,1)
    # dzr = dz.reshape(1,idz)
    # dPd2 = np.dot(dzr*opg*mult, dpvec)/idz
    # newpd2[t] = daily_precip[t] + dPd2



plt.plot(daily_precip, label='base')
plt.plot(newpd, label='adj')
plt.plot(newpd2, label='adj2')

plt.show()