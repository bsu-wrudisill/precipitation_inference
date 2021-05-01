import sys
sys.path.append("/Users/williamrudisill/Documents/Conceptual_Runoff_Model/fortran/")
sys.path.append("/Users/williamrudisill/Documents/Conceptual_Runoff_Model/")

from snowmodule17 import snowmodule17 as sm17
from ForcingModel import ForcingModel

import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import pandas as pd
from aso_process import *
from snow17 import *
mpl.use('Qt5Agg')
import pathlib
import xarray as xr
import scipy.optimize


##################################################################################################################
########                          Calibrate Snow17 against Snotel SWE (only)                              ########
########                                                                                                  ########
##################################################################################################################

start_date = "2010-10-01"
end_date = "2019-09-30"

data_base = "/Users/williamrudisill/Documents/Conceptual_Runoff_Model/data/"

# load data ...
frc = ForcingModel()
frc.read()

# read in forcings for the desired time period
daily_temp, daily_precip, q, day_len_hrs, daily_swe = frc(start_date, end_date)

# read these in (they are fixed in time)
elev = np.load(data_base + "elev.npy")
dz = np.load(data_base + "dz.npy")

# process dates
dates = pd.date_range(start_date, end_date, freq='D')
jdays = dates.dayofyear.values * 1.0



def snow17_caller_01(jdays, daily_precip, daily_temp, x):
  outflow,swe,rain,ptot = sm17.snow17driver(jdays,
                            daily_precip,
                            daily_temp,
                            nlayers=1,
                            dz=0,
                            dt=24,
                            rvs=1,
                            opg_method=1,
                            opg=.0021,
                            bias=x[0],
                            uadj=x[1],
                            mbase=x[2],
                            mfmax=x[3],
                            mfmin=x[4],
                            tipm=x[5],
                            nmf=x[6],
                            plwhc=x[7],
                            pxtemp=x[8],
                            pxtemp1=x[9],
                            pxtemp2=x[10])
  return outflow,swe,rain,ptot

# this packs in all of the forcing data.. so now the function is just of x
snow17_caller = lambda x: snow17_caller_01(jdays, daily_precip, daily_temp, x)

model_parameters = {"bias": 0.1,
                    "uadj":  .04,
                    "mbase": .3,
                    "mfmax": 4.0,
                    "mfmin": 1.0,
                    "tipm": .1,
                    "nmf": .40,
                    "plwhc": .04,
                    "pxtemp": 2.0,
                    "pxtemp1": -1.0,
                    "pxtemp2": 3.0}

x = np.fromiter(model_parameters.values(), dtype='float')


outflow, swe, rain,ptot = snow17_caller(x)
#plt.plot(m0.cumsum(), label='melt')
plt.plot(outflow.cumsum(), label='outflow')
plt.plot(rain.cumsum(), label='rain')
plt.plot(ptot.cumsum(), label='ptot')
plt.plot(swe[0], label='snow')
plt.legend()
plt.show()




