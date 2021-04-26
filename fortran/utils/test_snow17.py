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

start_date = "1990-10-01"
end_date = "2020-09-30"

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


def objective_function(snotel_swe, fx, x):
  # compute the error in the model based on snotel
  # and ASO overflights ....
  swe,m = fx(x)
  snotel_error = np.sqrt(np.mean((swe[0] - snotel_swe)**2))
  return snotel_error


def snow17_caller_01(jdays, daily_precip, daily_temp, x):
  swe,m = sm17.snow17driver(jdays,
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
  return swe,m

# this packs in all of the forcing data.. so now the function is just of x
snow17_caller_a = lambda x: snow17_caller_01(jdays, daily_precip, daily_temp, x)

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

x0 = np.fromiter(model_parameters.values(), dtype='float')

# we don't need to calibrate the undercatch
snow17objective_fun = lambda x: objective_function(daily_swe, snow17_caller_a, x)

# Using the ...Powell algorithm seems to produce really weird values (lots of negatives)

result1 = scipy.optimize.minimize(snow17objective_fun, x0, method='Nelder-Mead', options={"disp":True})
# result2 = scipy.optimize.minimize(snow17objective_fun, x0, method='Nelder-Mead', options={"maxiter":5000, "disp":True})

swe0,m0 = snow17_caller_a(x0)
swe1,m2 = snow17_caller_a(result1.x)
#swe2,m2 = snow17_caller(jdays, daily_precip, daily_temp, result2.x)

# save the parameters for the next calibration step ...
x1 = result1.x
# plt.plot(dates, daily_swe.values, linestyle='--', color='black', label='Butte Snotel')
# plt.plot(dates, swe1[0], color='purple', label="snow17-calibrated")#, label="calibrated -- Powel")
# sys.exit()


# # plt.plot(dates, swe2[0], color='orange', label="calibrated -- NM")
# # plt.plot(dates, swe0[0], color='blue', alpha=.75, label="uncalibrated")
# plt.legend()
# plt.savefig("snow17-snotel-calibrated")
# sys.exit()


##################################################################################################################
########                          Calibrate AGainst ASO Obs                                               ########
########                                                                                                  ########
##################################################################################################################
# Get forcing data for new dates
start_date = "2017-10-01"
end_date = "2019-09-30"
daily_temp, daily_precip, q, day_len_hrs, daily_swe = frc(start_date, end_date)
dates = pd.date_range(start_date, end_date, freq='D')
jdays = dates.dayofyear.values * 1.0



### Do the elevation computing
nlayers = 3
nelev = int(elev.shape[0])
intrv = int(nelev/nlayers)
#intrvD = int(np.argwhere(elev>3700)[0])

# make an array of the fraction of each bin
areas = np.array([intrv,
                  intrv,
                  intrv])/nelev

# array of the mean elevations of each bin
heights = np.array([np.mean(elev[0:intrv]),
                    np.mean(elev[intrv:intrv*2]),
                    np.mean(elev[intrv*2:]),
                    ])



aso_date01 = "2018-03-31"  # "2019-04-07"
aso_date02 = "2018-05-24"  # "2019-06-10"

aso_date03 = "2019-04-07"
aso_date04 = "2019-06-10"

# get aso data for first flight
aso_swe_01 = np.array([get_aso_snow(aso_date01, elev[0], elev[intrv]),
                       get_aso_snow(aso_date01, elev[intrv], elev[intrv*2]),
                       get_aso_snow(aso_date01, elev[intrv*2], 5000.)])

# aso data for 2nd flight
aso_swe_02 = np.array([get_aso_snow(aso_date02, elev[0], elev[intrv]),
                       get_aso_snow(aso_date02, elev[intrv], elev[intrv*2]),
                       get_aso_snow(aso_date02, elev[intrv*2], 5000.)])

#aso data for 3rd flight
aso_swe_03 = np.array([get_aso_snow(aso_date03, elev[0], elev[intrv]),
                       get_aso_snow(aso_date03, elev[intrv], elev[intrv*2]),
                       get_aso_snow(aso_date03, elev[intrv*2], 5000.)])

# aso data for 4th flight
aso_swe_04 = np.array([get_aso_snow(aso_date04, elev[0], elev[intrv]),
                       get_aso_snow(aso_date04, elev[intrv], elev[intrv*2]),
                       get_aso_snow(aso_date04, elev[intrv*2], 5000.)])



# indices for the flight dates
d1ind = (pd.to_datetime(aso_date01) - pd.to_datetime(start_date)).days
d2ind = (pd.to_datetime(aso_date02) - pd.to_datetime(start_date)).days
d3ind = (pd.to_datetime(aso_date03) - pd.to_datetime(start_date)).days
d4ind = (pd.to_datetime(aso_date04) - pd.to_datetime(start_date)).days






def snow17_caller_02(jdays, daily_precip, daily_temp, x0, x):
  swe,m = sm17.snow17driver(jdays,
                            daily_precip,
                            daily_temp,
                            nlayers=3,
                            dz=heights-3096.768,
                            dt=24,
                            rvs=1,
                            opg_method=1,
                            opg=x[0],
                            bias=x[1],
                            uadj=x0[1],
                            mbase=x0[2],
                            mfmax=x0[3],
                            mfmin=x0[4],
                            tipm=x0[5],
                            nmf=x0[6],
                            plwhc=x0[7],
                            pxtemp=x0[8],
                            pxtemp1=x0[9],
                            pxtemp2=x0[10])
  return swe,m


snow17_caller_b = lambda x: snow17_caller_02(jdays, daily_precip, daily_temp, x1, x)

def objective_function(snotel_swe, fx, x, alpha=.1, beta=.9):

  mod_swe,melt= fx(x)
  # compute the error in the model based on snotel
  # and ASO overflights ....
  #snotel_error = np.sqrt(np.mean((mod_swe[1] - snotel_swe)**2))

  # sum up errors for first flight. note indices are off since there is an extra layer for the snotel...
  f1_error = np.sqrt(np.mean((mod_swe[0][d1ind] - aso_swe_01[0])**2 +
                             (mod_swe[1][d1ind] - aso_swe_01[1])**2 +
                             (mod_swe[2][d1ind] - aso_swe_01[2])**2))

  f2_error = np.sqrt(np.mean((mod_swe[0][d2ind] - aso_swe_02[0])**2 +
                             (mod_swe[1][d2ind] - aso_swe_02[1])**2 +
                             (mod_swe[2][d2ind] - aso_swe_02[2])**2))

  f3_error = np.sqrt(np.mean((mod_swe[0][d3ind] - aso_swe_03[0])**2 +
                             (mod_swe[1][d3ind] - aso_swe_03[1])**2 +
                             (mod_swe[2][d3ind] - aso_swe_03[2])**2))


  f4_error = np.sqrt(np.mean((mod_swe[0][d4ind] - aso_swe_04[0])**2 +
                             (mod_swe[1][d4ind] - aso_swe_04[1])**2 +
                             (mod_swe[2][d4ind] - aso_swe_04[2])**2))


  # return alpha*snotel_error + (beta/4)*(f1_error + f2_error + f3_error + f4_error)
  return (beta/4)*(f1_error + f2_error + f3_error + f4_error)


# do a 2nd calibration w/ aso data
#snow17_caller_02x = lambda jdays, daily_precip, daily_temp, x: snow17_caller_02(jdays, daily_precip, daily_temp, x0, x1)

snow17objective_fun = lambda x: objective_function(daily_swe, snow17_caller_b, x)
result = scipy.optimize.minimize(snow17objective_fun, [.002, 0.], method='Powell', options={"maxiter":5000, "disp":True})


swe,m = snow17_caller_b(result.x)

# plot the results
plt.plot(dates, swe[0], label=heights[0], color='gray', alpha=.5)
plt.plot(dates, swe[1], label=heights[1], color='orange', alpha=.5)
plt.plot(dates, swe[2], label=heights[2], color='blue', alpha=.5)


#plot the 1st swe lidar flight
plt.scatter([pd.to_datetime(aso_date01)]*nlayers,
            aso_swe_01,
            c=['gray', 'orange', 'blue'],
            marker='x')


# plot the 2nd swe lidar flight
plt.scatter([pd.to_datetime(aso_date02)]*nlayers,
            aso_swe_02,
            marker='x',
            c=['gray', 'orange', 'blue'])


# plot the 2nd swe lidar flight
plt.scatter([pd.to_datetime(aso_date03)]*nlayers,
            aso_swe_03,
            marker='x',
            c=['gray', 'orange', 'blue'])

# plot the 2nd swe lidar flight
plt.scatter([pd.to_datetime(aso_date04)]*nlayers,
            aso_swe_04,
            marker='x',
            c=['gray', 'orange', 'blue'])



# SAVE THE FINAL PARAMETERS
plt.legend()
plt.show()
# Get the ASO mean swe for the different bands


