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


# data base
# start_date = "2017-10-01"
# end_date = "2018-09-30"

start_date = "2000-10-01"
end_date = "2020-09-30"

data_base = "/Users/williamrudisill/Documents/Conceptual_Runoff_Model/data/"

# load data ...
frc = ForcingModel()
frc.read()


# read these in (they are fixed in time)
elev = np.load(data_base + "elev.npy")
dz = np.load(data_base + "dz.npy")



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





def objective_function(mod_swe, snotel_swe, alpha=.3, beta=.8):
  # compute the error in the model based on snotel
  # and ASO overflights ....
  aso_error = np.sqrt(np.mean((mod_swe[1] - snotel_swe)**2))
  return alpha*aso_error #+ (beta/2)*(f1_error + f2_error)


def _objective_function_wrapper(x, jdays, daily_precip, daily_temp, daily_swe, fx):
  swe,m = fx(jdays, daily_precip, daily_temp, x)
  obj = objective_function(swe, daily_swe)
  return obj


def snow17_caller(jdays, daily_precip, daily_temp, x):
  swe,m = sm17.snow17driver(jdays,
                            daily_precip,
                            daily_temp,
                            nlayers=3,
                            dz=heights-3096.768,
                            dt=24,
                            rvs=1,
                            opg_method=0,
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

model_parameters = {"bias": .5,
                    "uadj":  .04,
                    "mbase": 1.0,
                    "mfmax": 1.9,
                    "mfmin": .6,
                    "tipm": .1,
                    "nmf": .40,
                    "plwhc": .04,
                    "pxtemp": 2.0,
                    "pxtemp1": -1.0,
                    "pxtemp2": 3.0}

def cross_validate(calibrate, validate, model_parameters=model_parameters):
    """Summary

    Args:
        calibrate (TYPE): Description
        validate (TYPE): Description
        model_parameters (TYPE, optional): Description
    """

    # RUN THE CALIBRATION PERIOD
    start_date, end_date = calibrate
    # read in forcings for the desired time period
    daily_temp, daily_precip, q, day_len_hrs, daily_swe = frc(start_date, end_date)
    dates = pd.date_range(start_date, end_date, freq='D')
    jdays = dates.dayofyear.values * 1.0
    # get the initial conditions
    x0 = np.fromiter(model_parameters.values(), dtype='float')

    # we don't need to calibrate the undercatch
    snow17objective_fun = lambda x: _objective_function_wrapper(x, jdays, daily_precip, daily_temp, daily_swe, snow17_caller)
    result = scipy.optimize.minimize(snow17objective_fun, x0, method='Nelder-Mead', options={"maxiter":5000})

    swe0,m0 = snow17_caller(jdays, daily_precip, daily_temp, x0)
    swe,m  = snow17_caller(jdays, daily_precip, daily_temp, result.x)

    # Compute/recompute the errors for intiial and calibrated model
    cal_init_error = objective_function(swe0,daily_swe)
    cal_post_error = objective_function(swe,daily_swe)


    # Run the VALIDATION PERIOD
    val_init_error_list = []
    val_post_error_list = []

    for d in validate:
        start_date, end_date = d
        daily_temp, daily_precip, q, day_len_hrs, daily_swe = frc(start_date, end_date)
        dates = pd.date_range(start_date, end_date, freq='D')
        jdays = dates.dayofyear.values * 1.0

        # make the function again
        snow17objective_fun = lambda x: _objective_function_wrapper(x, jdays, daily_precip, daily_temp, daily_swe, snow17_caller)
        result = scipy.optimize.minimize(snow17objective_fun, x0, method='Nelder-Mead', options={"maxiter":5000})

        swe0,m0 = snow17_caller(jdays, daily_precip, daily_temp, x0)
        swe,m  = snow17_caller(jdays, daily_precip, daily_temp, result.x)

        # Compute/recompute the errors for intiial and calibrated model
        val_init_error = objective_function(swe0,daily_swe)
        val_post_error = objective_function(swe,daily_swe)

        # append them
        val_init_error_list.append(val_init_error)
        val_post_error_list.append(val_post_error)

    # now finish
    return val_init_error_list, val_post_error_list, cal_init_error, cal_post_error


if __name__ == '__main__':
    calibration_period = ("2000-10-01", "2010-10-01")
    validation_period = [("2017-10-01", "2018-10-01"),
                         ("2018-10-01", "2019-10-01"),
                         ("2017-10-01", "2018-10-01"),
                         ("2017-10-01", "2018-10-01"),
                         ("2017-10-01", "2018-10-01"),
                         ("2017-10-01", "2018-10-01"),
                         ("2017-10-01", "2018-10-01"),
                         ("2017-10-01", "2018-10-01")]


    val_init_error_list, val_post_error_list, cal_init_error, cal_post_error = cross_validate(calibration_period, validation_period, model_parameters=model_parameters)



