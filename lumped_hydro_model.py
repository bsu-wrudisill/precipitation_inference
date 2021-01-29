#!/usr/bin/env python
# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import metpy.calc as mpcalc
from metpy.units import units
import sys


# #### Helper functions
class parms:
    frtdir = .09
    frtgw = .06
    smcap = 300
    ddpar = 1.
    etpar = .001
    tmelt = 1.0
    t_snow = 1.0
    t_melt = .8
    t_base = .9   # this for the ddpar calculation
    t_power = .09  # for ddpar calulation



def DayLength_Calc(dayOfYear, lat=40.0):
    """Computes the length of the day (the time between sunrise and
    sunset) given the day of the year and latitude of the location.
    Function uses the Brock model for the computations.
    For more information see, for example,
    Forsythe et al., "A model comparison for daylength as a
    function of latitude and day of year", Ecological Modelling,
    1995.
    Parameters
    ----------
    dayOfYear : int
        The day of the year. 1 corresponds to 1st of January
        and 365 to 31st December (on a non-leap year).
    lat : float
        Latitude of the location in degrees. Positive values
        for north and negative for south.
    Returns
    -------
    d : float
        Daylength in hours.
    """
    latInRad = np.deg2rad(lat)
    declinationOfEarth = 23.45*np.sin(np.deg2rad(360.0*(283.0+dayOfYear)/365.0))
    if -np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth)) <= -1.0:
        return 24.0
    elif -np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth)) >= 1.0:
        return 0.0
    else:
        hourAngle = np.rad2deg(np.arccos(-np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth))))
        return 2.0*hourAngle/15.0

    # https://gist.github.com/anttilipp/ed3ab35258c7636d87de6499475301ce




# Water  movement from top --> bottom layer
def Drainage(Wu, smcap):
    # D in the paper
    # Wu: water content of the upper layer
    # smcap: parameter
    if Wu > smcap:
        D = Wu - smcap
    else:
        D = 0
    # D is "drainage"
    return D


def PET_Calc(L, T, etpar):
    # Potential evapotranspiration calculator
    # L: Length of day
    # rho_v_sat: staturated absolute humidity
    # etpar: parameter
    temperature_c = T * units("degC")
    sat_vap_pres = mpcalc.saturation_vapor_pressure(temperature_c)

    # compute abs humidity
    abs_humidity = (2165 * sat_vap_pres.to('kPa').m)/temperature_c.to("K").m

    # compute PET from the equation in Clark
    return L*abs_humidity*etpar


###  Compute Snowmelt

def DDPar_Calc(T30, t_base, t_power):
    # T30: previous 30day mean temp
    # t_base: parameters
    # t_power: parameters

    k = 1.0 #mm/(d*C^2)
    if T30 > t_base:
        return (T30 - t_base)**t_power
    if T30 <= t_base:
        return t_base


def Ms(Td, Snow, tmelt, ddpar):
    # Td: temperature
    # Snow: amount in the bucket
    # ddpar: parameter (but not constant)
    # tmelt: parameter

    if Td > tmelt:
        melt = (Td-tmelt)*ddpar
        if melt > Snow: # check is there is enough snow to melt...
            melt = Snow
            return melt
        else:
            return melt
    else:
        return 0

###  Run the model...

snow_array_hide = []
def one_time_forward(Pd,               # Forcing (Precipitation)
                     Td,               # Forcing (Temperature)
                     T30,              # Forcing (30-Day Tempearture avg)
                     L,                # Forcing (Day Length)
                     Wu,               # Model State
                     Wb,               # Model State
                     Snow,             # Model State
                     frtdir,     # Parameter
                     frtgw,      # Parameter
                     smcap,      # Parameter
                     etpar,      # Parameter
                     tmelt,      # Parameter
                     t_snow,     # Parameter
                     t_melt,     # Parameter
                     t_base,     # Parameter
                     t_power):   # Parameter


    # Pd: Daily precipiation (mm/d)
    # Td: Daily temperature (T)
    # T30: Mean temperature from the previous 30 days
    # Wu: upper layer soil moisture (mm)
    # Wb: bottom layer soil moisture(mm)
    # Snow: Current snow bucket (mm)
    # L: length of day (hours)
    # parameters....

    # Determine precipitation phase
    # Add snow to snowpack if necessary
    if Td >= t_snow:
        Pr = Pd
    else:
        Pr = 0                   # no rain

    # compute daily snow accumulation
    Ps = Pd - Pr  # snow = total - rain

    # Compute the snowmelt
    # compute ddpar first
    ddpar = DDPar_Calc(T30, t_base, t_power)
    Md = Ms(Td, Snow, t_melt, ddpar)

    # Compute the direct runoff (Qd)
    Qd = frtdir * (Pr + Md)

    # Compute the Soil Water input
    Wi = (1-frtdir) * (Pr + Md) # Wi in paper

    # Compute drainage; 1st layer of soil --> 2nd layer of soil
    D = Drainage(Wu, smcap)

    # Compute baseflow (bottom layer to runoff)
    Qb = Wb * frtgw

    # Compute PET
    # Firt compute saturation specific humidity
    #rho_v_sat = .0001 # use mpcaclc to compute this
    PET = PET_Calc(L, Td, etpar)
#    PET = 8.0 # mm...

    # Compute AET
    E = PET*(Wu/smcap)

    # Compute change in the water balance of top layer
    dWudt = Wi - D - E

    # Compute change in water balance at bottom layer
    dWbdt = D - Qb

    # update values
    # soil moisture
    Wb = Wb + dWbdt
    Wu = Wu + dWudt

    # snow
    Snow = Snow + Ps #- Md # last snow, new snow, snow melt

    # Compute discharge
    Q = Qb + Qd

    # now return things..
    return Snow, Wb, Wu, Qb, Qd, E



def ForwardModel(DailyTemp,
                 DailyPrecip,
                 LenOfDayHr,
                 dz,
                 lapse_rate = -.0065,
                 orog_gradient = .002, # Parameter
                 M = 1.0,             # Parameter
                 frtdir = .0001,     # Parameter
                 frtgw = .06,      # Parameter
                 smcap = 100,      # Parameter
                 etpar = .005,     # Parameter
                 tmelt = .01,      # Parameter
                 t_snow = 0.0,     # Parameter
                 t_melt = 1.0,      # Parameter
                 t_base = 0,      # Parameter
                 t_power = .5):   # Parameter

    # Run the forward model for the given forcing data
    assert len(DailyTemp) == len(DailyPrecip) == len(LenOfDayHr)
    ntimes = len(DailyTemp)

    # Create the output/input arrays
    Wb = np.zeros(ntimes, dtype='float')
    Wu = np.zeros(ntimes, dtype='float')
    Snow = np.zeros(ntimes, dtype='float')
    Qb = np.zeros(ntimes, dtype='float')
    Qd = np.zeros(ntimes, dtype='float')
    ET = np.zeros(ntimes, dtype='float')

    # Loop through ntimes
    for t in range(1,ntimes):
        L = LenOfDayHr[t]
        # adjust precipitation
        Pd = M*np.mean(DailyPrecip[t] + DailyPrecip[t]*dz*orog_gradient)
        Td = np.mean(DailyTemp[t] + dz*lapse_rate)

        # get rid of nan values...
        # if np.nan in [Pd, Td]:
        #     Q_array[t] = Q_array[t-1]
        #     Wb[t] = Wb[t-1]
        #     Wu[t] = Wu[t-1]
        #     Snow[t] = Snow[t-1]
        #     break

        # compute T30; ignore temps from before the starting point...
        T30 = np.mean(DailyTemp[np.max([0, t-30]):t])

        # Run the model one timespep forward... save result
        Snow_t, Wb_t, Wu_t, Qb_t, Qd_t, ET_t = one_time_forward(Pd,
                                                                Td,
                                                                T30,
                                                                L,
                                                                Wu[t-1],
                                                                Wb[t-1],
                                                                Snow[t-1],
                                                                frtdir,
                                                                frtgw,
                                                                smcap,
                                                                etpar,
                                                                tmelt,
                                                                t_snow,
                                                                t_melt,
                                                                t_base,
                                                                t_power)

        # store the times
        Qb[t] = Qb_t
        Qd[t] = Qd_t
        Wb[t] = Wb_t
        Wu[t] = Wu_t
        Snow[t] = Snow_t
        ET[t] = ET_t

    # compute the total discharge and return it
    Q = Qb + Qd
    return Q



if __name__ == '__main__':
    # Import stuff to read the forcings
    from ForcingModel import *

    # Gather the forcings
    # start = "2010-10-01"
    # end = "2011-09-30"
    # f = ForcingModel()
    # f.read()
#    daily_temp, daily_precip, daily_q_observed, day_len_hrs = f(start, end)
    # dz = f.hypsometry - f.snotel_elev

    # np.save("daily_temp.npy", daily_temp)
    # np.save("daily_precip.npy", daily_precip)
    # np.save("day_len_hrs.npy", day_len_hrs)
    # np.save("dz.npy", dz)


    daily_temp = np.load("daily_temp.npy")
    daily_precip = np.load("daily_precip.npy")
    daily_q_observed = np.load("daily_q_observed.npy")
    day_len_hrs = np.load("day_len_hrs.npy")
    dz = np.load("dz.npy")
    Q = ForwardModel(daily_temp, daily_precip, day_len_hrs, dz, M=1.0)
    plt.plot(Q)
    plt.show()





