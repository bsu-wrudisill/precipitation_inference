#!/usr/bin/env python
"""Summary
"""
# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import metpy.calc as mpcalc
from metpy.units import units
import sys
import time


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
    t_base = 1.1   # this for the ddpar calculation
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
#@jit
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

#@jit
def PET_CalcFaster(L, t, etpar):
    # Potential evapotranspiration calculator
    # L: Length of day
    # rho_v_sat: staturated absolute humidity
    # etpar: parameter

    sat_vap_pres = (.6112)*np.exp((17.67*t)/(t + 243.5))

    # compute abs humidity
    abs_humidity = (2165 * sat_vap_pres)/(t + 273.15)

    # compute PET from the equation in Clark
    return L*abs_humidity*etpar


# def PET_Calc(L, T, etpar):
#     # Potential evapotranspiration calculator
#     # L: Length of day
#     # rho_v_sat: staturated absolute humidity
#     # etpar: parameter
#     temperature_c = T * units("degC")
#     sat_vap_pres = mpcalc.saturation_vapor_pressure(temperature_c)

#     # compute abs humidity
#     abs_humidity = (2165 * sat_vap_pres.to('kPa').m)/temperature_c.to("K").m

#     # compute PET from the equation in Clark
#     return L*abs_humidity*etpar


###  Compute Snowmelt
#@jit
def DDPar_Calc(T30, t_base, t_power):
    # T30: previous 30day mean temp
    # t_base: parameters
    # t_power: parameters

    k = 1.0 #mm/(d*C^2)
    if T30 > t_base:
        return (T30 - t_base)**t_power
    if T30 <= t_base:
        return t_base

#@jit
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

#
def one_time_forward(nLayers,         # number of layers for the snow model
                     Pd,               # Forcing (Precipitation)
                     Td,               # Forcing (Temperature)
                     T30,              # Forcing (30-Day Tempearture avg)
                     L,                # Forcing (Day Length)
                     Wu,               # Model State
                     Wb,               # Model State
                     Snow,             # Model State
                     frtdir,     # Parameter
                     frtgw,      # Parameter
                     smcap,      # Parameter
                     etpar,     # Parameter
                     tmelt,      # Parameter
                     t_snow,     # Parameter
                     t_melt,      # Parameter
                     t_base,      # Parameter
                     t_power,
                     ):   # Parameter
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

    # Compute the snowmelt and accumulation for each layer

    # Melt and rain accumulator variables
    Mdsum = 0
    Prsum = 0
    Tdavg = 0
    for k in range(nLayers):
        # Get the precipitation and temperature
        if Td[k] >= t_snow:
            Prk = Pd[k]
        else:
            Prk = 0                   # no rain

        # Compute daily snow accumulation
        Psk = Pd[k] - Prk

        # Compute the snowmelt
        # compute ddpar first
        print(T30[k], t_base, t_power)
        ddpar = DDPar_Calc(T30[k], t_base, t_power)
        Mdk = Ms(Td[k], Snow[k], t_melt, ddpar)

        # gather up the total melt
        Mdsum = Mdk + Mdsum
        Prsum = Prsum + Prk

        # melt/accumulate the snowpack
        Snow[k] = Snow[k] + Psk - Mdk # last snow, new snow, snow melt

        Tdavg = Tdavg + Td[k]

    Tdavg = Tdavg/nLayers
    # --------------------------------------------
    # Now compute the Soil runoff as one layer ...
    # --------------------------------------------


    # Compute the direct runoff (Qd)
    Qd = frtdir * (Prsum + Mdsum)

    # Compute the Soil Water input
    Wi = (1-frtdir) * (Prsum + Mdsum) # Wi in paper

    # Compute drainage; 1st layer of soil --> 2nd layer of soil
    D = Drainage(Wu, smcap)

    # Compute baseflow (bottom layer to runoff)
    Qb = Wb * frtgw

    # Compute PET
    # Firt compute saturation specific humidity
    #rho_v_sat = .0001 # use mpcaclc to compute this
    PET = PET_CalcFaster(L, Tdavg, etpar)
    # PET = 8.0 # mm...

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


    # Compute discharge
    Q = Qb + Qd

    # now return things..
    return Snow, Wb, Wu, Qb, Qd, E




def ForwardModelFasterLayers(DailyTemp,
                             DailyPrecip,
                             LenOfDayHr,
                             dz,
                             nLayers = 3,
                             lapse_rate = -.0065,
                             orog_gradient = .002, # Parameter
                             M = .30,             # Parameter
                             frtdir = .0001,     # Parameter
                             frtgw = .02,      # Parameter
                             smcap = 120,      # Parameter
                             etpar = .07,     # Parameter
                             tmelt = .01,      # Parameter
                             t_snow = 0.8,     # Parameter
                             t_melt = 0.1,      # Parameter
                             t_base = .9,      # Parameter
                             t_power = .9):   # Parameter

    # Run the forward model for the given forcing data
    #assert len(DailyTemp) == len(DailyPrecip) == len(LenOfDayHr)
#    ntimes = len(DailyTemp)
    ntimes = 365
    Wb = 20.0
    Wu = 20.0
    Qb = 0.0
    Qd = 0.0
    ET = 0.0
    Snow = np.array([0.0, 0.0, 0.0])
    Q = np.zeros(ntimes, dtype='float')
    Snowholder = np.zeros((3, ntimes))
    # Divide region into "nLayers" EQUAL segments -- maybe change this
    seg1 = int(len(dz)/nLayers)
    seg2 = seg1*2
    ldz = len(dz)

    tholder = np.zeros((3, ntimes))

    # Loop through ntimes
    for t in range(1,ntimes):
        L = LenOfDayHr[t]

        Pd = np.zeros(nLayers)
        Td = np.zeros(nLayers)
        T30 = np.zeros(nLayers)


        Pd[0] = M*np.mean(DailyPrecip[t] + DailyPrecip[t]*dz[0:seg1]*orog_gradient)
        Pd[1] = M*np.mean(DailyPrecip[t] + DailyPrecip[t]*dz[seg1:seg2]*orog_gradient)
        Pd[2] = M*np.mean(DailyPrecip[t] + DailyPrecip[t]*dz[seg2:ldz]*orog_gradient)


        Td[0] = np.mean(DailyTemp[t] + dz[0:seg1]*lapse_rate)
        Td[1] = np.mean(DailyTemp[t] + dz[seg1:seg2]*lapse_rate)
        Td[2] = np.mean(DailyTemp[t] + dz[seg2:]*lapse_rate)

        tholder[:,t] = Td

        # compute T30; ignore temps from before the starting point...
        # T30 = np.mean(DailyTemp[np.max([0, t-30]):t])

        T30[0] = np.mean([DailyTemp[x] + dz[:seg1]*lapse_rate for x in range(np.max([1, t-30]))])
        T30[1] = np.mean([DailyTemp[x] + dz[seg1:seg2]*lapse_rate for x in range(np.max([1, t-30]))])
        T30[2] = np.mean([DailyTemp[x] + dz[seg2:]*lapse_rate for x in range(np.max([1, t-30]))])


        # Run the model one timespep forward... save result
        Snow, Wb, Wu, Qb, Qd, ET = one_time_forward(3,
                                                    Pd,
                                                    Td,
                                                    T30,
                                                    L,
                                                    Wu,
                                                    Wb,
                                                    Snow,
                                                    frtdir,
                                                    frtgw,
                                                    smcap,
                                                    etpar,
                                                    tmelt,
                                                    t_snow,
                                                    t_melt,
                                                    t_base,
                                                    t_power)

        Q[t] = Qb + Qd
        Snowholder[:,t] = Snow
    # compute the total discharge and return it

    return Q, Snowholder, tholder


if __name__ == '__main__':
    # Import stuff to read the forcings

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


    daily_temp = np.load("./data/daily_temp.npy")
    daily_precip = np.load("./data/daily_precip.npy")
    daily_q_observed = np.load("./data/daily_q_observed.npy")
    day_len_hrs = np.load("./data/day_len_hrs.npy")
    dz = np.load("./data/dz_reduced.npy")

    q, snow, tholder = ForwardModelFasterLayers(daily_temp, daily_precip, day_len_hrs, dz, M=1.0)
    plt.plot(q)
    plt.plot(daily_q_observed)
    # plt.plot(snow[0,:])
    # plt.plot(snow[1,:])
    # plt.plot(snow[2,:])

    # plt.plot(tholder[0,:], label='0')
    # plt.plot(tholder[1,:], label='1')
    # plt.plot(tholder[2,:], label='2')

    # start = time.time()
    # ForwardModelFaster(daily_temp, daily_precip, day_len_hrs, dz, M=1.0)
    # print('It took', time.time()-start, 'seconds.')



