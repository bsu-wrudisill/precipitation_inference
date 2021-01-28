import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import metpy.calc as mpcalc
from metpy.units import units
import sys
from lumped_hydro_model import *
import xarray as xr

# parameters



parameters = parms()


# Gather the forcings
start = "2010-10-01"
end = "2011-09-30"
f = ForcingModel()
f.read()
daily_temp, daily_precip, day_len_hrs = f(start, end)


# Compute the temperature lapse rates..
lapse_rate = .0065 # degrees c/ m --> c/km
orog_gradient = 2/1000 # 2x precip per 1000m
dz = flat_ds_east_sort - snotel_elev


# Setup the Buckets
Wb = np.zeros_like(days, dtype='float')
Wu = np.zeros_like(days, dtype='float')
Snow = np.zeros_like(days, dtype='float')
Qb = np.zeros_like(days, dtype='float')
Qd = np.zeros_like(days, dtype='float')
ET = np.zeros_like(days, dtype='float')


# Run the model...
for t in range(1,ntimes):
    L = LenOfDayHr[t]

    # adjust precipitation
    Pd = np.mean(DailyPrecipForcing[t] + DailyPrecipForcing[t]*dz*orog_gradient)
    Td = np.mean(DailyTempForcing[t] + DailyTempForcing[t]*dz*lapse_rate)

    # get rid of nan values...
    if np.nan in [Pd, Td]:
        Q_array[t] = Q_array[t-1]
        Wb[t] = Wb[t-1]
        Wu[t] = Wu[t-1]
        Snow[t] = Snow[t-1]
        break

    # compute T30; ignore temps from before the starting point...
    T30 = np.nanmean(DailyTempForcing[np.max([0, t-30]):t])


    Snow_t, Wb_t, Wu_t, Qb_t, Qd_t, ET_t = one_time_forward(Pd, Td, T30, Wu[t-1], Wb[t-1], Snow[t-1], L, parameters)

    Qb[t] = Qb_t
    Qd[t] = Qd_t
    Wb[t] = Wb_t
    Wu[t] = Wu_t
    Snow[t] = Snow_t
    ET[t] = ET_t

# compute the total discharge ....
Q = Qb + Qd

fig,ax = plt.subplots(1)
plt.plot(usgs_df.discharge_unit_mm.values)
plt.plot(Qd+Qb)

plt.show()
