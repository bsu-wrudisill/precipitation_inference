import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import metpy.calc as mpcalc
from metpy.units import units
import sys
from lumped_hydro_model import *
import xarray as xr
from ForcingModel import ForcingModel

from scipy.optimize import curve_fit
from scipy import stats
# parameters
import matplotlib as mpl

mpl.use("Agg")

# columns = ['usgs', 'id', 'date', 'timezone', 'discharge', 'flag']
# usgs_df = pd.read_csv('./data/USGS_almont.txt', skiprows=41, sep='\t', names = columns)
# usgs_df['date'] = pd.to_datetime(usgs_df.date)
# usgs_df.set_index('date', inplace=True)
# usgs_df['month'] = usgs_df.index.month
# usgs_df['doy'] = usgs_df.index.dayofyear


# #summer = usgs_df.loc[(usgs_df.month == 7) | (usgs_df.month == 8) | (usgs_df.month == 9)]

# def sin_fit(x, freq, amplitude, phase, offset):
#     return np.sin(x * freq + phase) * amplitude + offset


# freqlist = []
# amplitude = []
# phase = []
# offset = []

# avg_list = []

# one_year = usgs_df.loc[(usgs_df.index.year == 2018)]

# for d in range(100,250):
#     onedayusgs = one_year.loc[(one_year.doy == d)]
#     # get the time..
#     avg_list.append(np.mean(onedayusgs.discharge.values))
#     y = onedayusgs.groupby(onedayusgs.index.hour).mean().discharge.values
#     N = len(y)
#     t = np.linspace(0,2*np.pi,N)

#     try:
#         popt, pcov = curve_fit(sin_fit, t, y)
#         freqlist.append(popt[0])
#         amplitude.append(popt[1])
#         phase.append(popt[2])
#         offset.append(popt[3])

#     except RuntimeError:
#         freqlist.append(np.nan)
#         amplitude.append(np.nan)
#         phase.append(np.nan)
#         offset.append(np.nan)


#     #

start = "1988-10-01"
end = "2019-09-30"
f = ForcingModel()
f.read()
daily_temp, daily_precip, q, day_len_hrs = f(start, end)


fig, ax = plt.subplots(1,2)

##----- Do the whole year  -----
year = []
qsumlist =  []
psumlist = []

for i in range(1989,2019):
    qsum = q.loc["%s-10-01"%(i):"%s-09-30"%(i+1)].sum()
    psum = daily_precip.loc["%s-10-01"%(i):"%s-09-30"%(i+1)].sum()
    #psum = daily_precip.loc["%s-10-01"%(i):"%s-05-01"%(i+1)].sum()
    year.append(i+1)
    qsumlist.append(qsum)
    psumlist.append(psum)


qsumlist =  np.array(qsumlist)
psumlist = np.array(psumlist)

slope, intercept, r_value, p_value, std_err = stats.linregress(psumlist, qsumlist)
ax[0].plot(psumlist, intercept + slope*psumlist, color='red')
ax[0].scatter(psumlist, qsumlist)
ax[0].set_xlim(100,1000)
ax[0].set_ylim(100,1000)

for i, txt in enumerate(year):
    ax[0].annotate(txt, (psumlist[i], qsumlist[i]))

ax[0].set_ylabel("Annual Stream Discharge (mm) \n USGS 09112500 @ Almont ")
ax[0].set_xlabel("Snotel precipitation (mm) \n NRCS 380 Butte")
ax[0].set_title("Total precipitation vs Discharge")
ax[0].annotate("R^2: %s"%(round(r_value,3)), (200,600))

print(r_value)
##----- Now do Oct-May for snotel precip
year = []
qsumlist =  []
psumlist = []

for i in range(1989,2019):
    qsum = q.loc["%s-10-01"%(i):"%s-09-30"%(i+1)].sum()
    #psum = daily_precip.loc["%s-10-01"%(i):"%s-09-30"%(i+1)].sum()
    psum = daily_precip.loc["%s-10-01"%(i):"%s-04-01"%(i+1)].sum()
    year.append(i+1)
    qsumlist.append(qsum)
    psumlist.append(psum)


qsumlist =  np.array(qsumlist)
psumlist = np.array(psumlist)

slope, intercept, r_value, p_value, std_err = stats.linregress(psumlist, qsumlist)
ax[1].plot(psumlist, intercept + slope*psumlist, color='red')
ax[1].scatter(psumlist, qsumlist)

for i, txt in enumerate(year):
    ax[1].annotate(txt, (psumlist[i], qsumlist[i]))

ax[1].set_xlim(100,1000)
ax[1].set_ylim(100,1000)
ax[1].set_ylabel("Annual Stream Discharge (mm) \n USGS 09112500 @ Almont ")
ax[1].set_xlabel("Snotel precipitation (mm) \n NRCS 380 Butte")
ax[1].set_title("Winter precipitation vs Discharge")
ax[1].annotate("R^2: %s"%(round(r_value,3)), (200,600))
print(r_value)



plt.show()




