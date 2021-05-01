import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import metpy.calc as mpcalc
from metpy.units import units
import sys
from scipy import optimize

## ---  Gather the Snotel Forcings  ----##
database = "/Volumes/Transcend/EastRiverClimatePaper/Snotel/CO_snotel.db"
dbengine = "sqlite:///{}".format(database)

# select the snotel station X
df = pd.read_sql_table('380', dbengine)
df_schofield = pd.read_sql_table('737', dbengine)

# butte snotel
df['Date'] = pd.to_datetime(df['Date'])
df.set_index('Date', inplace=True)

# scho snotel
df_schofield['Date'] = pd.to_datetime(df_schofield['Date'])
df_schofield.set_index('Date', inplace=True)

# get the number of times...
ntimes = len(df.index)
days = df.index.dayofyear.values

# comptue/fix the forcing data
DailyPrecipForcing = df.IncPrecip.fillna(0) * 25.4
DailyTempForcing = ((df.TavgF-32.0)*5/9)
DailyTempForcing = DailyTempForcing.fillna(DailyTempForcing.mean())

# butte, schofield
snotel_elev = 3096.768, 3261.36

SWE = df.SWE*25.4 # convert to mm

# remove bad dates
df = df.loc[pd.to_datetime("1986-01-01"):]
df_schofield = df_schofield.loc[pd.to_datetime("1986-01-01"):]

# Filter out the big jumps in temperature
df.TavgF = df.TavgF.where(abs(df.TavgF.diff()) < 60, np.nan)
df_schofield.TavgF = df_schofield.TavgF.where(abs(df_schofield.TavgF.diff()) < 50, np.nan)

# SUSPECT DATES --- BUTTE
bd1 = (pd.to_datetime("2002-07-01"), pd.to_datetime("2003-10-01"))
bd2 = (pd.to_datetime("1989-08-20"), pd.to_datetime("1989-08-21"))


# Remove the data from the bad date periods
df.TavgF.loc[bd1[0]:bd1[1]] = np.nan
df.TavgF.loc[bd2[0]:bd2[1]] = np.nan


####################################################################################
# 1) Adjust the schof. snotel by the adiabatic lapse rate
# The other snotel is at a different elevation. adjust temperature accordingly
####################################################################################

# adjust the schofield temp by the adiabatic lapse rate
lapse_rate_f = -.00356 # defgrees F per 1000m
df_schofield['TavgF_adjusted_schof'] = df_schofield.TavgF + lapse_rate_f * (3261.36 - 3096.768)


####################################################################################
# 2)  Replace bad Butte values with values from the adjusted Schofield !!!
####################################################################################

replace_nan = lambda a, b: b if np.isnan(a) else a
fixed = df.TavgF.combine(df_schofield['TavgF_adjusted_schof'], replace_nan)
df['fixed_v1'] = fixed

# # Not sure why we have to remove this one.... it looks fine in the
# df['fixed_v1'].loc[pd.to_datetime("1989-08-19"):pd.to_datetime("1989-08-20")] = np.nan

# Cut out the bad/missing data from the combined dataset (could do this eralier )


####################################################################################
# 3) Fit a sinusoidal function to fill in data for regions where niether station is good. should be small ish
####################################################################################

# interpolate bad values (must be done for curve fitting)
fixed_int = fixed.interpolate()

# take the weekly average
fixed_roll = fixed_int.rolling(7).mean() # 7 day rolling mean average. smooths out the sinish wave


# Remove the nans
valid = ~(np.isnan(fixed_roll))
ydata = fixed_roll.values[valid]

# create x data
xdata = np.linspace(0, len(ydata), len(ydata))

# fucntion to be fitted.
def sin_func(x, a, w, c, d):
    return a*np.sin(w*x + c) + d

# initial guess of the parameters
p0 = [30., 2*np.pi/365., 0., 30.]

# fit the curve
p, pcov = optimize.curve_fit(sin_func, xdata, ydata, p0=p0)


## make a plot to compare
## plt.plot(xdata, ydata)
## plt.plot(xdata, sin_func(xdata, *p))

df['fitted_data'] = np.nan # ugh more nans.
df.fitted_data.iloc[7:] = sin_func(xdata, *p) # we lose the first week when we do the rolling mean (and mask nans)
######################################
# 4) Sub in the fitted values
######################################
replace_nan = lambda a, b: b if np.isnan(a) else a
fixed_v2 = df.fixed_v1.combine(df.fitted_data, replace_nan)




#######################################
# 5) Create output timeseries
######################################
df = df.iloc[7:]
df['final_T2'] = fixed_v2
#

# df['TavgF'].plot(label='butte')
#df_schofield['TavgF'].plot(label='scho')
#df['fitted_data'].plot(label='sinusoidal')
df['final_T2'].plot(label='final')
df['fixed_v1'].plot(label='v1')

plt.legend()
plt.show()





