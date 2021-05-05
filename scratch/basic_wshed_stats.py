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
import seaborn as sb


mpl.use("Agg")

start = "1988-10-01"
end = "2019-09-30"
f = ForcingModel()
f.read()
daily_temp, daily_precip, q, day_len_hrs = f(start, end)


# compute mean
p_mean = daily_precip.groupby(daily_precip.index.year).sum()[1:].mean()
q_mean = q.groupby(q.index.year).sum()[1:].mean()

# compute std deviation
p_std = np.std(daily_precip.groupby(daily_precip.index.year).sum()[1:])
q_std = np.std(q.groupby(q.index.year).sum()[1:])


# compute cofficient of variation
print(p_mean, p_std, p_std/p_mean)
print(q_mean, q_std, q_std/q_mean)


sb.distplot(daily_precip.groupby(daily_precip.index.year).sum()[1:])
sb.distplot(q.groupby(q.index.year).sum()[1:])
plt.savefig('test')