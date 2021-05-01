# process the ASO data into elevation "bands" such
# that they can be comapred w/ the lumped snow model
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import pandas as pd
mpl.use('Qt5Agg')
import pathlib
import xarray as xr


def date_to_jday(date):
	pddate = pd.to_datetime(date)
	return pddate.dayofyear

def get_aso_snow(date, lower_lim, upper_lim):
	dem_base = pathlib.Path("/Volumes/Transcend/ASOdata/DEMs")
	dem = xr.open_rasterio(dem_base.joinpath("3mdem_upsample_50m_clipped_to_east.tif"))

	# PROCESS some ASO data
	base = pathlib.Path("/Volumes/Transcend/ASOdata/swe_data/50m_clipped_to_east")

	time_var = xr.Variable('time', ["2018-03-31", "2018-05-24", "2019-04-07", "2019-06-10"])
	ds = xr.concat([xr.open_rasterio(base.joinpath(x)) for x in base.glob("*USCOGE*.tif")], dim=time_var)


	# first layer
	asovals = ds.sel(time=date).values[0,:,:]
	asovals = np.where(asovals>=0, asovals, np.nan)
	demvals = dem.values[0,:,:]

	b1 = np.where((demvals >= lower_lim) & (demvals < upper_lim))
	return np.nanmean(asovals[b1])*1000.   # convert m --> mm


	# b1 = np.where((demvals > 2428.29) & (demvals < 2938.89))
	# b2 = np.where((demvals > 2938.89) & (demvals < 3264.06))
	# b3 = np.where((demvals > 3264.06) & (demvals < 3700.02))
	# b4 = np.where((demvals > 3700.02))



