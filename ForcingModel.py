import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import metpy.calc as mpcalc
from metpy.units import units
import sys
from lumped_hydro_model import *
import xarray as xr
# parameters

class ForcingModel:
    def __init__(self):
        self.discharge_unit_mm = None
        self.watershed_area = None
        self.hypsometry = None
        self.snotel_elev = None
        self.DailyTempForcing = None
        self.DailyPrecipForcing = None
        self.DayLength = None
        self.days = None

    def read(self):
        self.PrepareDEM()
        self.PrepareDischarge()
        self.PrepareDEM()
        self.PrepareSnowObs()
        self.ApplyDayLength(self.days)

    def PrepareDEM(self):
        ## ---  Get the basin as among other things ---##
        ds = xr.open_rasterio("/Volumes/Transcend/ASOdata/DEMs/3mdem_upsample_50m_clipped_to_east.tif").sel(band=1)
        flat_ds = ds.values.flatten()
        flat_ds_east = flat_ds[np.where(flat_ds > 0)]
        flat_ds_east_sort = np.sort(flat_ds_east)
        east_river_area_m2 = 748983000.0 # m2

        self.watershed_area = east_river_area_m2
        self.hypsometry = flat_ds_east_sort
        ds.close()

    def PrepareDischarge(self):
        ## --- Gather the USGS stream observations ----#
        columns = ['usgs', 'id', 'date', 'timezone', 'discharge', 'flag']
        usgs_df = pd.read_csv('usgs_data.txt', skiprows=41, sep='\t', names = columns)
        usgs_df['date'] = pd.to_datetime(usgs_df.date)
        usgs_df.set_index('date', inplace=True)

        # there might be missing data -- reinterpolate to correct time (adds timesteps)
        # go from hourly --> daily
        usgs_df = usgs_df.resample("D").mean()
        usgs_df = usgs_df.interpolate()

        del usgs_df['id']
        usgs_df['discharge_m3s'] = usgs_df['discharge']* 0.0283  # convert to m3/s
        usgs_df['discharge_unit_mm'] = usgs_df['discharge_m3s'] * (24*60*60)/self.watershed_area*1000
#       usgs_df = usgs_df.loc[start:end]
        self.discharge_unit_mm = usgs_df['discharge_unit_mm']

    def PrepareSnowObs(self):
        ## ---  Gather the Snotel Forcings  ----##
        database = "/Volumes/Transcend/EastRiverClimatePaper/Snotel/CO_snotel.db"
        dbengine = "sqlite:///{}".format(database)

        # select the snotel station X
        df = pd.read_sql_table('380', dbengine)
        df['Date'] = pd.to_datetime(df['Date'])
        df.set_index('Date', inplace=True)

        # get the number of times...
        ntimes = len(df.index)
        days = df.index.dayofyear.values

        # comptue/fix the forcing data
        DailyPrecipForcing = df.IncPrecip.fillna(0) * 25.4
        DailyTempForcing = ((df.TavgF-32.0)*5/9)
        DailyTempForcing = DailyTempForcing.fillna(DailyTempForcing.mean())

        self.snotel_elev = 3096.768 # meters
        self.DailyTempForcing = DailyTempForcing
        self.DailyPrecipForcing = DailyPrecipForcing
        self.days = days

    def ApplyDayLength(self, days):
        LenOfDayHr = [DayLength_Calc(d, lat=38.0) for d in days]
        self.DayLength = LenOfDayHr
        return LenOfDayHr

    def __call__(self, start, end):
        sub_dtf = self.DailyTempForcing.loc[start:end]
        sub_pcp = self.DailyPrecipForcing.loc[start:end]
        sub_q = self.discharge_unit_mm.loc[start:end]
        sub_doy = sub_pcp.index.dayofyear.values
        sub_day_len = self.ApplyDayLength(sub_doy)

        return sub_dtf, sub_pcp, sub_q, sub_day_len
        # sub_dlen = self.DailyTempForcing.loc[start:end]

if __name__ == "__main__":
    start = "2010-10-01"
    end = "2011-09-30"
    f = ForcingModel()
    f.read()
    daily_temp, daily_precip, day_len_hrs = f(start, end)



