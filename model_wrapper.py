import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import pandas as pd
import pickle
import sys
import pickle
import argparse
import pymc3 as pm
import theano
import theano.tensor as tt

sys.path.append("/Users/williamrudisill/Documents/Conceptual_Runoff_Model/fortran/")
sys.path.append("/Users/williamrudisill/Documents/Conceptual_Runoff_Model/")
from driver import driver as dr
from ForcingModel import ForcingModel



# parser = argparse.ArgumentParser()
# parser.add_argument("start_date", type=str)
# parser.add_argument("end_date", type=str)
# args = parser.parse_args()

# # load data ...
# start_date = args.start_date#"2019-10-01"
# end_date = args.end_date#"2020-09-30"
# data_base = "/Users/williamrudisill/Documents/Conceptual_Runoff_Model/data/"


class driverclass:

    def __init__(self, start_date, end_date):
        # load data ...
        frc = ForcingModel()
        frc.read()

        # read in forcings for the desired time period
        daily_temp, daily_precip, obs_q, day_len_hrs, daily_swe = frc(start_date, end_date)
        self.obs_q = obs_q

        # read these in (they are fixed in time)
        data_base = "/Users/williamrudisill/Documents/Conceptual_Runoff_Model/data/"
        elev = np.load(data_base + "elev.npy")
        dz = np.load(data_base + "dz.npy")

        # Load the snow17 parameters
        pkl_file = open('./parameter_files/snow17params.pkl', 'rb')
        snow17params = pickle.load(pkl_file)
        pkl_file.close()

        pkl_file = open('./parameter_files/driver_params_calibrated.pkl', 'rb')
        driver_params = pickle.load(pkl_file)
        pkl_file.close()

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

        # process dates
        dates = pd.date_range(start_date, end_date, freq='D')
        jdays = dates.dayofyear.values * 1.0

        ntimes=len(dates)

        # create PET forcing
        PET = np.zeros_like(jdays)

        i=0
        for l,t in zip(day_len_hrs, daily_temp):
           PET[i] = (dr.pet(l,t))
           i+=1


        ##### Run Snow17 Model ######

        self.SNOWOP = 0                       #         ! OPTION   SNOW option
        self.ETOP = 0                         #         ! OPTION   Percolation option
        self.DRAINOP = 0                      #         ! OPTION   Snow option
        self.SATOP = 0                        #         ! OPTION   Saturated Area Option
        self.BASEOP = 0                       #         ! OPTION   Surface Runoff option
        self.SMOP1 = 0                        #         ! OPTION   Soil Water Option -- Top
        self.SMOP2 = 0                        #         ! OPTION   Soil Water Option -- Bottom
        self.ntimes = ntimes                  #         ! FORCING   Number of model timesteps
        self.PET = PET                        #         ! FORCING   Potential Evapotranspiratio
        self.jday = jdays                     #         ! FORCING   Day of Year
        self.tair = daily_temp                #         ! FORCING   Air Temperature = 1 length N
        self.precip = daily_precip            #         ! FORCING   Precipitation = 1 length N
        self.nlayers = 3                      #         ! PARAMETER   SNOW17

        # Snow model parameters
        self.rvs = 1                           #         ! PARAMETER   SNOW17
        self.opg_method = 1                    #         ! PARAMETER   SNOW17
        self.dz  = heights - 3096.768          #         ! PARAMETER   SNOW17
        self.dt  = 24                          #         ! PARAMETER   SNOW17
        self.opg_clb = snow17params.get('opg')#.002                       #         ! PARAMETER   SNOW17
        self.bias_clb  = snow17params.get('bias') #.01                       #         ! PARAMETER   SNOW17
        self.uadj  = snow17params.get('uadj')#.04                       #         ! PARAMETER   SNOW17
        self.mbase = snow17params.get('mbase')#.3                       #         ! PARAMETER   SNOW17
        self.mfmax = snow17params.get('mfmax')#4.0                      #         ! PARAMETER   SNOW17
        self.mfmin =snow17params.get('mfmin')# 1.0                      #         ! PARAMETER   SNOW17
        self.tipm  = snow17params.get('tipm')#.1                        #         ! PARAMETER   SNOW17
        self.nmf   = snow17params.get('nmf')#.4                         #         ! PARAMETER   SNOW17
        self.plwhc = snow17params.get('plwhc')#04                      #         ! PARAMETER   SNOW17
        self.pxtemp  = snow17params.get('pxtemp')#2.                      #         ! PARAMETER   SNOW17
        self.pxtemp1 = snow17params.get('pxtemp1')#-1.                    #         ! PARAMETER   SNOW17
        self.pxtemp2 = snow17params.get('pxtemp2')#3.                     #         ! PARAMETER   SNOW17

        # hydrologic model parameters ...
        self.sm1max     = driver_params.get("sm1max")
        self.sm2max     = driver_params.get("sm2max")
        self.ku         = driver_params.get("ku")
        self.c          = driver_params.get("c")
        self.ks         = driver_params.get("ks")
        self.lowercasen = driver_params.get("lowercasen")
        self.beta       = driver_params.get("beta")
        self.nr = 2.
        self.kr = 2.9
        self.sm1i = 150
        self.sm2i = 250



    def run(self, x):
        output = dr.model_driver(snowop      = self.SNOWOP,     #         ! OPTION   SNOW option
                               etop       = self.ETOP,       #         ! OPTION   Percolation option
                               drainop    = self.DRAINOP,    #         ! OPTION   Snow option
                               satop      = self.SATOP,      #         ! OPTION   Saturated Area Option
                               baseop     = self.BASEOP,     #         ! OPTION   Surface Runoff option
                               smop1      = self.SMOP1,      #         ! OPTION   Soil Water Option -- Top
                               smop2      = self.SMOP2,      #         ! OPTION   Soil Water Option -- Bottom
                               ntimes     = self.ntimes,     #         ! FORCING   Number of model timesteps
                               pet        = self.PET,        #         ! FORCING   Potential Evapotranspiratio
                               jday       = self.jday,       #         ! FORCING   Day of Year
                               tair       = self.tair,       #         ! FORCING   Air Temperature, length N
                               precip     = self.precip,     #         ! FORCING   Precipitation, length N
                               nlayers    = self.nlayers,    #         ! PARAMETER   SNOW17
                               rvs        = self.rvs,        #         ! PARAMETER   SNOW17
                               opg_method = self.opg_method, #         ! PARAMETER   SNOW17
                               dz         = self.dz,         #         ! PARAMETER   SNOW17
                               dt         = self.dt,         #         ! PARAMETER   SNOW17
                               opg        = x[0],                       #         ! PARAMETER   SNOW17
                               bias       = x[1],                       #         ! PARAMETER   SNOW17
                               uadj       = self.uadj,       #         ! PARAMETER   SNOW17
                               mbase      = self.mbase,      #         ! PARAMETER   SNOW17
                               mfmax      = self.mfmax,      #         ! PARAMETER   SNOW17
                               mfmin      = self.mfmin,      #         ! PARAMETER   SNOW17
                               tipm       = self.tipm,       #         ! PARAMETER   SNOW17
                               nmf        = self.nmf,        #         ! PARAMETER   SNOW17
                               plwhc      = self.plwhc,      #         ! PARAMETER   SNOW17
                               pxtemp     = self.pxtemp,     #         ! PARAMETER   SNOW17
                               pxtemp1    = self.pxtemp1,    #         ! PARAMETER   SNOW17
                               pxtemp2    = self.pxtemp2,    #         ! PARAMETER   SNOW17
                               sm1i       = self.sm1i,
                               sm2i       = self.sm2i,
                               sm1max     = self.sm1max,     #            ! PARAMETER   ?
                               sm2max     = self.sm2max,     #            ! PARAMETER   ?
                               ku         = x[2],                       #            ! PARAMETER   PERCOLATION
                               c          = x[3],                       #            ! PARAMETER   PERCOLATION
                               sm1fmax    = 0.0,                        #            ! PARAMETER   PERCOLATION  --- OPTIONAL
                               psi        = 0.0,                        #            ! PARAMETER   PERCOLATION  --- OPTIONAL
                               alpha      = 0.0,                        #            ! PARAMETER   PERCOLATION  --- OPTIONAL
                               ks         = x[4],                       #ks,         ! PARAMETER   BASEFLOW
                               lam        = 0.0,                        #lam,        ! PARAMETER   BASEFLOW  --- OPTIONAL
                               lowercasen = x[5],                       #lowercasen, ! PARAMETER   BASEFLOW  --- OPTIONAL
                               beta       = x[6],                       #beta,       !  PARAMETER   SFROFF
                               nr         = self.nr,
                               kr         = self.kr)

        qtot, qchan, qb, qsx, eVec, qin, sm1, sm2 =  output
        return qchan

