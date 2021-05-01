import sys
sys.path.append("/Users/williamrudisill/Documents/Conceptual_Runoff_Model/fortran/")
sys.path.append("/Users/williamrudisill/Documents/Conceptual_Runoff_Model/")

from driver import driver as dr
from ForcingModel import ForcingModel

import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import pandas as pd
from snow17 import *
mpl.use('Qt5Agg')
import scipy.optimize
from scipy.special import gamma
import pickle

# load data ...
start_date = "2010-10-01"
end_date = "2020-09-30"
data_base = "/Users/williamrudisill/Documents/Conceptual_Runoff_Model/data/"

# load data ...
frc = ForcingModel()
frc.read()

# read in forcings for the desired time period
daily_temp, daily_precip, q, day_len_hrs, daily_swe = frc(start_date, end_date)

# read these in (they are fixed in time)
elev = np.load(data_base + "elev.npy")
dz = np.load(data_base + "dz.npy")

# Load the snow17 parameters
pkl_file = open('snow17params.pkl', 'rb')
snow17params = pickle.load(pkl_file)
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
etpar = .09
i=0
for l,t in zip(day_len_hrs, daily_temp):
   PET[i] = (dr.pet(l,t,etpar))
   i+=1


##### Run Snow17 Model ######

SNOWOP = 0                       #         ! OPTION   SNOW option
ETOP = 0                         #         ! OPTION   Percolation option
DRAINOP = 0                      #         ! OPTION   Snow option
SATOP = 0                        #         ! OPTION   Saturated Area Option
BASEOP = 0                       #         ! OPTION   Surface Runoff option
SMOP1 = 0                        #         ! OPTION   Soil Water Option -- Top
SMOP2 = 0                        #         ! OPTION   Soil Water Option -- Bottom
ntimes = ntimes                 #         ! FORCING   Number of model timesteps
PET = PET                        #         ! FORCING   Potential Evapotranspiratio
jday = jdays                     #         ! FORCING   Day of Year
tair = daily_temp                #         ! FORCING   Air Temperature = 1 length N
precip = daily_precip            #         ! FORCING   Precipitation = 1 length N
nlayers = 3                      #         ! PARAMETER   SNOW17

# Snow model parameters
rvs = 1                          #         ! PARAMETER   SNOW17
opg_method = 1                   #         ! PARAMETER   SNOW17
dz  = heights - 3096.768                    #         ! PARAMETER   SNOW17
dt  = 24                          #         ! PARAMETER   SNOW17
opg = snow17params.get('opg')#.002                       #         ! PARAMETER   SNOW17
bias  = snow17params.get('bias') #.01                       #         ! PARAMETER   SNOW17
uadj  = snow17params.get('uadj')#.04                       #         ! PARAMETER   SNOW17
mbase = snow17params.get('mbase')#.3                       #         ! PARAMETER   SNOW17
mfmax = snow17params.get('mfmax')#4.0                      #         ! PARAMETER   SNOW17
mfmin =snow17params.get('mfmin')# 1.0                      #         ! PARAMETER   SNOW17
tipm  = snow17params.get('tipm')#.1                        #         ! PARAMETER   SNOW17
nmf   = snow17params.get('nmf')#.4                         #         ! PARAMETER   SNOW17
plwhc = snow17params.get('plwhc')#04                      #         ! PARAMETER   SNOW17
pxtemp  = snow17params.get('pxtemp')#2.                      #         ! PARAMETER   SNOW17
pxtemp1 = snow17params.get('pxtemp1')#-1.                    #         ! PARAMETER   SNOW17
pxtemp2 = snow17params.get('pxtemp2')#3.                     #         ! PARAMETER   SNOW17

# hydrologic model parameters ...
sm1max = 200.                    #         ! PARAMETER   ?
sm2max = 3000.                   #         ! PARAMETER   ?
ku = .01                         #         ! PARAMETER   PERCOLATION
c = 1.                           #         ! PARAMETER   PERCOLATION
sm1Fmax = 240                    #         ! PARAMETER   PERCOLATION  --- OPTIONAL
psi = .1                         #         ! PARAMETER   PERCOLATION  --- OPTIONAL
alpha = .1                       #         ! PARAMETER   PERCOLATION  --- OPTIONAL
ks = .4                          #         ! PARAMETER   BASEFLOW
lam = .1                         #         ! PARAMETER   BASEFLOW  --- OPTIONAL
lowercasen = 1.                  #         ! PARAMETER   BASEFLOW  --- OPTIONAL
beta = 1.5                       #         ! PARAMETER   SFROFF
nr = 2.
kr = 2.9

foo = dr.model_driver(snowop   = SNOWOP,     #         ! OPTION   SNOW option
                      etop     = ETOP,       #         ! OPTION   Percolation option
                      drainop  = DRAINOP,    #         ! OPTION   Snow option
                      satop    = SATOP,      #         ! OPTION   Saturated Area Option
                      baseop   = BASEOP,     #         ! OPTION   Surface Runoff option
                      smop1    = SMOP1,      #         ! OPTION   Soil Water Option -- Top
                      smop2    = SMOP2,      #         ! OPTION   Soil Water Option -- Bottom
                      ntimes   = ntimes,     #         ! FORCING   Number of model timesteps
                      pet        = PET,        #         ! FORCING   Potential Evapotranspiratio
                      jday       = jday,       #         ! FORCING   Day of Year
                      tair       = tair,       #         ! FORCING   Air Temperature, length N
                      precip     = precip,     #         ! FORCING   Precipitation, length N
                      nlayers    = nlayers,    #         ! PARAMETER   SNOW17
                      rvs        = rvs,        #         ! PARAMETER   SNOW17
                      opg_method = opg_method, #         ! PARAMETER   SNOW17
                      dz         = dz,         #         ! PARAMETER   SNOW17
                      dt         = dt,         #         ! PARAMETER   SNOW17
                      opg        = opg,        #         ! PARAMETER   SNOW17
                      bias       = bias,       #         ! PARAMETER   SNOW17
                      uadj       = uadj,       #         ! PARAMETER   SNOW17
                      mbase      = mbase,      #         ! PARAMETER   SNOW17
                      mfmax      = mfmax,      #         ! PARAMETER   SNOW17
                      mfmin      = mfmin,      #         ! PARAMETER   SNOW17
                      tipm       = tipm,       #         ! PARAMETER   SNOW17
                      nmf        = nmf,        #         ! PARAMETER   SNOW17
                      plwhc      = plwhc,      #         ! PARAMETER   SNOW17
                      pxtemp     = pxtemp,     #         ! PARAMETER   SNOW17
                      pxtemp1    = pxtemp1,    #         ! PARAMETER   SNOW17
                      pxtemp2    = pxtemp2,     #           ! PARAMETER   SNOW17
                      sm1max     = sm1max,        #         ! PARAMETER   ?
                      sm2max     = sm2max,       #         ! PARAMETER   ?
                      ku         = ku,           #         ! PARAMETER   PERCOLATION
                      c          = c,            #         ! PARAMETER   PERCOLATION
                      sm1fmax    = sm1Fmax,      #         ! PARAMETER   PERCOLATION  --- OPTIONAL
                      psi        = psi,          #         ! PARAMETER   PERCOLATION  --- OPTIONAL
                      alpha      = alpha,        #         ! PARAMETER   PERCOLATION  --- OPTIONAL
                      ks         = ks,           #         ! PARAMETER   BASEFLOW
                      lam        = lam,      #        ! PARAMETER   BASEFLOW  --- OPTIONAL
                      lowercasen = lowercasen,   #         ! PARAMETER   BASEFLOW  --- OPTIONAL
                      beta       = beta,         #         ! PARAMETER   SFROFF
                      nr         = nr,
                      kr         = kr)

# # routing function...
# def ht(t, k=3.5, N=4):
#     a = (t/k)**(N-1) * np.exp(-t/k)/k*gamma(N)
#     return a/a.sum()

# t=np.linspace(0,len(jdays)-1,len(jdays))
# uht = ht(t)
# nx = len(foo)
# nb = len(foo)
# ny = len(foo)*2 - 1
# c=dr.convolve(nx=nx, nb=nb, ny=ny, bb=uht, xx=foo)

#plt.plot(foo[0])
plt.plot(foo[1])
plt.plot(q.values)
plt.show()

#np.convolve(uht,foo)
# plt.plot(foo)
# plt.plot(q.values)
# plt.show()



