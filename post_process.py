import seaborn as sb
import numpy as np
import pymc3 as pm
import scipy.stats as stats
import theano
import theano.tensor as tt
#from lumped_hydro_model import ForwardModel
import sys
import arviz as az
import matplotlib as mpl
import matplotlib.pyplot as plt
from attempt_3 import *

# import some cool fortran modules that are precompiled and fast
libPath = './fortran/'
sys.path.insert(0,libPath)
from driver import driver as dr
from ForcingModel import ForcingModel
import pickle
mpl.use('Qt5Agg')

# FUNCTIONS
print('foo')

def kge(mod, obs):
    #
    mobs = np.mean(obs)
    sobs = np.std(obs)

    # mean ratio
    b = np.mean(mod) / mobs
    # std
    a = np.std(mod) / sobs
    # corr coeff
    r = np.corrcoef(mod, obs)[0, 1]  # corrcoef returns the correlation matrix...
    # the diagonals are 1, the off-diags are the 'r'
    # value that we want
    kgeval = 1-np.sqrt((r - 1.)**2 + (a - 1.)**2 + (b - 1)**2)
    return kgeval


def bdk(pet_p_ratio, n=1.5):
    return pet_p_ratio/((pet_p_ratio)**n + 1.0)**(1/n)


# BEGIN MAIN
with open('./trace_2018-10-01_2019-09-30.pkl', 'rb') as buff:
    trace = pickle.load(buff)

# with open('./posterior_predictive.pkl', 'rb') as buff:
#     poster = pickle.load(buff)


# load data ...
start_date = "2018-10-01"
end_date = "2019-09-30"
data_base = "/Users/williamrudisill/Documents/Conceptual_Runoff_Model/data/"

# load data ...
frc = ForcingModel()
frc.read()

# read in forcings for the desired time period
daily_temp, daily_precip, qObs, day_len_hrs, daily_swe = frc(start_date, end_date)

# read these in (they are fixed in time)
elev = np.load(data_base + "elev.npy")
dz = np.load(data_base + "dz.npy")

# Load the snow17 parameters
pkl_file = open('./parameter_files/snow17params.pkl', 'rb')
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
i=0
for l,t in zip(day_len_hrs, daily_temp):
   PET[i] = (dr.pet(l,t))
   i+=1

##### Run Snow17 Model ######
print('foo')

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
opg_clb = snow17params.get('opg')#.002                       #         ! PARAMETER   SNOW17
bias_clb  = snow17params.get('bias') #.01                       #         ! PARAMETER   SNOW17
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

output_trace = np.empty((ntimes, 7000))
output_post = np.empty((ntimes, 7000))
bias = 0
kge_vals = np.empty(14000)

def runr(x):
    nr = 2.
    kr = 2.9
    sm1i = 150
    sm2i = 250
    output = dr.model_driver(snowop     = SNOWOP,       #         ! OPTION   SNOW option
                             etop       = ETOP,       #         ! OPTION   Percolation option
                             drainop    = DRAINOP,    #         ! OPTION   Snow option
                             satop      = SATOP,      #         ! OPTION   Saturated Area Option
                             baseop     = BASEOP,     #         ! OPTION   Surface Runoff option
                             smop1      = SMOP1,      #         ! OPTION   Soil Water Option -- Top
                             smop2      = SMOP2,      #         ! OPTION   Soil Water Option -- Bottom
                             ntimes     = ntimes,     #         ! FORCING   Number of model timesteps
                             pet        = PET,        #         ! FORCING   Potential Evapotranspiratio
                             jday       = jday,       #         ! FORCING   Day of Year
                             tair       = tair,       #         ! FORCING   Air Temperature, length N
                             precip     = precip,     #         ! FORCING   Precipitation, length N
                             nlayers    = nlayers,    #         ! PARAMETER   SNOW17
                             rvs        = rvs,        #         ! PARAMETER   SNOW17
                             opg_method = opg_method, #         ! PARAMETER   SNOW17
                             dz         = dz,         #         ! PARAMETER   SNOW17
                             dt         = dt,         #         ! PARAMETER   SNOW17
                             opg        = x[0],        #         ! PARAMETER   SNOW17
                             bias       = x[1],       #         ! PARAMETER   SNOW17
                             uadj       = uadj,       #         ! PARAMETER   SNOW17
                             mbase      = mbase,      #         ! PARAMETER   SNOW17
                             mfmax      = mfmax,      #         ! PARAMETER   SNOW17
                             mfmin      = mfmin,      #         ! PARAMETER   SNOW17
                             tipm       = tipm,       #         ! PARAMETER   SNOW17
                             nmf        = nmf,        #         ! PARAMETER   SNOW17
                             plwhc      = plwhc,      #         ! PARAMETER   SNOW17
                             pxtemp     = pxtemp,     #         ! PARAMETER   SNOW17
                             pxtemp1    = pxtemp1,    #         ! PARAMETER   SNOW17
                             pxtemp2    = pxtemp2,    #         ! PARAMETER   SNOW17
                             sm1i       = sm1i,
                             sm2i       = sm2i,
                             sm1max     = sm1max,      #      ! PARAMETER   ?
                             sm2max     = sm2max,      #      ! PARAMETER   ?
                             ku         = x[2],  #ku,         ! PARAMETER   PERCOLATION
                             c          = x[3],  #c,          ! PARAMETER   PERCOLATION
                             sm1fmax    = 0.0,   #            ! PARAMETER   PERCOLATION  --- OPTIONAL
                             psi        = 0.0,   #            ! PARAMETER   PERCOLATION  --- OPTIONAL
                             alpha      = 0.0,   #            ! PARAMETER   PERCOLATION  --- OPTIONAL
                             ks         = x[4],  #ks,         ! PARAMETER   BASEFLOW
                             lam        = 0.0,   #lam,        ! PARAMETER   BASEFLOW  --- OPTIONAL
                             lowercasen = x[5],  #lowercasen,   #         ! PARAMETER   BASEFLOW  --- OPTIONAL
                             beta       = x[6],  #beta,         #         ! PARAMETER   SFROFF
                             nr         = nr,
                             kr         = kr)

    qtot, qchan, qb, qsx, eVec, qin, sm1, sm2 =  output
    return qchan

for i in range(7000):
    x=[trace.get_values("opg")[i],
       trace.get_values("bias")[i],
       trace.get_values("ku")[i],
       trace.get_values("c")[i],
       trace.get_values("ks")[i],
       trace.get_values("lowercasen")[i],
       trace.get_values("beta")[i]]

    qchan = runr(x)
    output_post[:, i] = qchan
    kge_vals[i] = kge(qchan,qObs)

model_values = np.save("output_post", output_post)
np.save('qObs', qObs)
# kge_vals = np.where(kge_vals > 1.0, -.999)
# best_loc = np.argwhere(kge_vals == kge_vals.max())

# best_hydrograph = output_post[:, best_loc][:,0,0]


