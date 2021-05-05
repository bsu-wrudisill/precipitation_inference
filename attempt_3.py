"""Summary
"""
import numpy as np
import pymc3 as pm
import theano
import theano.tensor as tt
#from lumped_hydro_model import ForwardModel
import sys
import arviz as az
import pickle
import sys
sys.path.append("/Users/williamrudisill/Documents/Conceptual_Runoff_Model/fortran/")
sys.path.append("/Users/williamrudisill/Documents/Conceptual_Runoff_Model/")
from driver import driver as dr
from ForcingModel import ForcingModel

import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import pandas as pd
#from snow17 import *
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
daily_temp, daily_precip, obs_q, day_len_hrs, daily_swe = frc(start_date, end_date)

# read these in (they are fixed in time)
elev = np.load(data_base + "elev.npy")
dz = np.load(data_base + "dz.npy")

# Load the snow17 parameters
pkl_file = open('fortran/utils/snow17params.pkl', 'rb')
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


# x=[200.,
#    3000.,
#    .01,
#    1.,
#    240,
#    .1,
#    .1,
#    .4,
#    .1,
#    1.,
#    1.5,
#    2.,
#    2.9]

def my_loglike(x, obs_data, sigma):
    qr,qchan,snow = dr.model_driver(snowop   = SNOWOP,     #         ! OPTION   SNOW option
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
                                    opg        = x[10],        #         ! PARAMETER   SNOW17
                                    bias       = x[11],       #         ! PARAMETER   SNOW17
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
                                    sm1max     = x[0],#sm1max,        #         ! PARAMETER   ?
                                    sm2max     = x[1],#sm2max,       #         ! PARAMETER   ?
                                    ku         = x[2],#ku,           #         ! PARAMETER   PERCOLATION
                                    c          = x[3],#c,            #         ! PARAMETER   PERCOLATION
                                    sm1fmax    = sm1Fmax,      #         ! PARAMETER   PERCOLATION  --- OPTIONAL
                                    psi        = x[4],#psi,          #         ! PARAMETER   PERCOLATION  --- OPTIONAL
                                    alpha      = x[5],#alpha,        #         ! PARAMETER   PERCOLATION  --- OPTIONAL
                                    ks         = x[6],#ks,           #         ! PARAMETER   BASEFLOW
                                    lam        = x[7],#lam,      #        ! PARAMETER   BASEFLOW  --- OPTIONAL
                                    lowercasen = x[8],#lowercasen,   #         ! PARAMETER   BASEFLOW  --- OPTIONAL
                                    beta       = x[9],#beta,         #         ! PARAMETER   SFROFF
                                    nr         = nr,
                                    kr         = kr)

    return -(0.5/sigma**2)*np.sum((obs_data - qchan)**2)



class LogLike(tt.Op):
    itypes = [tt.dvector] # expects a vector of parameter values when called
    otypes = [tt.dscalar] # outputs a single scalar value (the log likelihood)

    def __init__(self, loglike, obs_data, sigma):

        # add inputs as class attributes
        self.likelihood = loglike
        self.obs_data = obs_data
        self.sigma = sigma

    def perform(self, node, inputs, outputs):
        # the method that is used when calling the Op
        parameters, = inputs  # this will contain my variables

        # call the log-likelihood function
        logl = self.likelihood(parameters, self.obs_data, self.sigma)
        outputs[0][0] = np.array(logl)



def main():
    # Get forcing data
    #obs_data = np.load("./data/daily_fakeq.npy")

    # Put the forcings into a list ...
    # CREATE HTE MODEL LIKLIHOOD FUNCTION
    # sigma for the log liklihood function ...
    sigma = 2.0

    # Create the logliklihood
    logl = LogLike(my_loglike, obs_q, sigma)

    with pm.Model() as model:
        sm1max     = pm.Normal("sm1max", mu=150., sigma=20.)
        sm2max     = pm.Normal("sm2max", mu=600., sigma=30.)
        ku         = pm.Normal("ku", mu=.01, sigma=.001)
        c          = pm.Normal("c", mu=1., sigma=.01)
        psi        = pm.Normal("psi", mu=.1, sigma=.001)
        alpha      = pm.Normal("alpha", mu=.1, sigma=.001)
        ks         = pm.Normal("ks", mu=.4, sigma=.001)
        lam        = pm.Normal("lam", mu=.1, sigma=.001)
        lowercasen = pm.Normal("lowercasen", mu=.1, sigma=.01)
        beta       = pm.Normal("beta", mu=1.5, sigma=.01)
        opg        = pm.Normal("opg", mu=opg_clb, sigma=.001)
        bias       = pm.Normal("bias", mu=bias_clb, sigma=.1)

        parameters = tt.as_tensor_variable([sm1max,sm2max,ku,c,psi,alpha,ks,lam,lowercasen,beta, opg, bias])

        pm.Potential('loglike', logl(parameters))

        prior = pm.sample_prior_predictive()

        with model:
            step1 = pm.Metropolis()
            trace = pm.sample(step = step1, chains=14, tune=2000)#tune=nburn, discard_tuned_samples=True)# start={'m':0.4, 'c':3})
            # posterior_predictive = pm.sample_posterior_predictive(trace, random_seed=)

            with open('trace.pkl', 'wb') as buff:
                pickle.dump(trace, buff)

            # with open('posterior_predictive.pkl', 'wb') as buff:
            #     pickle.dump(trace, buff)

            with open('prior.pkl', 'wb') as buff:
                pickle.dump(prior, buff)




if __name__ == '__main__':
    print("pymc3 version: ", pm.__version__)
    main()


