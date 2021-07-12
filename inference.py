"""Summary
"""
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

# import arviz as az
sys.path.append("/Users/williamrudisill/Documents/Conceptual_Runoff_Model/fortran/")
sys.path.append("/Users/williamrudisill/Documents/Conceptual_Runoff_Model/")
from driver import driver as dr
from ForcingModel import ForcingModel

#from snow17 import *
mpl.use('Qt5Agg')


def my_loglike(x, obs_data, sigma):
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



def main(dict_of_parameters, obs_q):
    # Get forcing data
    #obs_data = np.load("./data/daily_fakeq.npy")

    # Put the forcings into a list ...
    # CREATE HTE MODEL LIKLIHOOD FUNCTION
    # sigma for the log liklihood function ...
    sigma = 2.0
    tracename = "/Users/williamrudisill/Documents/Conceptual_Runoff_Model/mcmc_output_files/trace_%s_%s.pkl"%(start_date,end_date)
    priorname = "/Users/williamrudisill/Documents/Conceptual_Runoff_Model/mcmc_output_files/prior_%s_%s.pkl"%(start_date,end_date)

    # Create the logliklihood
    logl = LogLike(my_loglike, obs_q, sigma)

    with pm.Model() as model:
        # loop through the parameters and create the tt variables ...
        list_of_parameters = []
        for key in dict_of_parameters.keys():
            p0 = dict_of_parameters.get(key)
            p1 = pm.Normal(key, mu=p0[1], sigma=p0[2])
            list_of_parameters.append(p1)


        parameters = tt.as_tensor_variable(list_of_parameters)#[sm1max,sm2max,ku,c,ks,lowercasen,beta, opg, beta])
        pm.Potential('loglike', logl(parameters))
        prior = pm.sample_prior_predictive()

        with model:
            step1 = pm.Metropolis()
            trace = pm.sample(step = step1, chains=14, tune=2000, discard_tuned_samples=True)# start={'m':0.4, 'c':3})
            # posterior_predictive = pm.sample_posterior_predictive(trace, random_seed=)

            with open(tracename, 'wb') as buff:
                pickle.dump(trace, buff)

            # with open('posterior_predictive.pkl', 'wb') as buff:
            #     pickle.dump(trace, buff)

            with open(priorname, 'wb') as buff:
                pickle.dump(prior, buff)



# BEGIN #
#-------#
parser = argparse.ArgumentParser()
parser.add_argument("start_date", type=str)
parser.add_argument("end_date", type=str)
args = parser.parse_args()


# load data ...
start_date = args.start_date#"2019-10-01"
end_date = args.end_date#"2020-09-30"
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
#etpar = .09
i=0
for l,t in zip(day_len_hrs, daily_temp):
   PET[i] = (dr.pet(l,t))
   i+=1

##### Run Snow17 Model ######

SNOWOP = 0                       #         ! OPTION   SNOW option
ETOP = 0                         #         ! OPTION   Percolation option
DRAINOP = 0                      #         ! OPTION   Snow option
SATOP = 0                        #         ! OPTION   Saturated Area Option
BASEOP = 0                       #         ! OPTION   Surface Runoff option
SMOP1 = 0                        #         ! OPTION   Soil Water Option -- Top
SMOP2 = 0                        #         ! OPTION   Soil Water Option -- Bottom
ntimes = ntimes                  #         ! FORCING   Number of model timesteps
PET = PET                        #         ! FORCING   Potential Evapotranspiratio
jday = jdays                     #         ! FORCING   Day of Year
tair = daily_temp                #         ! FORCING   Air Temperature = 1 length N
precip = daily_precip            #         ! FORCING   Precipitation = 1 length N
nlayers = 3                      #         ! PARAMETER   SNOW17

# Snow model parameters
rvs = 1                           #         ! PARAMETER   SNOW17
opg_method = 1                    #         ! PARAMETER   SNOW17
dz  = heights - 3096.768          #         ! PARAMETER   SNOW17
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

sm1max     = driver_params.get("sm1max")
sm2max     = driver_params.get("sm2max")
ku         = driver_params.get("ku")
c          = driver_params.get("c")
ks         = driver_params.get("ks")
lowercasen = driver_params.get("lowercasen")
beta       = driver_params.get("beta")
nr = 2.
kr = 2.9
sm1i = 150
sm2i = 250


if __name__ == '__main__':
    print("pymc3 version: ", pm.__version__)

    dict_of_parameters = {"opg":  ["normal", opg_clb, .0009],
                          "bias": ["normal", .5, .05],
                          "ku":   ["normal", ku, 2.0],
                          "c":    ["normal", c, 2.0],
                          "ks":    ["normal", ks, 2.0],
                          "lowercasen": ["normal", lowercasen, .5],
                          "beta": ["normal", beta, .5]}

    # now run the "main"
    main(dict_of_parameters, obs_q)


