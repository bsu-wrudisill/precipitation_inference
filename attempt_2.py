"""Summary
"""
import numpy as np
import pymc3 as pm
import theano
import theano.tensor as tt
#from lumped_hydro_model import ForwardModel
import sys
import arviz as az
import matplotlib.pyplot as plt
import pickle


# import some cool fortran modules that are precompiled and fast
libPath = './fortran/'
sys.path.insert(0,libPath)
from myfastmodule import fast_foward_model as ffm

theano.config.exception_verbosity='high'


# def line(theta, x):
#     return theta[1] + x * theta[0]



def my_loglike(parameters, forcings, obs_data, sigma):
    """
    A Gaussian log-likelihood function for a model with parameters given in theta
    """
    frtdir = parameters[0]     # Parameter
    frtgw  = parameters[1]     # Parameter
    smcap  = parameters[2]     # Parameter
    etpar  = parameters[3]     # Parameter
    t_snow = parameters[4]     # Parameter
    t_melt = parameters[5]    # Parameter
    t_base = parameters[6]    # Parameter
    t_power = parameters[7]    # Parame
    # bias = parameters[8]
    opg  = parameters[8]


    # frtdir = .09
    # frtgw  = .06
    # smcap = 100,      # Parameter
    # etpar = .005,     # Parameter
    # tmelt = .01,      # Parameter
    # t_snow = 0.0,     # Parameter
    # t_melt = 1.0,      # Parameter
    # t_base = 0,      # Parameter
    # t_power = .5   # Parame
    # ogp  = .02
    bias = 0
    PdVec, TdVec, LVec, dz = forcings



    model = ffm.fwd(PdVec, TdVec, LVec, dz, frtdir, frtgw, smcap, etpar, t_snow, t_melt, t_base, t_power, bias, opg)

    # Compute the liklihood function
    return -(0.5/sigma**2)*np.sum((obs_data - model)**2)


class LogLike(tt.Op):
    itypes = [tt.dvector] # expects a vector of parameter values when called
    otypes = [tt.dscalar] # outputs a single scalar value (the log likelihood)

    def __init__(self, loglike, forcings, obs_data, sigma):

        # add inputs as class attributes
        self.likelihood = loglike
        self.forcings = forcings
        self.obs_data = obs_data
        self.sigma = sigma

    def perform(self, node, inputs, outputs):
        # the method that is used when calling the Op
        parameters, = inputs  # this will contain my variables

        # call the log-likelihood function
        logl = self.likelihood(parameters, self.forcings, self.obs_data, self.sigma)
        outputs[0][0] = np.array(logl)


def main():
    # Get forcing data
    TdVec = np.load("daily_temp.npy")
    PdVec = np.load("daily_precip.npy")
    LVec = np.load("day_len_hrs.npy")
    dz = np.load("dz_reduced.npy")
    obs_data = np.load("daily_q_observed.npy")

    # Put the forcings into a list ...
    forcings = [PdVec, TdVec, LVec, dz]


    # CREATE HTE MODEL LIKLIHOOD FUNCTION
    # sigma for the log liklihood function ...
    sigma = 10.0

    # Create the logliklihood
    logl = LogLike(my_loglike, forcings, obs_data, sigma)

    with pm.Model() as model:

        # Create the paramters that the model uses
        frtdir = pm.Normal("frtdir", mu=.09, sigma= .03)
        frtgw  = pm.Normal("frtgw",  mu=.06, sigma=.03)
        smcap  = pm.Normal("smcap",  mu=250, sigma=75)
        etpar  = pm.Normal("etpar",  mu=.05, sigma=.2)
        t_snow = pm.Normal("t_snow",  mu=0., sigma=1.0)
        t_melt = pm.Normal("t_melt",  mu=0., sigma=1.0)
        t_base = pm.Normal("t_base",  mu=.06, sigma=.05)
        t_power = pm.Normal("t_power",  mu=.5, sigma=.03)
        # bias    = pm.Normal("bias",  mu=0, sigma=.05)
        ogp     = pm.distributions.continuous.Uniform("ogp", lower=0, upper=.0025)


        # convert m and c to a tensor vector
        parameters = tt.as_tensor_variable([frtdir, frtgw, smcap, etpar, t_snow, t_melt, t_base, t_power, ogp])

        pm.Potential('loglike', logl(parameters))

        prior = pm.sample_prior_predictive()

        with model:
            step1 = pm.Metropolis()
            trace = pm.sample(step = step1, chains=14, tune=2000)#tune=nburn, discard_tuned_samples=True)# start={'m':0.4, 'c':3})
            posterior_predictive = pm.sample_posterior_predictive(trace)

            with open('trace.pkl', 'wb') as buff:
                pickle.dump(trace, buff)

            with open('posterior_predictive.pkl', 'wb') as buff:
                pickle.dump(trace, buff)

            with open('prior.pkl', 'wb') as buff:
                pickle.dump(prior, buff)




if __name__ == '__main__':
    print("pymc3 version: ", pm.__version__)
    main()


