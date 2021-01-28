"""Summary
"""
import numpy as np
import pymc3 as pm
import theano
import theano.tensor as tt
from lumped_hydro_model import ForwardModel
import sys

theano.config.exception_verbosity='high'


# def line(theta, x):
#     return theta[1] + x * theta[0]



def my_loglike(parameters, forcings, obs_data, sigma):
    """
    A Gaussian log-likelihood function for a model with parameters given in theta
    """
    orog_gradient = parameters[0]
    Multiplier = parameters[1]
    frtdir = parameters[2]
    frtgw = parameters[3]

    daily_temp, daily_precip, day_len_hrs, dz = forcings

    model = ForwardModel(
                         daily_temp,
                         daily_precip,
                         day_len_hrs,
                         dz,
                         M = Multiplier,
                         orog_gradient = orog_gradient,
                         frtdir = frtdir,
                         frtgw = frtgw)


    # Compute the liklihood function
    return -(0.5/sigma**2)*np.sum((obs_data - model)**2)


class LogLike(tt.Op):
    itypes = [tt.dvector] # expects a vector of parameter values when called
    otypes = [tt.dscalar] # outputs a single scalar value (the log likelihood)

    def __init__(self, loglike, obs_data, x, sigma):

        # add inputs as class attributes
        self.likelihood = loglike
        self.forcings = forcings
        self.obs_data=obs_data
        self.sigma=sigma

    def perform(self, node, inputs, outputs):
        # the method that is used when calling the Op
        parameters, = inputs  # this will contain my variables

        # call the log-likelihood function
        logl = self.likelihood(parameters, self.forcings, self.obs_data, self.sigma)
        outputs[0][0] = np.array(logl)


if __name__ == '__main__':

    # Get forcing data
    daily_temp = np.load("daily_temp.npy")
    daily_precip = np.load("daily_precip.npy")
    day_len_hrs = np.load("day_len_hrs.npy")
    dz = np.load("dz.npy")
    # daily_q_observed = np.load("daily_q_observed.npy")

    # Put the forcings into a list ...
    forcings = [daily_temp, daily_precip, day_len_hrs, dz]

    Multiplier_true    = 1.5
    orog_gradient_true = 1./100
    frtdir_true        = .0001
    frtgw_true         = .06


    truemodel = ForwardModel(*forcings,
                             M = Multiplier_true,
                             orog_gradient = orog_gradient_true,
                             frtdir = frtdir_true,
                             frtgw = frtgw_true)

    # create the "observed" data
    obs_data = 2.0*np.random.randn(len(truemodel)) + truemodel

    # sigma for the log liklihood function ...
    sigma = 2.0

    # Create the logliklihood
    logl = LogLike(my_loglike, obs_data, forcings, sigma)

    ndraws = 100  # number of draws from the distribution
    nburn = 10   # number of "burn-in points" (which we'll discard)

    with pm.Model():

        # Create the paramters that the model uses
        orog_gradient = pm.Normal("orog_gradient", mu=.002, sigma = .5)
        Multiplier    = pm.Normal("M", mu=1, sigma = .25)
        frtdir        = pm.Normal("frtdir", mu=.09, sigma= .05)
        frtgw         = pm.Normal("frtgw",  mu=.06, sigma=.05)


        # convert m and c to a tensor vector
        parameters = tt.as_tensor_variable([orog_gradient, Multiplier, frtdir, frtgw])

        # use a DensityDist (use a lamdba function to "call" the Op)
        pm.DensityDist('likelihood', lambda v: logl(v), observed={'v': parameters})

        step1 = pm.Metropolis()
#        trace = pm.sample(ndraws, tune=nburn, discard_tuned_samples=True)# start={'m':0.4, 'c':3})
        trace = pm.sample(ndraws, step=step1)
   # _ = pm.traceplot(trace, lines={"m": mtrue, "c": ctrue})





