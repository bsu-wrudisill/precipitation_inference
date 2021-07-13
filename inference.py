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
from model_wrapper import *



#from snow17 import *
mpl.use('Qt5Agg')

# def wls_like(func, parameters, data):
#     qpred = func(parameters[:-1])
#     sigma = parameters[-1]
#     w = 0.30014989
#     b = 0.72752388
#     return (-.5 / sigma**2)*(np.sum((w*data+b)*(qpred - data)**2))

def wls_like(func, parameters, data):
    model = func(parameters[:-1])
    sigma = parameters[-1]
#    return (-.5 / sigma**2)*(np.sum(qpred - data)**2)
    return -np.sum(0.5 * (np.log(2 * np.pi * sigma ** 2) + ((data - model) / sigma) ** 2))




class likelihood_evaluate(tt.Op):
    itypes = [tt.dvector] # expects a vector of parameter values when called
    otypes = [tt.dscalar] # outputs a single scalar value (the log likelihood)

    def __init__(self, func, loglike, obs_data):
        # add inputs as class attributes
        self.func = func
        self.likelihood = loglike
        self.obs_data = obs_data

    def perform(self, node, inputs, outputs):
        # the method that is used when calling the Op
        parameters, = inputs  # this will contain my variables

        # call the log-likelihood function
        logl = self.likelihood(self.func, parameters, self.obs_data)
        outputs[0][0] = np.array(logl)



def main(func, dict_of_parameters, obs_q):
    # Get forcing data
    #obs_data = np.load("./data/daily_fakeq.npy")

    # Put the forcings into a list ...
    # CREATE HTE MODEL LIKLIHOOD FUNCTION
    # sigma for the log liklihood function ...
#    sigma = 2.0
    tracename = "/Users/williamrudisill/Documents/Conceptual_Runoff_Model/trace.pkl"
    priorname = "/Users/williamrudisill/Documents/Conceptual_Runoff_Model/prior.pkl"

    # Create the logliklihood
    logl = likelihood_evaluate(func, wls_like, obs_q)

    with pm.Model() as model:
        # loop through the parameters and create the tt variables ...
        list_of_parameters = []
        print(list_of_parameters)
        for key in dict_of_parameters.keys():
            p0 = dict_of_parameters.get(key)
            p1 = pm.Normal(key, mu=p0[1], sigma=p0[2])
            list_of_parameters.append(p1)

        list_of_parameters.append(pm.HalfNormal("sigma", sigma=2.0))
        parameters = tt.as_tensor_variable(list_of_parameters)#[sm1max,sm2max,ku,c,ks,lowercasen,beta, opg, beta])
        pm.Potential('loglike', logl(parameters))
        prior = pm.sample_prior_predictive()

        with model:
            step1 = pm.Metropolis()
            trace = pm.sample(step = step1, chains=14, tune=2000, discard_tuned_samples=True)# start={'m':0.4, 'c':3})
            posterior_predictive = pm.sample_posterior_predictive(trace)

            with open(tracename, 'wb') as buff:
                pickle.dump(trace, buff)

            with open('posterior_predictive.pkl', 'wb') as buff:
                pickle.dump(trace, buff)

            with open(priorname, 'wb') as buff:
                pickle.dump(prior, buff)

            pm.traceplot(trace)
            plt.show()
            print(pm.summary(post_pred))


def run_post(model, params):
    prior = "prior.pkl"
    trace = "trace.pkl"
    post_pred = "posterior_predictive.pkl"
    # with open(prior, 'rb') as buff:
    #     prior_dict = pickle.load(buff)
    with open(trace, 'rb') as buff:
        post_pred = pickle.load(buff)

    num_iterations = post_pred[post_pred.varnames[0]].shape[0]
    matrix = np.empty((num_iterations,model.ntimes))
    for i in range(num_iterations):
        matrix[i,:] = model.run([post_pred[p][i] for p in params.keys()])

    return matrix


if __name__ == '__main__':
    print("pymc3 version: ", pm.__version__)


    model = driverclass("2010-09-01", "2011-09-01")
    real_values = [0.00250, .1, 87.582, 0.956, 4.404, 1.567, 1.5298]
    #qobs = model.run(real_values)

    dict_of_parameters = {"opg":  ["normal", model.opg_clb, .0009],
                          "bias": ["normal", .5, .05],
                          "ku":   ["normal", model.ku, 2.0],
                          "c":    ["normal", model.c, 2.0],
                          "ks":    ["normal", model.ks, 2.0],
                          "lowercasen": ["normal", model.lowercasen, .5],
                          "beta": ["normal", model.beta, .5]}#
                          # "sigma": ["normal", 1.0, 10.]}

    main(model.run, dict_of_parameters, model.obs_q)
    # matrix = run_post(model, dict_of_parameters)
    # prior = "prior.pkl"
    # trace = "trace.pkl"
    # post_pred = "posterior_predictive.pkl"
    # # with open(prior, 'rb') as buff:
    # #     prior_dict = pickle.load(buff)
    #with open(trace, 'rb') as buff:
    #     post_pred = pickle.load(buff)


    #pm.traceplot(post_pred)
    # plt.show()