"""Summary
"""
import seaborn as sb
import numpy as np
import pymc3 as pm
import scipy.stats as stats
import theano
import theano.tensor as tt
#from lumped_hydro_model import ForwardModel
import sys
import arviz as az
import matplotlib.pyplot as plt
from attempt_2 import *

# import some cool fortran modules that are precompiled and fast
libPath = './fortran/'
sys.path.insert(0,libPath)
from myfastmodule import fast_foward_model as ffm
import pickle



with open('trace.pkl', 'rb') as buff:
    trace = pickle.load(buff)

with open('prior.pkl', 'rb') as buff:
    prior = pickle.load(buff)

with open('posterior_predictive.pkl', 'rb') as buff:
    posterior_predictive = pickle.load(buff)

with pm.Model() as model:
    fig,ax = plt.subplots(3)
    az.plot_posterior(posterior_predictive, var_names=['ogp', 'etpar', 't_melt'], ax=ax, color='red')
    az.plot_posterior(prior, var_names=['ogp', 'etpar', 't_melt'], ax=ax, color='black')
