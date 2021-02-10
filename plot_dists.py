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


PdVec = np.load("./data/daily_precip.npy")
dz = np.load("./data/dz_reduced.npy")



with open('trace.pkl', 'rb') as buff:
    trace = pickle.load(buff)

with open('prior.pkl', 'rb') as buff:
    prior = pickle.load(buff)

with open('posterior_predictive.pkl', 'rb') as buff:
    posterior_predictive = pickle.load(buff)


ppvec = posterior_predictive.get_values('ogp').reshape(14000, 1)
dzvec = dz.reshape(1, 599)
avgvec = np.ones((1, 599))/599


precip = PdVec.sum() + PdVec.sum()*ppvec.dot(dzvec).dot(avgvec.T)

# fig,ax = plt.subplots()
# sb.histplot(posterior_predictive.get_values('ogp'), kde=True, ax=ax)
# #sb.histplot(prior.get('ogp'), kde=True, ax=ax, color='red')
# plt.show()