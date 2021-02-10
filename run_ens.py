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
    return kgeval, r, b, a




with open('trace.pkl', 'rb') as buff:
    trace = pickle.load(buff)

with open('posterior_predictive.pkl', 'rb') as buff:
	poster = pickle.load(buff)


TdVec = np.load("./data/daily_temp.npy")
PdVec = np.load("./data/daily_precip.npy")
LVec = np.load("./data/day_len_hrs.npy")
dz = np.load("./data/dz_reduced.npy")
obs_data = np.load("./data/daily_q_observed.npy")


param_names = ['frtdir', 'frtgw', 'smcap', 'etpar', 't_snow', 't_melt', 't_base', 't_power', 'ogp']

output_trace = np.empty((365, 14000))
output_post = np.empty((365, 14000))

bias = 0

kge_vals = np.empty(14000)

# # do the trace
# for i in range(14000):
# 	frtdir = trace.get_values('frtdir')[i]
# 	frtgw = trace.get_values('frtgw')[i]
# 	smcap = trace.get_values('smcap')[i]
# 	etpar = trace.get_values('etpar')[i]
# 	t_snow = trace.get_values('t_snow')[i]
# 	t_melt = trace.get_values('t_melt')[i]
# 	t_base = trace.get_values('t_base')[i]
# 	t_power = trace.get_values('t_power')[i]
# #	bias = trace.get_values('bias')[i]
# 	opg = trace.get_values('ogp')[i]
# 	model = ffm.fwd(4, PdVec, TdVec, LVec, dz, frtdir, frtgw, smcap, etpar, t_snow, t_melt, t_base, t_power, bias, opg)
# 	output_trace[:, i] = model


for i in range(14000):
	frtdir =  poster.get_values('frtdir')[i]
	frtgw =   poster.get_values('frtgw')[i]
	smcap =   poster.get_values('smcap')[i]
	etpar =   poster.get_values('etpar')[i]
	t_snow =  poster.get_values('t_snow')[i]
	t_melt =  poster.get_values('t_melt')[i]
	t_base =  poster.get_values('t_base')[i]
	t_power = poster.get_values('t_power')[i]
#	bias = trace.get_values('bias')[i]
	opg = trace.get_values('ogp')[i]
	model = ffm.fwd(4, PdVec, TdVec, LVec, dz, frtdir, frtgw, smcap, etpar, t_snow, t_melt, t_base, t_power, bias, opg)
		[:, i] = model
	kgev, r,b,a = kge(model, obs_data)
	kge_vals[i] = kgev







