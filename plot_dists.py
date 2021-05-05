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
from basic_model import fast_foward_model as ffm
import pickle
from ForcingModel import ForcingModel

# begin
start_date = "2018-10-01"
end_date = "2019-09-30"
data_base = "/Users/williamrudisill/Documents/Conceptual_Runoff_Model/data/"

# load data ...
frc = ForcingModel()
frc.read()

# read in forcings for the desired time period
daily_temp, PdVec, qObs, day_len_hrs, daily_swe = frc(start_date, end_date)
dz = frc.hypsometry - 3096.0

obs_precip_total = PdVec.sum()
# with open('trace.pkl', 'rb') as buff:
#     trace = pickle.load(buff)
with open('prior.pkl', 'rb') as buff:
    prior_dict = pickle.load(buff)

with open('trace 11.pkl', 'rb') as buff:
    posterior_predictive = pickle.load(buff)

# Compute Precipitation -- make this a separate plot
fig,ax = plt.subplots(1)
sb.histplot(prior_dict['opg'], color='blue', kde=True, ax=ax, stat="probability", label='prior')
sb.histplot(posterior_predictive.get_values('opg'), color='red', kde=True, ax=ax, stat="probability", label='prior')
plt.show()

# prior_precip_list = []
# inferred_precip_list = []
# #precip_total = PdVec.sum()

# for ogp, bias in zip(posterior_predictive.get_values('opg'), posterior_predictive.get_values('bias')):
# 	PdVec_bias= np.where(PdVec > 0., PdVec + bias, 0)
# 	precip_total = PdVec.sum()
# 	inferred_precip_list.append(precip_total + np.mean(precip_total*ogp*dz))

# for ogp, bias in zip(prior_dict['opg'], prior_dict['bias']):
# 	PdVec_bias = np.where(PdVec > 0., PdVec + bias, 0)
# 	precip_total = PdVec_bias.sum()
# 	inferred_precip_list.append(precip_total + np.mean(precip_total*ogp*dz))



# fig,ax = plt.subplots(1)
# #sb.histplot(inferred_precip_list, color='red', kde=True, ax=ax, stat="probability", label='posterior')
# sb.histplot(prior_precip_list, color='blue', kde=True, ax=ax, stat="probability", label='prior')
# ax.set_xlabel("Precipitation mm")

# ax.axvline(obs_precip_total, color='black', ls='--')
# fig.legend()
# plt.show()


sys.exit()
# plt.savefig("inferred_precip_mm_2018", dpi=600)
# make all of the subplots

keylist = list(prior_dict.keys())
keylist.remove('loglike')
fig,ax = plt.subplots(2)

k=0
for key in ["opg", "bias"]:#, "t_power", "frtdir"]:
	posterior = posterior_predictive.get_values(key)
	prior = prior_dict.get(key)
	# make plots
	prior=sb.histplot(prior, kde=True, ax=ax[0,k], stat='density')
	poster=sb.histplot(posterior, kde=True, ax=ax[0,k], color='red', stat='density')
	ax[k].set_title(key)
	k+=1


# k=0
# for key in ["frtgw", "t_base", "smcap"]:
# 	posterior = posterior_predictive.get_values(key)
# 	prior = prior_dict.get(key)
# 	# make plots
# 	prior=sb.histplot(prior, kde=True, ax=ax[1,k], stat='density')
# 	poster=sb.histplot(posterior, kde=True, ax=ax[1,k], color='red', stat='density')
# 	ax[1,k].set_title(key)
# 	k+=1

ax[1,k].set_visible(False)

#labels=["prior"]
labels=["prior", "posterior"]

fig.legend([prior],     # The line objects
           labels=labels,   # The labels for each line
           loc="center right",   # Position of legend
           borderaxespad=0.1,    # Small spacing around legend box
           )

plt.tight_layout()
#plt.savefig('wy2018_params_prior_post', dpi=400)
plt.show()