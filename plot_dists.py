import numpy as np
import sys
import seaborn as sb
import scipy.stats as stats
import matplotlib.pyplot as plt
import pickle
import argparse
import pymc3 as pm
# import some cool fortran modules that are precompiled and fast
libPath = '/Users/williamrudisill/Documents/Conceptual_Runoff_Model/'
sys.path.insert(0,libPath)
libPath = '/Users/williamrudisill/Documents/Conceptual_Runoff_Model/fortran/'
sys.path.insert(0,libPath)
from attempt_3 import *
from ForcingModel import ForcingModel


# parser = argparse.ArgumentParser()
# parser.add_argument("trace", type=str)
# parser.add_argument("prior", type=str)
# args = parser.parse_args()

## trace = "trace_2013-10-01_2014-09-30.pkl"
## prior = "prior_2013-10-01_2014-09-30.pkl"

# start_date = args.prior.split("_")[1]
# end_date =  args.prior.split("_")[2].split(".")[0]

yearlist = range(20018, 2019)
num_plots = len(yearlist)
fig,axx = plt.subplots(1)

# loop through years
# for year,axx in zip(yearlist, ax):
for year in [2018]:

    # starting and ending
    start, end = "%s-10-01"%(year), "%s-09-30"%(year+1)
    prior = "prior_%s_%s.pkl"%(start,end)
    trace = "trace_%s_%s.pkl"%(start,end)

    with open(prior, 'rb') as buff:
        prior_dict = pickle.load(buff)

    with open(trace, 'rb') as buff:
        posterior_predictive = pickle.load(buff)
        pm.traceplot(posterior_predictive)

    # data_base = "/Users/williamrudisill/Documents/Conceptual_Runoff_Model/data/"

    # # load data ...
    # frc = ForcingModel()
    # frc.read()
    # daily_temp, daily_precip, obs_q, day_len_hrs, daily_swe = frc(start, end)
    # dz = frc.hypsometry - 3096.0
    # obs_precip_total = daily_precip.sum()

    # prior_precip_list = []
    # inferred_precip_list = []


    # for ogp, bias in zip(posterior_predictive.get_values('opg'), posterior_predictive.get_values('bias')):
    #     daily_precip_bias= np.where(daily_precip > 0., daily_precip + bias, 0)
    #     precip_total = daily_precip.sum()
    #     inferred_precip_list.append(precip_total + np.mean(precip_total*ogp*dz))

    # for ogp, bias in zip(prior_dict['opg'], prior_dict['bias']):
    #     daily_precip_bias = np.where(daily_precip > 0., daily_precip + bias, 0)
    #     precip_total = daily_precip_bias.sum()
    #     prior_precip_list.append(precip_total + np.mean(precip_total*ogp*dz))



    # sb.histplot(inferred_precip_list, color='red', kde=True, ax=axx, stat="density", label='posterior', alpha=.2)
    # sb.histplot(prior_precip_list, color='blue', kde=True, ax=axx, stat="density", label='prior', alpha=.2)
    # axx.set_xlabel("WY %s Precipitation (mm)"%(year))
    # axx.axvline(obs_precip_total, color='black', ls='--')
#    axx.set_ylim(0,.2)
#    axx.set_xlim(450.,800.)

    #fig.legend()

plt.show()


#sys.exit()
# plt.savefig("inferred_precip_mm_2018", dpi=600)
# make all of the subplots

# keylist = list(prior_dict.keys())
# keylist.remove('loglike')
# fig,ax = plt.subplots(2)

# k=0
# for key in ["opg", "bias"]:#, "t_power", "frtdir"]:
# 	posterior = posterior_predictive.get_values(key)
# 	prior = prior_dict.get(key)
# 	# make plots
# 	prior=sb.histplot(prior, kde=True, ax=ax[0,k], stat='density')
# 	poster=sb.histplot(posterior, kde=True, ax=ax[0,k], color='red', stat='density')
# 	ax[k].set_title(key)
# 	k+=1


# # k=0
# # for key in ["frtgw", "t_base", "smcap"]:
# # 	posterior = posterior_predictive.get_values(key)
# # 	prior = prior_dict.get(key)
# # 	# make plots
# # 	prior=sb.histplot(prior, kde=True, ax=ax[1,k], stat='density')
# # 	poster=sb.histplot(posterior, kde=True, ax=ax[1,k], color='red', stat='density')
# # 	ax[1,k].set_title(key)
# # 	k+=1

# ax[1,k].set_visible(False)

# #labels=["prior"]
# labels=["prior", "posterior"]

# fig.legend([prior],     # The line objects
#            labels=labels,   # The labels for each line
#            loc="center right",   # Position of legend
#            borderaxespad=0.1,    # Small spacing around legend box
#            )

# plt.tight_layout()
# #plt.savefig('wy2018_params_prior_post', dpi=400)
# plt.show()