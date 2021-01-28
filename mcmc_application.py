import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pymc3 as pm
import theano.tensor as tt

from pymc3 import Model, Normal, Slice, sample, traceplot
from pymc3.distributions import Interpolated
from scipy import stats
from theano import as_op
from lumped_hydro_model import *


# Load the model forcings and other stuff....
daily_temp = np.load("daily_temp.npy")
daily_precip = np.load("daily_precip.npy")
daily_q_observed = np.load("daily_q_observed.npy")
day_len_hrs = np.load("day_len_hrs.npy")
dz = np.load("dz.npy")


if __name__ == '__main__':

	plt.style.use("seaborn-darkgrid")
	print(f"Running on PyMC3 v{pm.__version__}")


	# Initialize random number generator
	np.random.seed(93457)


	basic_model = Model()

	with basic_model:
		orog_gradient = Normal("orog_gradient", mu=.002, sigma = .5)
#    	Multiplier    = Normal("M", mu=1, sigma = .5)
		frtdir        = Normal("frtdir", mu=.09, sigma= .05)
		frtgw         = Normal("frtgw",  mu=.06, sigma=.05)


		# Expected value of outcome
		mu = ForwardModel(daily_temp,
			              daily_precip,
			              day_len_hrs,
			              dz,
						  orog_gradient = orog_gradient,
         			      frtdir = frtdir,
		                  frtgw = frtgw)


		# Likelihood (sampling distribution) of observations
		Y_obs = Normal("Y_obs", mu=mu, sigma=1, observed=daily_q_observed)

		# draw 1000 posterior samples
#		trace = sample(1000)

#	traceplot(trace)
#	plt.show()
