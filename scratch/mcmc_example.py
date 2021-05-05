import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pymc3 as pm
import theano.tensor as tt

from pymc3 import Model, Normal, Slice, sample, traceplot
from pymc3.distributions import Interpolated
from scipy import stats
from theano import as_op


if __name__ == '__main__':

	plt.style.use("seaborn-darkgrid")
	print(f"Running on PyMC3 v{pm.__version__}")


	# Initialize random number generator
	np.random.seed(93457)

	# True parameter values
	alpha_true = 5
	beta0_true = 7
	beta1_true = 13

	# Size of dataset
	size = 100

	# Predictor variable
	X1 = np.random.randn(size)
	X2 = np.random.randn(size) * 0.2

	# Simulate outcome variable
	Y = alpha_true + beta0_true * X1 + beta1_true * X2 + np.random.randn(size)
	#Y0 = alpha_true + beta0_true * X1 + beta1_true * X2
	basic_model = Model()

	with basic_model:

		# Priors for unknown model parameters
		alpha = Normal("alpha", mu=0, sigma=1)
		beta0 = Normal("beta0", mu=12, sigma=1)
		beta1 = Normal("beta1", mu=18, sigma=1)

		# Expected value of outcome
		mu = alpha + beta0 * X1 + beta1 * X2

		# Likelihood (sampling distribution) of observations
		Y_obs = Normal("Y_obs", mu=mu, sigma=1, observed=Y)

		# draw 1000 posterior samples
		trace = sample(1000)

	traceplot(trace)
	plt.show()
