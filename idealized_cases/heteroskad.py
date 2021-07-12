import pymc3 as pm
import arviz as az
import matplotlib.pyplot as plt
import numpy as np
print(f"Running on PyMC3 v{pm.__version__}")

qObs = np.load('qObs.npy')



if __name__ == "__main__":
	size = 200
	true_intercept = 1
	true_slope = 2

	x = np.linspace(0, 10, size)

	# y = a + b*x
	true_regression_line = true_intercept + true_slope * x

	# add noise
	sigma_true = 2
	y = true_regression_line + true_regression_line * np.random.normal(scale=sigma_true,size=size)

	plt.plot(x,y,'x')
	plt.plot(x,true_regression_line)
	plt.show()

	with pm.Model() as model:
		alpha = pm.HalfNormal("alpha", 10.)
		beta = pm.Normal("beta", 0., 2.)
		ep = pm.Normal("ep", 0., 2.)

		sigma = alpha * true_regression_line + beta


		likelihood= pm.Normal("ypred", mu = true_intercept + true_slope*sigma*x, sd=ep, observed=y)
		trace=pm.sample()
		print(pm.summary(trace))





