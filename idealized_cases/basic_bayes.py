import pymc3 as pm
import arviz as az
import matplotlib.pyplot as plt
import numpy as np
print(f"Running on PyMC3 v{pm.__version__}")

qObs = np.load('qObs.npy')



if __name__ == "__main__":
	# size = 200
	# true_intercept = 1
	# true_slope = 2

	# x = np.linspace(0, 1, size)
	# # y = a + b*x
	# true_regression_line = true_intercept + true_slope * x
	# # add noise
	# y = true_regression_line + np.random.normal(scale=0.5, size=size)

	# data = dict(x=x, y=y)


	# # fig = plt.figure(figsize=(7, 7))
	# # ax = fig.add_subplot(111, xlabel="x", ylabel="y", title="Generated data and underlying model")
	# # ax.plot(x, y, "x", label="sampled data")
	# # ax.plot(x, true_regression_line, label="true regression line", lw=2.0)
	# # plt.legend(loc=0);
	# # plt.show()


	# with pm.Model() as model:  # model specifications in PyMC3 are wrapped in a with-statement
	#     # Define priors
	#     sigma = pm.HalfCauchy("sigma", beta=10, testval=1.0)
	#     intercept = pm.Normal("Intercept", 0, sigma=20)
	#     x_coeff = pm.Normal("x", 0, sigma=20)

	#     # Define likelihood
	#     likelihood = pm.Normal("y", mu=intercept + x_coeff * x, sigma=sigma, observed=y)

	#     # Inference!
	#     trace = pm.sample(3000, cores=2)  # draw 3000 posterior samples using NUTS sampling

	with pm.Model() as model:
		b1 = pm.normal(mean=0, sigma=.25)
		b2 = pm.normal(mean=0, sigma=.25)
		mu =
		sig = pm.Normal(mean=0, sigma=)




