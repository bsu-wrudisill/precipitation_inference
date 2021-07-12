import numpy as np
import pymc3 as pm
import theano
import theano.tensor as tt
theano.config.exception_verbosity='high'
import matplotlib.pyplot as plt

if __name__ == '__main__':

    N = 50
    mtrue = 2.1
    btrue = 1.1

    # real data
    x = np.linspace(0,1,N)
    y_obs = x*mtrue + btrue + np.random.randn(N)*.2

    # Noise model
    mu_true = .9
    b_true = .01

    with pm.Model() as model:
        x_coeff = pm.Normal("x", 1.8, .8)
        intercept = pm.Normal("Intercept", .5, 1.)
        sigma = pm.HalfNormal("sigma", 2.)

        # this is the ... simulated data...?

        # make the likelihood function
        likelihood = pm.Normal("y", mu=x*x_coeff+intercept, sigma=sigma, observed=y_obs)

        trace=pm.sample()
    # pm.traceplot(trace)
    # plt.show()

    plt.figure(figsize=(7, 7))
    plt.plot(x, y_obs, "x", label="data")
    pm.plot_posterior_predictive_glm(trace, samples=100, label="posterior predictive regression lines")
    plt.plot(x, x*mu_true+b_true,  label="true regression line", lw=3.0, c="y")

    plt.title("Posterior predictive regression lines")
    plt.legend(loc=0)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()
