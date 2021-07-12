import pymc3 as pm
import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import theano
import theano.tensor as tt
print(f"Running on PyMC3 v{pm.__version__}")



def ssl_like(func, parameters, data):
    model = func(parameters[:-1])
    sigma0 = parameters[-1]
    return -np.sum(0.5 * (np.log(2 * np.pi * sigma0 ** 2) + ((data - model) / sigma0) ** 2))

class likelihood_evaluate(tt.Op):
    itypes = [tt.dvector] # expects a vector of parameter values when called
    otypes = [tt.dscalar] # outputs a single scalar value (the log likelihood)

    def __init__(self, func, loglike, obs_data):
        # add inputs as class attributes
        self.func = func
        self.likelihood = loglike
        self.obs_data = obs_data

    def perform(self, node, inputs, outputs):
        # the method that is used when calling the Op
        parameters, = inputs  # this will contain my variables

        # call the log-likelihood function
        logl = self.likelihood(self.func, parameters, self.obs_data)
        outputs[0][0] = np.array(logl)



if __name__ == "__main__":
    size = 25

    x = np.linspace(1, 7, size)
    G = lambda m: m[0]*np.exp(m[1]*x) + m[2]*x*np.exp(m[3]*x)
    m = [1.0, -0.5, 1.0, -0.75]
    noise = np.random.randn(25)*(.1)
    data = G(m) + noise

    logl = likelihood_evaluate(G, ssl_like, data)

    with pm.Model() as model:
        m1 = pm.Uniform("m1", 0., 2.)
        m2 = pm.Uniform("m2", -1., 0.)
        m3 = pm.Uniform("m3", 0., 2.)
        m4 = pm.Uniform("m4", -1., 0.)
        sigma = pm.HalfNormal('sigma', sigma=.5)

        parameters = tt.as_tensor_variable([m1, m2, m3, m4, sigma])
        pm.Potential('loglike', logl(parameters))

        with model:
            step1 = pm.Metropolis(scaling=.005, tune=100000, discard_tuned_samples=True, chains=4)
            trace = pm.sample(step = step1)



        # pm.traceplot(trace)
        # plt.show()
        print(pm.summary(trace))


