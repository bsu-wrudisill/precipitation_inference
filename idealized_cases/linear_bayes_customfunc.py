import pymc3 as pm
import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import theano
import theano.tensor as tt
import sys
print(f"Running on PyMC3 v{pm.__version__}")

#qObs = np.load('qObs.npy')


# def ssl_like(func, parameters, data):
#     qpred = func(parameters[:-1])
#     sigma = parameters[-1]
#     return (-.5 / sigma**2)*(np.sum(qpred - data)**2)

def my_loglike(theta, sigma0, x, data, sigma):
    model = my_model(theta, x)
    data0 = data
    return -np.sum(0.5 * np.log(2 * np.pi * sigma0 ** 2) + ((data0 - model) / sigma0) ** 2)

def my_loglike_wls(theta, alpha, beta, x, data, sigma):
    model = my_model(theta, x)
    data0 = data
    weights = data0 * alpha + beta
    return -np.sum(0.5 * np.log(2 * np.pi * sigma0 ** 2) + ((weights*data0 - model) / sigma0) ** 2)


class likelihood_evaluate(tt.Op):
    itypes = [tt.dvector] # expects a vector of parameter values when called
    otypes = [tt.dscalar] # outputs a single scalar value (the log likelihood)

    def __init__(self, func, obs_data):
        # add inputs as class attributes
        self.func = func
        self.data = data

    def perform(self, node, inputs, outputs):
        # the method that is used when calling the Op
        parameters, = inputs  # this will contain my variables
        sigma0 = parameters[-1]
        # call the log-likelihood function
#        log_likelihood = (.5 /sigma **2) * (np.sum(self.func(parameters) - self.data)**2)
        model = self.func(parameters[:-1])
        log_likelihood = -np.sum(0.5 * (np.log(2 * np.pi * sigma0 ** 2) + ((self.data - model) / sigma0) ** 2))
        outputs[0][0] = np.array(log_likelihood)



if __name__ == "__main__":
    size = 100

    x = np.linspace(0, 1, size)
    G = lambda m: m[0]*x + m[1]
    m = [2.0, 1.0]

    # make some random noise... this simulates the error in the real world
    noise = np.random.randn(size)*.2
    true = G(m)
    data = G(m) + noise*x
    logl = likelihood_evaluate(G, data)

    plt.plot(x, true)
    plt.plot(x, data, 'x')
    plt.show()

    # with pm.Model() as model:  # model specifications in PyMC3 are wrapped in a with-statement
    #     # Define priors
    #     intercept = pm.Uniform("Intercept", 0., 3.)
    #     x_coeff = pm.Uniform("x", 0., 3.)
    #     sigma = pm.HalfCauchy("sigma", beta=10, testval=1.0)

        # Define likelihood
        # likelihood = pm.Normal("y", mu=intercept + x_coeff * x, sigma=sigma, observed=data)

        # # Inference!
        # #step1 = pm.Metropolis()
        # trace = pm.sample(cores=2)  # draw 3000 posterior samples using NUTS sampling
        # print(pm.summary(trace))


    #     pm.traceplot(trace)
    #   plt.show()

    with pm.Model() as model:
        m1 = pm.Uniform("m1", 0., 3)
        m2 = pm.Uniform("m2", 0., 3.)
        sigma = pm.HalfNormal('sigma', sigma=1)
        beta = pm.HalfNormal('beta', sigma=1)

        parameters = tt.as_tensor_variable([m1, m2,sigma, beta])

        pm.Potential('loglike', logl(parameters))

        with model:
            step1 = pm.Metropolis(tune=100000, discard_tuned_samples=True, chains=4)
            trace = pm.sample(step = step1)

            #pm.traceplot(trace)
            #plt.show()
            print(pm.summary(trace))

        #plot up the lines...




        for i in range(1000):
            plt.plot(x, trace["m1"][i]*x + trace["m2"][i], color='blue', alpha = .02)
            # plt.plot(x, trace["x"][i]*x + trace["Intercept"][i], color='blue', alpha = .02)
        plt.show()
        plt.plot(x,data,'x', alpha=.001)
        plt.plot(x,G(m))
        plt.plot()
        # # now plot confidence intervals...
        # plt.plot(x, G(m) + .2 * 1.96, color='red')
        # plt.plot(x, G(m) - .2 * 1.96, color='red')

        # plt.show()

