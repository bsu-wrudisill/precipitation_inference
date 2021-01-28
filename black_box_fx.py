import numpy as np
import pymc3 as pm
import theano
import theano.tensor as tt
theano.config.exception_verbosity='high'


def line(theta, x):
    return theta[1] + x * theta[0]


def my_loglike(theta, x, data, sigma):
    """
    A Gaussian log-likelihood function for a model with parameters given in theta
    """

    model = line(theta, x)
    return -(0.5/sigma**2)*np.sum((data - model)**2)


class LogLike(tt.Op):
    itypes = [tt.dvector] # expects a vector of parameter values when called
    otypes = [tt.dscalar] # outputs a single scalar value (the log likelihood)

    def __init__(self, loglike, data, x, sigma):

        # add inputs as class attributes
        self.likelihood = loglike
        self.x = x
        self.data=data
        self.sigma=sigma

    def perform(self, node, inputs, outputs):
        # the method that is used when calling the Op
        theta, = inputs  # this will contain my variables

        # call the log-likelihood function
        logl = self.likelihood(theta, self.x, self.data, self.sigma)
        outputs[0][0] = np.array(logl)


if __name__ == '__main__':

    # set up our data
    N = 10  # number of data points
    sigma = 1.  # standard deviation of noise
    x = np.linspace(0., 9., N)

    mtrue = 0.4  # true gradient
    ctrue = 3.   # true y-intercept

    truemodel = line([mtrue, ctrue], x)

    # make data
    np.random.seed(716742)  # set random seed, so the data is reproducible each time
    data = sigma*np.random.randn(N) + truemodel

    ndraws = 3000  # number of draws from the distribution
    nburn = 1000   # number of "burn-in points" (which we'll discard)

    logl = LogLike(my_loglike, data, x, sigma)

    with pm.Model():
        # your external function takes two parameters, a and b, with Uniform priors
        m = pm.Uniform('m', lower=-10., upper=10.)
        c = pm.Uniform('c', lower=-10., upper=10.)

        # convert m and c to a tensor vector
        theta = tt.as_tensor_variable([m, c])

        # use a DensityDist (use a lamdba function to "call" the Op)
        pm.DensityDist('likelihood', lambda v: logl(v), observed={'v': theta})
        trace = pm.sample(ndraws, tune=nburn, discard_tuned_samples=True, start={'m':0.4, 'c':3})

   # _ = pm.traceplot(trace, lines={"m": mtrue, "c": ctrue})





