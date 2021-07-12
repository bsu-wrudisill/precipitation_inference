import seaborn as sb
import numpy as np
import scipy.stats as stats
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy

mpl.use('Qt5Agg')
output_post = np.load('output_post.npy')
qObs = np.load('qObs.npy')
qObs = qObs.reshape(1, 365)

# res =np.subtract(output_post, qObs.T)

# for i in range(0,10):
# 	resi = res[:,i]
# 	ztest = (resi - resi.mean())/resi.std()
# 	theoretical = theoretical = np.random.randn(365)*resi.std() + resi.mean()
# 	plt.plot(sorted(ztest), sorted(theoretical))


	#wls error model...
#	sigma_t = alpha *
# plt.plot([-5,5], [-5,5])
# plt.xlim(-5,5)
# plt.ylim(-5,5)
# plt.xlabel("theoretical (normal) quantiles")
# plt.ylabel("residual quantiles")


def f(x, x0, alpha, beta):
	# wls model
    sigma_t = alpha * x0 + beta
    res = abs(x - x0)
    inside = (sigma_t - res)/x0
    return np.sqrt(np.sum(inside**2))



fmin = lambda params: f(output_post[:,1], qObs[0,:], params[0], params[1])

foo =scipy.optimize.minimize(fmin, [.1,1.])

# apply it ...
minq = qObs[0,:] + foo.x[0]*qObs[0,:] + foo.x[1]
maxq = qObs[0,:] - (foo.x[0]*qObs[0,:] + foo.x[1])
plt.plot(minq, label='minq')
plt.plot(maxq, label='maxq')
plt.plot(qObs[0,:], label='q')
plt.legend()
