import scipy as sc
import numpy as np

def cdf(u, v, x):
    erf = sc.special.erf((x-u) / (v * numpy.sqrt(2)))
    return 0.5 * (1 + erf)

def pdf(u, v, x):
    return np.exp(-(x - u)**2 / (2 * v**2)) / (np.sqrt(2 * v**2 * np.pi))