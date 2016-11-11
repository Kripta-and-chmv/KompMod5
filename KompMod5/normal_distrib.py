import scipy as sc
import scipy.stats as scst
import numpy as np


def cdf(x, u, v):
    erf = sc.special.erf((x-u) / (v * np.sqrt(2)))
    return scst.norm.cdf(x, u, v)
    return 0.5 * (1 + erf)

def pdf(x, u, v):
    return scst.norm.pdf(x, u, v)
    return np.exp(-(x - u)**2 / (2 * v**2)) / (v * np.sqrt(2 * np.pi))