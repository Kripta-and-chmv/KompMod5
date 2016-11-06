import scipy
import scipy.integrate
import normal_distrib as n_d
import random
import numpy as np

def bounds(u, v):
    prob_b = 1
    prob_a = 1    

    while prob_a > 0.01:
        a = int(random.uniform(-5, 0))
        prob_a, err = scipy.integrate.quad(n_d.pdf, -np.inf, a, args = (u, v))

    while prob_b > 0.01:
        b = int(random.uniform(a+1, 5))
        prob_b, err = scipy.integrate.quad(n_d.pdf, b, np.inf, args = (u, v))
    
    return a, b

def generate(u, v):
    a, b = bounds(u, v)
    M = []
    for x in range(a, b+1):
        M.append(n_d.pdf(x, u, v))
    M = max(M)

    x = a + random.random() * (b - a)
    y = random.random() * M
    
    while y >= n_d.pdf(x, u, v):
        x = a + random.random() * (b - a)
        y = random.random() * M    

    return x