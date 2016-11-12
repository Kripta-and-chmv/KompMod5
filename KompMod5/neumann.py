import scipy
import scipy.stats
import normal_distrib as n_d
import random
import numpy as np
import math

def bounds(u, v):
    
    a = (scipy.stats.norm.ppf(0.01, u, v))
    b = (scipy.stats.norm.ppf(0.99, u, v))

    return a, b

def generate(u, v):
    a, b = bounds(u, v)
    arr_M = []
    step = (b - a) / 100
    i = a
    #for x in range(a, b+1):
    #    arr_M.append(n_d.pdf(x, u, v))
    #M = max(arr_M)
    while i < b :
        arr_M.append(n_d.pdf(i, u, v))
        i += step
    M = max(arr_M)        

    x = a + random.random() * (b - a)
    y = random.random() * M
    
    while y >= n_d.pdf(x, u, v):
        x = a + random.random() * (b - a)
        y = random.random() * M    

    return x