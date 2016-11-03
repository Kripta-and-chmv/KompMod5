import scipy as sc
import sys
import numpy as np
import normal_distrib as norm_d

def smirnov(seq, v, u, alpha):
    def calc_d_plus(seq, v, u):
        d = []
        lng = len(seq)
        for i, x in zip(range(1, lng+1), seq):
            el = i/lng - norm_d.cdf(v, u, x)
            d.append(el)
        return max(d)

    def calc_d_minus(seq, v, u):
        d = []
        lng = len(seq)
        for i, x in zip(range(1, lng+1), seq):
            el = norm_d.cdf(v, u, x) - (i - 1)/lng
            d.append(el)
        return max(d)

    def calc_dn(seq, v, u):
        d_min = calc_d_minus(seq, v, u)
        d_plus = calc_d_plus(seq, v, u)
        return max(d_min, d_plus)

    def calc_s_star(seq, v, u):
        dn = calc_dn(seq, v, u)
        lng = len(seq)
        s_star = (6 * lng * dn + 1)**2 / (9 * np.sqrt(lng))
        return s_star

    def calc_prob_s_grtr_sstr(s_star):
        return np.exp(-s_star / 2)

    

    seq.sort()
    
    s_star = calc_s_star(seq, v, u)
    

    prob_s = calc_prob_s_grtr_sstr(s_star)
    

    
    hit = prob_s > alpha
    
    return hit

def chisqr_test(sequence, alpha, v, u):
    

    mod = max(sequence)
    len_seq = len(sequence)
    
    intervals_amount = int(5 * sc.log10(len_seq))
    K = intervals_amount
    lngth = mod/K   
    intervals = [x * lngth for x in range(0, K+1)]
    
    
    hits_amount = []    
    for a, b in zip(intervals[:-1], intervals[1:]):
            count = sum([a <= x < b for x in sequence])
            hits_amount.append(count)

    
    def calc_probs(intervals):
        return [norm_d.cdf(v, u, x) - norm_d.cdf(v, u, y) for x, y in zip(intervals[1:], intervals[:-1])]

    probabils = calc_probs(intervals)

    
    addition = 0
    for hits, probs in zip(hits_amount, probabils):
        if probs == 0: continue
        addition += (hits / len_seq - probs)**2 / probs

    s_star = len(sequence) * addition
    

    
    r = intervals_amount - 1
    

    def integrand(x, r):
        return x ** (r / 2 - 1) * sc.exp(-x / 2)

    prob_s = sc.integrate.quad(integrand, s_star, np.inf, args = (r))
    prob_s = prob_s[0] / (2 ** (r / 2) * sc.special.gamma(r / 2))

    hit = prob_s > alpha
    
    return hit