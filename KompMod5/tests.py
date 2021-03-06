import scipy as sc
import sys
import numpy as np
import normal_distrib as norm_d
import matplotlib.pyplot as plt

def smirnov(seq, u, v, alpha):
    def calc_d_plus(seq, u, v):
        """Вычисление D+"""
        d = []
        lng = len(seq)
        for i, x in zip(range(1, lng+1), seq):
            el = i/lng - norm_d.cdf(x, u, v)
            d.append(el)
        return max(d)

    def calc_d_minus(seq, u, v):
        """Вычисление D-"""
        d = []
        lng = len(seq)
        for i, x in zip(range(1, lng+1), seq):
            el = norm_d.cdf(x, u, v) - (i - 1)/lng
            d.append(el)
        return max(d)

    def calc_dn(seq, u, v):
        """Вычисление Dn"""
        d_min = calc_d_minus(seq, u, v)
        d_plus = calc_d_plus(seq, u, v)
        return max(d_min, d_plus)

    def calc_s_star(seq, u, v):
        """Вычисление S*"""
        dn = calc_dn(seq, u, v)
        lng = len(seq)
        s_star = (6 * lng * dn + 1)**2 / (9 * lng)
        return s_star

    def calc_prob_s_grtr_sstr(s_star):
        return np.exp(-s_star / 2)

    print("Тест Смирнова:")

    seq.sort()
    
    s_star = calc_s_star(seq, u, v)
    print("\tЗначение статистики - {}".format(s_star))

    prob_s = calc_prob_s_grtr_sstr(s_star)
    print("\tP(S* > S) - {}".format(prob_s))
    
    hit = prob_s > alpha
    print("\tРезультат прохождения теста - {}\n".format(hit))
    return hit

def chisqr_test(sequence, alpha, u, v):
    """Тест Хи-квадрат"""
    print("Тест хи квадрат:")

    max_in_seq = max(sequence)
    min_in_seq = min(sequence)
    len_seq = len(sequence)
    # разбиваем отрезок на интервалы
    intervals_amount = int(5 * sc.log10(len_seq))
    K = intervals_amount
    lngth = (max_in_seq - min_in_seq)/K
    intervals = [x * lngth + min_in_seq for x in range(0, K+1)]
    
    #определяем количество попаданий в интервалы
    hits_amount = []    
    for a, b in zip(intervals[:-1], intervals[1:]):
            count = sum([a <= x < b for x in sequence])
            hits_amount.append(count)

    emper_prob = [x / len_seq for x in hits_amount]

    # Вычисляется вероятность попадания слчайной величины в заданные
    # интервалы
    def calc_probs(intervals):
        return [norm_d.cdf(x, u, v) - norm_d.cdf(y, u, v) for x, y in zip(intervals[1:], intervals[:-1])]

    probabils = calc_probs(intervals)
    kk = sum(probabils)
    kkk = sum(emper_prob)
    def graph(intervals, probabils, emper_prob, width):
        #width = intervals[len(intervals) - 1] / (len(intervals) - 1)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.bar(intervals[:len(intervals) - 1], probabils, width, color="white", label = u'Theoretical')
        ax.bar(intervals[:len(intervals) - 1], emper_prob, width, alpha=0.5, color="black", label = u'Emperical')
        ax.legend(loc = 'best', frameon = True)
        plt.title('Chi2 Histogram')
        plt.xlabel('intervals')
        plt.ylabel('hits amount')
        plt.xticks(intervals)
        plt.show()

    graph(intervals, probabils, emper_prob, lngth)
    # вычисляется статистика
    addition = 0
    for hits, probs in zip(hits_amount, probabils):
        if probs == 0: continue
        addition += (hits / len_seq - probs)**2 / probs

    s_star = len(sequence) * addition
    print("\tЗначение статистики - {}".format(s_star))

    # вычисляется P(S*>S)
    r = intervals_amount - 1
    print("\tКоличество степеней свободы - {}".format(r))

    def integrand(x, r):
        return x ** (r / 2 - 1) * sc.exp(-x / 2)

    prob_s = sc.integrate.quad(integrand, s_star, np.inf, args = (r))
    prob_s = prob_s[0] / (2 ** (r / 2) * sc.special.gamma(r / 2))

    print("\tP(S*>S) - {}".format(prob_s))
    print("\tПрохождение теста хи квадрат - {}\n".format(prob_s > alpha))
