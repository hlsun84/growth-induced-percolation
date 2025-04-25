import math
from scipy.special import comb, perm
from scipy.stats import binom
import numpy as np
from scipy.special import zeta

Accuracy = 1e-12

def gen_possion_distribution(ave_k, k_max):
    p_k = {}
    p_k[0] = math.exp(-ave_k)
    for kk in range(1, k_max + 1):
        p_k[kk] = p_k[kk - 1] * ave_k / kk
    return p_k

def gen_powerlaw_distribution(k_min, k_max, gamma):
    pk = {}
    sum = 0.0
    for k in np.arange(k_min, k_max + 1):
        #print(str(pow(k,-gamma)))
        temp = math.pow(k,-gamma)/zeta(gamma)
        pk[k] = temp
        sum += temp
    for k in pk:
        pk[k] = pk[k] / sum

    average_k = cal_averageK(pk)
    qk = pk2qk(pk, average_k)
    return pk, qk, average_k

def pk2qk(pk, average_k):
    qk = {}
    for k in pk:
        qk[k-1] = pk[k] * k / average_k
    try:
        del qk[-1]
    except:
        pass
    return qk

def degree_dis_del_nodes(pk, remain_p):
    pk2 = {}
    for k in range(max(pk.keys())+1):
        pk2[k] = 0.0
    for k in pk:
        temp_pk = binom.pmf(range(k+1),k,remain_p)
        for k2 in range(k+1):
            pk2[k2] += temp_pk[k2] * pk[k]
    average_k = cal_averageK(pk2)
    qk = pk2qk(pk2, average_k)
    return pk2, qk, average_k

def cal_averageK(pk):
    return sum([ k * pk[k] for k in pk])

def G0x(pk, average_k, x):
    return sum([pk[k] * pow(x, k) for k in pk])
