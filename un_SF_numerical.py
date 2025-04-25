import matplotlib.pyplot as plt
import numpy as np
import math
import function_general
import numpy
from scipy.special import comb, perm
import networkx as nx
import os

Accuracy = 1e-15
k_min=1
k_max=1000


def get_xy(q,p_k,q_k,para_m):
    x_new = q
    y_new = 0
    delta = 1
    while delta > Accuracy:
        x = x_new
        y = y_new
        x_new = 1 - (1 - q) * sum([q_k[k] * (1 - y_new) ** (k) for k in q_k])
        temp_sum = 0
        for k in range(k_max):
            temp = sum([comb(k,ss) * pow(1-x,k-ss) * (pow(x,ss) - (1-q)*pow(x-y,ss)) for ss in range(min(para_m-1,k)+1)])
            temp_sum += temp * q_k[k]
        y_new = x_new - temp_sum
        delta = max(abs(x-x_new), abs(y-y_new))
    return x_new, y_new


def get_abr(q,x,y,p_k,q_k,para_m):
    a_new = 1
    b_new = x
    r_new = y
    delta = 1
    while delta > Accuracy:
        a = a_new
        b = b_new
        r = r_new

        temp_sum = 0
        for k in range(k_max):
            temp = sum([comb(k,ss) * pow(x-b,ss) * (pow(1-x,k-ss) - pow(1-x-a+b,k-ss)) for ss in range(min(para_m-2,k)+1)])
            temp_sum += temp * q_k[k]
        temp_sum_a = sum([q_k[k] * pow(1-a,k) for k in q_k])
        a_new = 1 - temp_sum_a - temp_sum

        temp_sum = 0
        for k in range(k_max):
            temp = sum([comb(k,ss) * (pow(x-b,ss) - pow(x-y-b+r,ss)) * (pow(1-x,k-ss) - pow(1-x-a+b,k-ss))
                        for ss in range(min(para_m-2,k)+1)])
            temp_sum += temp * q_k[k]
        temp_sum_b = 0
        for k in range(k_max):
            temp = 1-pow(1-y,k)-pow(1-a,k)+pow(1-y-a+r,k)
            temp_sum_b += temp * q_k[k]
        b_new = q * a_new + (1 - q) * (temp_sum_b - temp_sum)

        temp_sum = 0
        for k in range(k_max):
            temp = sum([comb(k,ss) * ( pow(x,ss) * pow(1-x,k-ss) - pow(x-b,ss) * pow(1-x-a+b,k-ss)
                        - (1-q) * (pow(x-y,ss) * pow(1-x,k-ss) - pow(x-y-b+r,ss) * pow(1-x-a+b,k-ss)) )
                        for ss in range(min(para_m-1,k)+1)])
            temp_sum += temp * q_k[k]
        temp_sum_r = 0
        for k in range(k_max):
            temp = 1 - pow(1-a,k) - (1-q) * pow(1-y,k) + (1-q) * pow(1-y-a+r,k)
            temp_sum_r += temp * q_k[k]
        r_new = temp_sum_r - temp_sum

        delta = max(abs(a-a_new), abs(b-b_new), abs(r-r_new))
    return a_new, b_new, r_new


def get_P_inf_2(q,x,y,a,b,r,p_k,para_c,para_m):
    temp_sum = 0
    for k in range(1,k_max + 1):
        temp = sum([comb(k,ss) * pow(x-b,ss) * (pow(1-x,k-ss) - pow(1-x-a+b,k-ss)) for ss in range(min(para_m-1,k)+1)])
        temp_sum += temp * p_k[k]

    temp_sum_a = 0
    for k in range(1,k_max + 1):
        temp = 1 - pow(1-a,k)
        temp_sum_a += temp * p_k[k]
    a_new = temp_sum_a - temp_sum

    temp_sum_b = 0
    for k in range(1,k_max + 1):
        temp = 1 - pow(1 - y, k) - pow(1 - a, k) + pow(1 - y - a + r, k)
        temp_sum_b += temp * p_k[k]

    temp_sum = 0
    for k in range(1,k_max + 1):
        temp = sum([comb(k,ss) * (pow(x-b,ss) - pow(x-y-b+r,ss)) * (pow(1-x,k-ss) - pow(1-x-a+b,k-ss))
                    for ss in range(min(para_m-1,k)+1)])
        temp_sum += temp * p_k[k]

    P_inf = q * a_new + (1-q) * (temp_sum_b - temp_sum)
    return P_inf

def calc_un_SF_gOUT(qq,gamma_x,mm):
    ## qq 是比例为1的用户数
    ## mm 是渗流指数
    ## gamma_x 是SF的度分布指数
    pk, qk, average_k = function_general.gen_powerlaw_distribution(k_min, k_max, gamma_x)
    print("averagek"+str(average_k))
    x,y = get_xy(qq, pk, qk, mm)
    print("x,y"+str(x)+str(y))
    alpha_1, beta_1, gamma_1 = get_abr(qq, x, y, pk, qk, mm)
    P_GOUT = get_P_inf_2(qq, x, y, alpha_1, beta_1, gamma_1, pk, average_k, mm)
    print(P_GOUT)
    return P_GOUT,average_k

if __name__=="__main__":
    mm = 2
    gamma_x = 3.5
    for qq in np.arange(0.01, 0.5, 0.01):
        calc_un_SF_gOUT(qq,gamma_x,mm)