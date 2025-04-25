import time
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.special import comb
import function_general
import networkx as nx
import os

Accuracy = 1e-12
k_min = 1
k_max = 100

###
###@author Hong-Liang SUN
###@date 20220105
###
def calc_XYA(avg_k, p_k, m, q):
    x2 = q
    y2 = 0.0
    delta = 1.0
    while delta > Accuracy:
        x1 = x2
        y1 = y2
        x2 = 1-(1-q)*sum([p_k[k] * (1-y1)**k for k in p_k])
        dxy = sum([p_k[k] * sum([comb(k, s) * ((x1 ** s - (1-q)*(x1-y1)**s))*(1-x1)**(k-s) for s in range(min(m, k+1))]) for k in p_k])
        y2 = x2 - dxy
        A = x2
        delta = abs(x1 - x2) + abs(y1 - y2)
    return A,x2,y2

def calc_XINF_YINF_PINF(x2,y2,avg_k,p_k,m,q):
    x_inf2 = x2 ###不知道
    y_inf2 = y2
    delta = 1.0
    while delta>Accuracy:
        x_inf1 = x_inf2
        y_inf1 = y_inf2
        sum1 = sum([p_k[k]*(1-x_inf1)**k for k in p_k])
        sum2 = sum([p_k[k]*((1 - y2) ** k - (1 - y2 - x_inf1 + y_inf1) ** k) for k in p_k])
        x_inf2 = 1-sum1-(1-q)*sum2
        sum3 = sum([p_k[k]  * sum([
                comb(k, s) * (1 - x2) ** (k - s) * (x2 ** s - (x2 - x_inf1) ** s) for s in range(min(m, k + 1))]) for k in p_k])
        sum4 = sum([p_k[k] * sum([
                comb(k, s) * (1 - x2) ** (k - s) * ((x2 - y2) ** s - (x2 - y2 - x_inf1 + y_inf1) ** s) for s in
                range(min(m, k + 1))]) for k in p_k])
        dxy = sum3 - (1-q)*sum4
        y_inf2 = x_inf2-dxy
        p_inf = x_inf2
        delta = abs(x_inf2-x_inf1)+abs(y_inf2-y_inf1)
    return p_inf,x_inf2,y_inf2

def calc_SF_gOUT(qq,gamma_x,mm):
    ## qq 是比例为1的用户数
    ## mm 是渗流指数
    ## gamma_x 是SF的度分布指数
    print("kmin"+str(k_min))
    print("kmax" + str(k_max))
    print("gammax" + str(gamma_x))
    pk, qk, average_k = function_general.gen_powerlaw_distribution(k_min, k_max, gamma_x)
    A, x, y = calc_XYA(average_k, pk, mm, qq)
    P_GOUT, x_inf, y_inf = calc_XINF_YINF_PINF(x, y, average_k, pk, mm, qq)
    return P_GOUT

def calc_SF_gOUT_q(gamma_x,mm):
    ## mm 是渗流指数
    ## gamma_x 是SF的度分布指数
    #print("kmin"+str(k_min))
    #print("kmax" + str(k_max))
    #print("gammax" + str(gamma_x))
    GOUT_ret=0
    for qq in np.linspace(1e-3,2e-2,num=200):
        pk, qk, average_k = function_general.gen_powerlaw_distribution(k_min, k_max, gamma_x)
        A, x, y = calc_XYA(average_k, pk, mm, qq)
        P_GOUT, x_inf, y_inf = calc_XINF_YINF_PINF(x, y, average_k, pk, mm, qq)
        if P_GOUT>1e-6:
            GOUT_ret = P_GOUT
            break
    return GOUT_ret

def incuded_q_0513(para_m,gamma_list):
    for mm in range(2,11,2):
        for gamma in gamma_list:
            gamma = round(gamma,4)
            print("gamma:" + str(gamma))
            index_f = str(mm)+str('_')+str(gamma)+str('.txt')
            out_file_2 = open('/Users/hongliangsun/Desktop/workspace/percolation/shl/data/SF-1008/numerical/GOUT_'+index_f, 'w')
            kk_min = 1
            kk_max = 100
            for qq in np.arange(0.001,0.2,0.001):
                qq=round(qq,4)
                pk, qk, average_k = function_general.gen_powerlaw_distribution(kk_min, kk_max, gamma)
                print("qq:" + str(qq))
                print("average_k:"+str(average_k))
                A,x,y = calc_XYA(average_k, pk, mm, qq)
                P_GOUT,x_inf,y_inf = calc_XINF_YINF_PINF(x,y,average_k,pk,mm,qq)
                print("GOUT:"+str(P_GOUT))
                out_file_2.write('%f\t%.8f\t%.8f\n' % (average_k,qq,P_GOUT))
                out_file_2.flush()

if __name__=="__main__":
    gamma_list=[1.8]
    para_m = 2
    incuded_q_0513(para_m,gamma_list)
