import time
import numpy as np
import math
from scipy import stats
import matplotlib.pyplot as plt
from scipy.special import comb
import numpy

Accuracy = 1e-14
k_max = 100
###
###@author Hong-Liang SUN
###@date 20220117
###

def gen_possion_distribution(ave_k, k_max):
    p_k = {}
    p_k[1] = math.exp(-ave_k)
    for kk in range(2, k_max + 1):
        p_k[kk] = p_k[kk - 1] * ave_k / kk
    return p_k

def func_AAA(x,s,k):
    sum =  comb(k-1,s)*(x**s)*(1-x)**(k-1-s)
    return sum

def func_BBB(x,y,t,s):
    sum = comb(s,t)*(y/x)**(t)*(1-y/x)**(s-t)
    return sum

def selfXY(ave_k, p_k, m, q):
    x2 = q
    y2 = 0.0
    delta = 1.0
    while delta > Accuracy:
        x1 = x2
        y1 = y2
        x2 = 1 - (1 - q) * sum([k*p_k[k]/ave_k*(1-y1)**(k-1) for k in p_k])
        summ = 0
        for k in p_k:
            sum1 = q* sum([func_AAA(x1,s,k) for s in range(m,k)])
            sum2 = (1-q)*sum([func_AAA(x1,s,k)*(1-(1-y1/x1)**s) for s in range(m,k)])
            summ += k*p_k[k]/ave_k*(sum1+sum2)
        y2 = summ
        delta = max(abs(y1-y2),abs(x1-x2))
    return x2,y2

def calc_alpha(x2,y2,alpha_2,beta_2,avg_k,p_k,m,q):
    temp_sum = 0.0
    for k in p_k:
        temp_sum_1 = sum([comb(k-1,s)*(1-x2)**(k-1-s)*(x2**s-(x2-beta_2)**s) for s in range(m-1)])
        temp_sum_2 = sum([comb(k-1,s)*(x2**s*(1-x2)**(k-1-s)-(x2-beta_2)**s*(1-x2-alpha_2+beta_2)**(k-1-s)) for s in range(min(m-1,k-1),k)])
        temp_sum += k*p_k[k]/avg_k * (temp_sum_1+temp_sum_2)
    return temp_sum

def calc_beta(x2,y2,alpha_2,beta_2,gamma_2,p_k,m,q,avg_k):
    temp_sum = 0.0
    for k in p_k:
        temp_sum1 = q*sum([comb(k-1,s)*(1-x2)**(k-1-s)*(x2**s-(x2-beta_2)**s) for s in range(m-1)])
        temp_sum2 = (1-q)*sum([comb(k-1,s)*(1-x2)**(k-1-s)*(x2**s-(x2-y2)**s+(x2-beta_2)**s-(x2-y2+beta_2-gamma_2)**s) for s in range(m-1)])
        temp_sum3 = q*sum([comb(k-1,s)*((x2**s)*(1-x2)**(k-1-s)-(x2-beta_2)**s*(1-x2-alpha_2+beta_2)**(k-1-s)) for s in range(min(m-1,k-1),k)])
        temp_sum4 = (1-q)*sum([comb(k-1,s)*(x2**s*(1-x2)**(k-1-s)-(x2-y2)**s*(1-x2)**(k-1-s)-(x2-beta_2)**s*(1-x2-alpha_2+beta_2)**(k-1-s)+(x2-y2+beta_2-gamma_2)**s*(1-x2-alpha_2+beta_2)**(k-1-s)) for s in range(min(m-1,k-1),k)])
        temp_sum += p_k[k]*k/avg_k*(temp_sum1 + temp_sum2 + temp_sum3 + temp_sum4)
    return temp_sum

def calc_gamma(x2,y2,alpha2,beta_2,gamma_2,p_k,m,q,avg_k):
    temp_sum = 0.0
    for k in p_k:
        temp_sum1 = q*sum([comb(k-1,s)*((x2**s*(1-x2)**(k-1-s))-(x2-beta_2)**s*(1-x2-alpha2+beta_2)**(k-1-s)) for s in range(min(m,k-1),k)])
        temp_sum2 = (1-q)*sum([(x2**s*(1-x2)**(k-1-s)-(x2-beta_2)**s*(1-x2-alpha2+beta_2)**(k-1-s)-(x2-y2)**s*(1-x2)**(k-1-s)+(x2-y2+beta_2-gamma_2)**s*(1-x2-alpha2+beta_2)**(k-1-s)) for s in range(min(m-1,k-1),k)])
        temp_sum += p_k[k]*k/avg_k*(temp_sum1+temp_sum2)
    return temp_sum

def calc_ABGTP(x2,y2,avg_k,p_k,m,q):
    alpha_2 = 1
    beta_2 = x2
    gamma_2 = y2
    delta = 1
    alpha_1 = alpha_2
    beta_1 = beta_2
    gamma_1 = gamma_2
    while delta > Accuracy:
        alpha_1 = alpha_2
        beta_1 = beta_2
        gamma_1 = gamma_2
        alpha_2 = calc_alpha(x2, y2, alpha_1, beta_1, avg_k, p_k, m, q)
        beta_2 = calc_beta(x2, y2, alpha_1,beta_1, gamma_1, p_k, m, q,avg_k)
        gamma_2 = calc_gamma(x2, y2,alpha_1, beta_1, gamma_1, p_k, m, q,avg_k)
        delta = max(abs(alpha_1-alpha_2),abs(beta_1-beta_2),abs(gamma_1-gamma_2))
    return alpha_2,beta_2,gamma_2

def cals_GOUT(x2, y2, alpha_1,beta_1, gamma_1, p_k, m, q):
    temp_sum = 0.0
    for k in p_k:
        temp_sum_1 = q*sum([comb(k,s)* (1-x2)**(k-s)*(x2**s-(x2-beta_1)**s) for s in range(m)])
        temp_sum_2 = (1-q)*sum ([comb(k,s)*(1-x2)**(k-s)*(x2**s-(x2-y2)**s-(x2-beta_1)**s+(x2-y2+beta_1-gamma_1)**s) for s in range(m)])
        temp_sum_3 = q*sum([comb(k,s)*(x2**s*(1-x2)**(k-s)-(x2-beta_1)**s*(1-x2-alpha_1+beta_1)**(k-s)) for s in range(min(m,k),k+1)])
        temp_sum_4 = (1-q)*sum([comb(k,s)*(x2**s*(1-x2)**(k-s)-(x2-y2)**s*(1-x2)**(k-s)-(x2-beta_1)**s*(1-x2-alpha_1+beta_1)**(k-s)+(x2-y2+beta_1-gamma_1)**s*(1-x2-alpha_1+beta_1)**(k-s)) for s in range(min(m,k),k+1)])
        temp_sum += p_k[k]*(temp_sum_1+temp_sum_2+temp_sum_3+temp_sum_4)
    return temp_sum

def test_induced_q():
    for mm in range(3,4):
        for qq in np.linspace(0.01,0.9901,num=99):
            qq = round(qq,4)
            index_f = str(mm)+str('_')+str(qq)+str('.txt')
            out_file_2 = open('/Users/hongliangsun/Desktop/workspace/percolation/shl/data/ER_GOUT_un_shl_0822/GOUT_'+index_f, 'w')
            for average_k in numpy.arange(3.6, 3.7, 0.01):
                p_k = gen_possion_distribution(average_k, k_max)
                x,y = selfXY(average_k, p_k, mm, qq)
                alpha_1, beta_1, gamma_1 = calc_ABGTP(x,y,average_k,p_k,mm,qq)
                P_GOUT = cals_GOUT(x, y, alpha_1, beta_1,gamma_1, p_k, mm, qq)
                print('==============================================')
                print(average_k,P_GOUT)
                out_file_2.write('%f\t%f\n' % (average_k,P_GOUT))
    return

def test_induced_avg_k():
    qq_list = [0.12,0.19,0.292,0.42,0.72]
    qq_list_2 = [0.282]
    for mm in range(3,4):
        for qq in qq_list_2:
            qq = round(qq,4)
            print('==============================================')
            print(qq)
            index_f = str(mm)+str('_')+str(qq)+str('.txt')
            out_file_2 = open('/Users/hongliangsun/Desktop/workspace/percolation/shl/data/ER_GOUT_un_shl_0911/GOUT_'+index_f, 'w')
            for average_k in numpy.linspace(1, 9.9901,num=300):
                average_k = round(average_k, 4)
                p_k = gen_possion_distribution(average_k, k_max)
                x,y = selfXY(average_k, p_k, mm, qq)
                alpha_1, beta_1, gamma_1 = calc_ABGTP(x,y,average_k,p_k,mm,qq)
                P_GOUT = cals_GOUT(x, y, alpha_1, beta_1,gamma_1, p_k, mm, qq)
                print('==============================================')
                print(average_k,P_GOUT)
                out_file_2.write('%f\t%f\n' % (average_k,P_GOUT))
                out_file_2.flush()
    return

if __name__=="__main__":
    test_induced_avg_k()