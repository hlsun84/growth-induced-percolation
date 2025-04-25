import time
import numpy as np
import math
from scipy import stats
import matplotlib.pyplot as plt
from scipy.special import comb

Accuracy = 1e-14
k_max = 100

###
###@author Hong-Liang SUN
###@date 20240105
###

def selfXY(ave_k, p_k, m, q):
    x2 = q
    y2 = 0.0
    delta = 1.0
    while delta > Accuracy:
        x1 = x2
        y1 = y2
        x2 = 1 - (1-q) * math.exp(- ave_k * y1)
        dxy = 0.0
        for k in p_k:
            temp_2 = sum(comb(k,s) * (1-x1)**(k-s) * (x1**s - (1-q)*(x1-y1)**s) for s in range(min(m,k+1)))
            dxy += p_k[k] * temp_2
        y2 = x2 - dxy
        delta = abs(x1 - x2) + abs(y1 - y2)
    A = x2
    return A, x2, y2

def gen_possion_distribution(ave_k, k_max):
    p_k = {}
    p_k[0] = math.exp(-ave_k)
    for kk in range(1, k_max + 1):
        p_k[kk] = p_k[kk - 1] * ave_k / kk
    return p_k

def calc_XYA(avg_k, p_k, m, q):
    x2 = q
    y2 = 0.0
    delta = 1.0
    while delta>Accuracy:
        x1 = x2
        y1 = y2
        x2 = 1 - (1 - q) * math.exp(- avg_k * y1)
        dxy = 0.0
        for k in p_k:
            sum1 = sum( comb(k,s)*((x1**s-(1-q)*(x1-y1)**s))*((1-x1)**(k-s))  for s in range(min(m,k+1)))
            dxy += p_k[k]*sum1
        y2 = x2 -dxy
        delta = abs(x1-x2)+abs(y1-y2)
    A = x2
    return A,x2,y2

def calc_XINF_YINF_PINF(x2,y2,avg_k,p_k,m,q):

    x_inf2 = x2 ###不知道
    y_inf2 = y2
    delta = 1.0
    while delta>Accuracy:
        x_inf1 = x_inf2
        y_inf1 = y_inf2

        x_inf2 = 1- math.exp(-avg_k*x_inf1)-(1-q)*math.exp(-avg_k*y2)+(1-q)*math.exp(-avg_k*(y2+x_inf1-y_inf1))
        sum2 = 0.0
        sum1 = 0.0
        for k in p_k:
            sum1 = sum(comb(k,s)*((1-x2)**(k-s))*(x2**s - (x2-x_inf1)**s - (1-q)*(x2-y2)**s +(1-q)*(x2-y2-x_inf1 +y_inf1)**s) for s in range(min(m,k+1)))
            sum2 += p_k[k]*sum1
        y_inf2 = x_inf2-sum2
        delta = abs(x_inf2-x_inf1)+abs(y_inf2-y_inf1)
    p_inf = x_inf2
    return p_inf,x_inf2,y_inf2

def calc_gOUT(q,avg_k,m):
    ## q 是比例为1的用户数
    ## avg_k 是平均的度
    ## m 是渗流指数

    start_time = time.time()
    comb = math.comb
    pk = gen_possion_distribution(avg_k,k_max)
    A,x,y = calc_XYA(avg_k,pk,m,q)
    P_GOUT,x_inf,y_inf = calc_XINF_YINF_PINF(x,y,avg_k,pk,m,q)
    return A, P_GOUT

def test_incuded_k_m(m_para,avg_k_para):
    average_k = round(avg_k_para,4)
    index_f = str(m_para)+str('_')+str(average_k)+str('.txt')
    out_file_2 = open('/Users/hongliangsun/Desktop/workspace/percolation/shl/code/20250414/data/GOUT_'+index_f, 'w')
    for qq in np.linspace(0.01, 0.2, num=100):
        p_k = gen_possion_distribution(average_k, k_max)
        A,x,y = calc_XYA(average_k, p_k, m_para, qq)
        P_GOUT,x_inf,y_inf = calc_XINF_YINF_PINF(x,y,average_k,p_k,m_para,qq)
        if P_GOUT>0.999 and average_k>10.0:
            print(qq, P_GOUT)
            out_file_2.write('%f\t%.8f\t%.8f\n' % (average_k,qq,P_GOUT))
            out_file_2.flush()
            break
        else :
            print(qq, P_GOUT)
            out_file_2.write('%f\t%.8f\t%.8f\n' % (average_k,qq,P_GOUT))
            out_file_2.flush()

def test_incuded_q_m(m_para,q_para):
    q_para = round(q_para,4)
    index_f = str(m_para)+str('_')+str(q_para)+str('.txt')
    out_file_2 = open('/Users/hongliangsun/Desktop/workspace/percolation/shl/data/24-0323/direct4/GOUT_'+index_f, 'w')
    for avg_k in np.linspace(1.0, 10, num=900):
        p_k = gen_possion_distribution(avg_k, k_max)
        A,x,y = calc_XYA(avg_k, p_k, m_para, q_para)
        P_GOUT,x_inf,y_inf = calc_XINF_YINF_PINF(x,y,avg_k,p_k,m_para,q_para)
        print(q_para, P_GOUT)
        out_file_2.write('%f\t%.8f\t%.8f\n' % (avg_k,q_para,P_GOUT))
        out_file_2.flush()

def test_induced_k_test():
    m_para=2
    avg_k_para=2.115
    test_incuded_k_m(m_para,avg_k_para)
    m_para = 3
    avg_k_para = 2.843
    test_incuded_k_m(m_para, avg_k_para)
    m_para = 4
    avg_k_para = 3.544
    test_incuded_k_m(m_para, avg_k_para)
    return

def test_0320():
    k=2.3
    m=2
    test_incuded_k_m(m,k)

if __name__=="__main__":
    test_0320()