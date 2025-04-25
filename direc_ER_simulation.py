import time
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import function_general
import math
k_max=100
N=500000
###
###@author Hong-Liang SUN
###@date 20220105
###

def m_index(g,g_stat,i):
    m = 0
    j_list = list(g.predecessors(i))
    k_list = []
    m_list = []
    for j in j_list:
        tmp_k = 0
        if g_stat[j] == 1:
            k_list = list(g.predecessors(j))
            for k in k_list:
                if g_stat[k]==1:
                    tmp_k += 1
        m_list.append(tmp_k)
    if len(m_list)> 0:
        m = max(m_list)
    else:
        m = 0
    return int(m)

def g_active_num(g_stat):
    num = 0
    for (v,s) in g_stat.items():
        if s ==1:
            num +=1
    return num

def gen_possion_distribution(ave_k, k_max):
    p_k = {}
    p_k[0] = math.exp(-ave_k)
    for kk in range(1, k_max + 1):
        p_k[kk] = p_k[kk - 1] * ave_k / kk
    return p_k

def simulation(g,m,q,k):
    ##变量区
    n = N
    ## 状态初始化
    np.random.seed(123)
    active = np.random.choice(n, int(n * q), replace=False)
    g_stat = {}
    for v in g.nodes():
        g_stat[v] = 0
    for k in active:
        g_stat[k] = 1

    ## 诱导
    loop = 1
    while loop:
        num1 = g_active_num(g_stat)
        for v in range(0, len(g)):
            if g_stat[v] == 0:
                v_m = m_index(g,g_stat,v)
                if v_m >= m:
                    g_stat[v] = 1
        num2 = g_active_num(g_stat)
        loop = abs(num2 - num1)

    ## 计算Active Number的数量
    A = g_active_num(g_stat)

    for v in range(0, n):
        if g_stat[v] == 0:
            g.remove_node(v)
    ##
    g_list = sorted(list(nx.strongly_connected_components(g)), key=len, reverse=True)
    gSCC = list(g_list[0])
    gout_list = list(nx.descendants(g, source=gSCC[0]))
    GOUT = len(gout_list)/n
    return GOUT,A

def simulation_second(g,m,q,k):
    ##变量区
    n = N
    ## 状态初始化
    np.random.seed(123)
    active = np.random.choice(n, int(n * q), replace=False)
    g_stat = {}
    for v in g.nodes():
        g_stat[v] = 0
    for k in active:
        g_stat[k] = 1

    ## 诱导
    loop = 1
    loop_num = 0
    while loop:
        num1 = g_active_num(g_stat)
        for v in range(0, len(g)):
            if g_stat[v] == 0:
                v_m = m_index(g,g_stat,v)
                if v_m >= m:
                    g_stat[v] = 1
        num2 = g_active_num(g_stat)
        loop_num = loop_num + 1
        loop = abs(num2 - num1)

    ## 计算Active Number的数量
    A = g_active_num(g_stat)

    for v in range(0, n):
        if g_stat[v] == 0:
            g.remove_node(v)
    ##
    g_list = sorted(list(nx.strongly_connected_components(g)), key=len, reverse=True)
    gOUT_1 = list(g_list[0])
    gSCC_2 = list(g_list[1])
    gout_list_1 = list(nx.descendants(g, source=gOUT_1[0]))
    #gout_list_2 = list(nx.descendants(g, source=gSCC_2[0]))

    GOUT_1 = len(gout_list_1)/n
    GSCC_2 = len(gSCC_2)/n
    return GOUT_1,GSCC_2,loop_num

def simmulation_ER_config_model(m,q,avg_k):
    n = N
    pk= gen_possion_distribution(avg_k, k_max)
    llist = []
    sum1 = 0
    sum2 = 0
    for k in pk:
        num = round(n * pk[k])
        sum1 += num
        for i in range(num):
            llist.append(k)
            sum2 += k
    if (sum2 % 2 == 1):
        llist[0] = 2

    g = nx.directed_configuration_model(llist, llist)
    g.remove_edges_from(nx.selfloop_edges(g))
    n = len(g)
    GOUT,A = simulation(g,m,q,k)
    return GOUT

def simulation_fast(n,m,q,k):
    p = k / (n - 1)
    g = nx.fast_gnp_random_graph(n,p,seed=123,directed=True)
    GOUT,A = simulation(g,m,q,k)
    return GOUT,A

def simulation_k(m_para, k_list):
    n=N
    m = m_para
    for avg_k in k_list:
        index_f = str(m) + '_' + str(avg_k) + '.txt'
        out_file_1 = open('/Users/hongliangsun/Desktop/workspace/percolation/shl/code/20250414/data/simu_k_' + index_f, 'w')
        p = avg_k / n
        g = nx.fast_gnp_random_graph(n, p, directed=True)
        for q in np.linspace(0.01, 0.2, num=100):
            gg = g.copy()
            GOUT = simulation(gg,m, q, avg_k)[0]
            print("average_k:" + str(avg_k))
            print("GOUT:" + str(GOUT))
            out_file_1.write('%f\t%f\n' % (q, GOUT))
            out_file_1.flush()
            del gg
        del g

def simulation_k_second(m_para, k_list):
    n=N
    m = m_para
    for avg_k in k_list:
        index_f = str(m) + '_' + str(avg_k) + '.txt'
        out_file_1 = open('/Users/hongliangsun/Desktop/workspace/percolation/shl/code/20250414/data/simu_k_' + index_f, 'w')
        p = avg_k / n
        g = nx.fast_gnp_random_graph(n, p,seed=123, directed=True)
        for q in np.linspace(0.1, 0.4, num=100):
            gg = g.copy()
            GOUT1,GOUT2,loop_num = simulation_second(gg,m, q, avg_k)
            print("average_k:" + str(avg_k))
            print("GOUT1:" + str(GOUT1))
            print("GOUT2:" + str(GOUT2))
            print("loop number:" + str(loop_num))
            #if GOUT > threshold:
            out_file_1.write('%f\t%f\t%f\t%f\n' % (q, GOUT1,GOUT2,loop_num))
            out_file_1.flush()
            del gg
        del g

def simulation_q_second(m_para, q_list):
    n=N
    m = m_para
    for q in q_list:
        index_f = str(m) + '_' + str(q) + '.txt'
        out_file_1 = open('/Users/hongliangsun/Desktop/workspace/percolation/shl/code/20250414/data/simu_q_' + index_f, 'w')
        for avg_k in np.linspace(1, 5, num=100):
            p = avg_k / n
            gg = nx.fast_gnp_random_graph(n, p,seed=123, directed=True)
            GOUT1,GOUT2,loop_num = simulation_second(gg,m, q, avg_k)
            print("average_k:" + str(avg_k))
            print("GOUT1:" + str(GOUT1))
            print("GOUT2:" + str(GOUT2))
            print("loop number:" + str(loop_num))
            #if GOUT > threshold:
            out_file_1.write('%f\t%f\t%f\t%f\n' % (avg_k, GOUT1,GOUT2,loop_num))
            out_file_1.flush()

def simulation_q(m_para, q_list):
    n=N
    m = m_para
    for qq in q_list:
        index_f = str(m) + '_' + str(qq) + '.txt'
        out_file_1 = open('/home/hlsun/percolation/data/simu_' + index_f, 'w')
        for avg_k in np.linspace(1.01, 10, num=900):
            p =  avg_k/ n
            g = nx.fast_gnp_random_graph(n, p, directed=True,seed=123)
            GOUT = simulation(g,m, qq, avg_k)[0]
            print("qq:" + str(avg_k))
            print("GOUT:" + str(GOUT))
            out_file_1.write('%f\t%f\n' % (avg_k,GOUT))
            out_file_1.flush()
            del g

def test_0430():
    m_para=2
    k_list=[1.5,1.9,2.0,2.1141,2.3]
    simulation_k(m_para, k_list)
    m_para = 4
    k_list = [2.2, 2.7, 3.4, 3.542, 3.7]
    simulation_k(m_para, k_list)
    m_para = 5
    k_list = [3.4, 3.7, 4.1, 4.244, 4.7]
    simulation_k(m_para, k_list)
    m_para = 6
    k_list = [3.6, 4.1, 4.7, 4.956, 5.6]
    simulation_k(m_para, k_list)

if __name__=="__main__":
    test_0430()