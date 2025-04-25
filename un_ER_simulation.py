import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

###
###@author Hong-Liang SUN
###@date 2023 0610
###

def m_index(g,g_stat,i):
    m = 0
    j_list = list(g.neighbors(i))
    k_list = []
    m_list = []
    for j in j_list:
        tmp_k = 0
        if g_stat[j] == 1:
            k_list = list(g.neighbors(j))
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

def test_m_index():
    n = 10
    p = 0.4
    q = 0.4
    # m = n*(n-1)*p/2
    # g = nx.gnm_random_graph(n,m,directed=True)
    g = nx.erdos_renyi_graph(n, p, directed=False)
    # g = g1.to_directed()
    # m = g.number_of_edges()
    # print("m:"+str(m))
    ## 状态初始化

    active = np.random.choice(n, int(n * q), replace=False)
    g_stat = {}
    for v in g.nodes():
        g_stat[v] = 0
    for k in active:
        g_stat[k] = 1
    for v in g.nodes():
        mm = m_index(g,g_stat,v)
        print(str(mm))
    node_labels = dict([(v,g_stat[v]) for v in g.nodes()])
    position = nx.circular_layout(g)
    nx.draw_networkx_labels(g,pos=position,labels=node_labels)
    nx.draw_networkx(g,pos=position,with_labels=False)
    plt.show()
    return

def simulation(m,q,k):
    ##变量区
    n = 500000
    p = k/(n-1)
    g = nx.fast_gnp_random_graph(n,p,directed=False)
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

    g_list = sorted(list(nx.connected_components(g)), key=len, reverse=True)
    gSCC = list(g_list[0])
    GOUT = len(gSCC)/n
    return GOUT,A

def test_simulation():
    q_list=[0.282]
    for m in range(3,4):
        for qq in q_list:
            q = round(qq,4)
            index_f = str(m) + '_' + str(q) + '.txt'
            out_file_1 = open('/Users/hongliangsun/Desktop/workspace/percolation/shl/data/simu_u4/simu_' + index_f, 'w')
            for average_k in np.linspace(1, 9.9901, num=300):
                GOUT = simulation(m, q, average_k)[0]
                print("average_k:" + str(average_k))
                print("GOUT:" + str(GOUT))
                if average_k < 10:
                    out_file_1.write('%f\t%f\n' % (average_k, GOUT))
                    out_file_1.flush()
                else:
                    out_file_1.write('%f\t%f' % (average_k, GOUT))
                    out_file_1.flush()

def simulation_second(m,q,k):
    ##变量区
    n = 500000
    p = k/(n-1)
    g = nx.fast_gnp_random_graph(n,p,seed=123,directed=False)
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
        loop_num = loop_num +1
        loop = abs(num2 - num1)

    ## 计算Active Number的数量
    A = g_active_num(g_stat)

    for v in range(0, n):
        if g_stat[v] == 0:
            g.remove_node(v)

    g_list = sorted(list(nx.connected_components(g)), key=len, reverse=True)
    #gSCC = list(g_list[0])
    if len(g_list)>1:
        gSCC_1 = list(g_list[0])
        gSCC_2 = list(g_list[1])

        GSCC_1 = len(gSCC_1)/n
        GSCC_2 = len(gSCC_2)/n
    else:
        GSCC_1 = 0
        GSCC_2 = 0
    return GSCC_1,GSCC_2,loop_num

def test_simulation_second_q():
    q_list=[0.19]
    for m in range(3,4):
        for qq in q_list:
            q = round(qq,4)
            index_f = str(m) + '_' + str(q) + '.txt'
            out_file_1 = open('/Users/hongliangsun/Desktop/workspace/percolation/shl/code/20250414/data/un_simu_' + index_f, 'w')
            for avg_k in np.linspace(1, 5, num=100):
                GSCC1,GSCC2,loop_num = simulation_second(m, qq, avg_k)
                print("average_k:" + str(avg_k))
                print("GSCC1:" + str(GSCC1))
                print("GSCC2:" + str(GSCC2))
                print("loop_num:" + str(loop_num))
                if avg_k < 10:
                    out_file_1.write('%f\t%f\t%f\t%f\n' % (avg_k, GSCC1,GSCC2,loop_num))
                    out_file_1.flush()
                else:
                    out_file_1.write('%f\t%f\t%f\t%f\n' % (avg_k, GSCC1,GSCC2,loop_num))
                    out_file_1.flush()


if __name__=="__main__":
    test_simulation_second_q()