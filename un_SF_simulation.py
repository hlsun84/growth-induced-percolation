import numpy as np
import function_general
import networkx as nx
import json
from scipy.special import zeta
import time

k_min = 1
k_max = 100
n_node = 5000000

def gen_un_SF(gamma):
    np.random.seed(123)
    degrees = np.random.choice(range(k_min, k_max + 1), n_node,
                               p=[(i ** -gamma)/zeta(gamma)/ sum((np.arange(k_min, k_max + 1) ** -gamma/zeta(gamma))) for i in range(k_min, k_max + 1)])
    if sum(degrees) % 2 == 1:
        degrees[0] += 1
    G = nx.configuration_model(degrees, create_using=nx.Graph,seed=0)
    G.remove_edges_from(nx.selfloop_edges(G))
    return G

def get_gout(g, state):
    index = np.where(state==0)
    g.remove_nodes_from(index[0])
    g_list = sorted(list(nx.connected_components(g)), key=len, reverse=True)
    gSCC = len(list(g_list[0]))/n_node
    return gSCC

def get_second_gout(g, state):
    index = np.where(state==0)
    g.remove_nodes_from(index[0])
    g_list = sorted(list(nx.connected_components(g)), key=len, reverse=True)
    gSCC = len(list(g_list[1]))/n_node
    return gSCC

def run_simulation(g, state, m):
    m_each_node = np.zeros(n_node, dtype=int) # m[i]: the number of state 1 neighbors of node i
    for node_now in range(n_node):
        neis = list(g.neighbors(node_now))
        m_each_node[node_now] = sum([state[nei] for nei in neis])
    new_strong_nodes = np.where( (m_each_node>=m) & (state==1))
    new_strong_nodes = list(new_strong_nodes[0])

    flag = 1
    while flag:
        flag = 0
        new_activite_nodes = list()
        for node_now in new_strong_nodes:
            neis = list(g.neighbors(node_now))
            for nei in neis:
                if state[nei] == 0:
                    state[nei] = 1
                    new_activite_nodes.append(nei)
                    flag = 1
        
        new_strong_nodes = list()
        for node_now in new_activite_nodes:
            if m_each_node[node_now] >= m:
                new_strong_nodes.append(node_now)
            neis = list(g.neighbors(node_now))
            for nei in neis:
                m_each_node[nei] += 1
                if (m_each_node[nei] == m) & (state[nei] == 1):
                    new_strong_nodes.append(nei)
                    flag = 1

    return get_gout(g, state)


def run_second_simulation(g, m,initialization_ran,qq):
    state = np.zeros(n_node, dtype=int)
    index = np.where(initialization_ran <= qq)
    state[index[0]] = 1
    m_each_node = np.zeros(n_node, dtype=int)  # m[i]: the number of state 1 neighbors of node i
    for node_now in range(n_node):
        neis = list(g.neighbors(node_now))
        m_each_node[node_now] = sum([state[nei] for nei in neis])
    new_strong_nodes = np.where((m_each_node >= m) & (state == 1))
    new_strong_nodes = list(new_strong_nodes[0])

    flag = 1
    loop_num = 0
    while flag:
        loop_num = loop_num + 1
        flag = 0
        new_activite_nodes = list()
        for node_now in new_strong_nodes:
            neis = list(g.neighbors(node_now))
            for nei in neis:
                if state[nei] == 0:
                    state[nei] = 1
                    new_activite_nodes.append(nei)
                    flag = 1

        new_strong_nodes = list()
        for node_now in new_activite_nodes:
            if m_each_node[node_now] >= m:
                new_strong_nodes.append(node_now)
            neis = list(g.neighbors(node_now))
            for nei in neis:
                m_each_node[nei] += 1
                if (m_each_node[nei] == m) & (state[nei] == 1):
                    new_strong_nodes.append(nei)
                    flag = 1

    return get_gout(g, state),get_second_gout(g,state),loop_num

def scan_gamma(param,gamma_list):
    start_time = time.time()
    m = param
    for gamma in gamma_list:
        g = gen_un_SF(gamma)
        initialization_ran = np.random.random(n_node)
        filename = '/Users/hongliangsun/Desktop/workspace/percolation/shl/code/20250414/data/un_big_SF_simu_'+str(m)+"_"+str(gamma)+'.txt'
        out_file = open(filename, 'w')
        for qq in np.arange(0.001, 0.1, 0.001):
            qq = round(qq,4)
            g_temp = g.copy()
            P_GOUT,P_second_GOUT,loop_num = run_second_simulation(g_temp,m,initialization_ran,qq)
            del g_temp
            print(qq, P_GOUT)
            out_file.write('%e\t%e\t%e\t%e\n' %(qq, P_GOUT,P_second_GOUT,loop_num))
            out_file.flush()
            end_time = time.time()
            print(end_time - start_time)

def scan_q():
    m = 2
    gamma = 2.5
    g = gen_un_SF(gamma)
    initialization_ran = np.random.random(n_node)
    out_file = open('/Users/hongliangsun/Desktop/workspace/percolation/shl/data/SF0507/simu_xie/result.txt', 'w')
    for qq in np.arange(0.001, 0.02, 0.001):
        state = np.zeros(n_node, dtype=int)
        index = np.where(initialization_ran <= qq)
        state[index[0]] = 1
        g_temp = g.copy()
        gSCC = run_simulation(g_temp, state, m)
        P_GOUT = gSCC / n_node
        print(qq, P_GOUT)
        out_file.write('%e\t%e\n' %(qq, P_GOUT))

def test_gamma_second_GOUT():
    param_m = 6
    gamma_list = [2.8]
    scan_gamma(param_m,gamma_list)

if __name__=="__main__":
    test_gamma_second_GOUT()