import numpy as np
import function_general
import networkx as nx
import json
import time
from scipy.special import zeta
k_min = 1
k_max = 100
n_node = 1000000
import function_general as func

def power_law_distribution(k_min, k_max, n_node, gamma):
    pk, qk, average_k = func.gen_powerlaw_distribution(k_min,k_max,gamma)
    degrees_sequence =  np.random.choice(range(k_min, k_max + 1), n_node,
                               p=[((i ** -gamma)/zeta(gamma)) / sum((np.arange(k_min, k_max + 1) ** -gamma/zeta(gamma))) for i in pk])
    return degrees_sequence

def gen_direct_SF(gammax):
    degrees_in = power_law_distribution(k_min, k_max, n_node, gammax)
    # Compute the total sum of the first distribution
    degrees_out = degrees_in.copy()
    np.random.shuffle(degrees_out)
    G = nx.directed_configuration_model(degrees_in,degrees_out,create_using=nx.DiGraph,seed=0)
    G.remove_edges_from(nx.selfloop_edges(G))
    return G

def get_gout(g, state):
    index = np.where(state==0)
    g.remove_nodes_from(index[0])
    g_list = sorted(list(nx.strongly_connected_components(g)), key=len, reverse=True)
    gSCC = list(g_list[0])
    gOUT_list = list(nx.descendants(g, source=gSCC[0]))
    P_GOUT = len(gOUT_list) / n_node
    return P_GOUT

def get_second_gout(g, state):
    index = np.where(state==0)
    g.remove_nodes_from(index[0])
    g_list = sorted(list(nx.strongly_connected_components(g)), key=len, reverse=True)
    gSCC = list(g_list[1])
    gOUT_list = list(nx.descendants(g, source=gSCC[0]))
    P_GOUT = len(gOUT_list) / n_node
    return P_GOUT

def run_simulation(g, m,initialization_ran,qq):
    state = np.zeros(n_node, dtype=int)
    index = np.where(initialization_ran <= qq)
    state[index[0]] = 1
    m_each_node = np.zeros(n_node, dtype=int) # m[i]: the number of state 1 in-neighbors of node i
    for node_now in g.nodes():
        neis = list(g.predecessors(node_now))
        m_each_node[node_now] = sum([state[nei] for nei in neis])
    new_strong_nodes = np.where( (m_each_node>=m) & (state==1))
    new_strong_nodes = list(new_strong_nodes[0])

    flag = 1
    while flag:
        flag = 0
        new_activite_nodes = list()
        for node_now in new_strong_nodes:
            neis = list(g.successors(node_now))
            for nei in neis:
                if state[nei] == 0:
                    state[nei] = 1
                    new_activite_nodes.append(nei)
                    flag = 1
        
        new_strong_nodes = list()
        for node_now in new_activite_nodes:
            if m_each_node[node_now] >= m:
                new_strong_nodes.append(node_now)
            neis = list(g.successors(node_now))
            for nei in neis:
                m_each_node[nei] += 1
                if (m_each_node[nei] == m) & (state[nei] == 1):
                    new_strong_nodes.append(nei)
                    flag = 1
    return get_gout(g, state)

def run_second_simulation(g, m, initialization_ran, qq):
    state = np.zeros(n_node, dtype=int)
    index = np.where(initialization_ran <= qq)
    state[index[0]] = 1
    m_each_node = np.zeros(n_node, dtype=int)  # m[i]: the number of state 1 in-neighbors of node i
    for node_now in g.nodes():
        neis = list(g.predecessors(node_now))
        m_each_node[node_now] = sum([state[nei] for nei in neis])
    new_strong_nodes = np.where((m_each_node >= m) & (state == 1))
    new_strong_nodes = list(new_strong_nodes[0])

    flag = 1
    loop_num = 0
    while flag:
        loop_num = loop_num+1
        flag = 0
        new_activite_nodes = list()
        for node_now in new_strong_nodes:
            neis = list(g.successors(node_now))
            for nei in neis:
                if state[nei] == 0:
                    state[nei] = 1
                    new_activite_nodes.append(nei)
                    flag = 1

        new_strong_nodes = list()
        for node_now in new_activite_nodes:
            if m_each_node[node_now] >= m:
                new_strong_nodes.append(node_now)
            neis = list(g.successors(node_now))
            for nei in neis:
                m_each_node[nei] += 1
                if (m_each_node[nei] == m) & (state[nei] == 1):
                    new_strong_nodes.append(nei)
                    flag = 1
    return get_gout(g, state),get_second_gout(g, state),loop_num

def scan_gamma_1(param,gamma_list):
    start_time = time.time()
    m = param
    for gamma in gamma_list:
        g = gen_direct_SF(gamma)
        initialization_ran = np.random.random(n_node)
        filename = '/Users/hongliangsun/Desktop/workspace/percolation/shl/code/20250414/data/SF_simu_'+str(m)+'_'+str(gamma)+'.txt'
        out_file = open(filename, 'w')
        for qq in np.arange(0.001, 0.2, 0.001):
            qq = round(qq,4)
            g_temp = g.copy()
            P_GOUT,P_second_GOUT,loop_num = run_second_simulation(g_temp,m,initialization_ran,qq)
            del g_temp
            print(qq, P_GOUT)
            out_file.write('%e\t%e\t%e\t%e\n' %(qq, P_GOUT,P_second_GOUT,loop_num))
            out_file.flush()
            end_time = time.time()
            print(end_time - start_time)

def scan_q(para_m,gamma_list):
    m = para_m
    out_file_path="result_simu.txt"
    for gamma in gamma_list:
        g = gen_direct_SF(gamma)
        initialization_ran = np.random.random(n_node)
        out_file = open(out_file_path, 'w')
        for qq in np.arange(0.001, 0.02, 0.001):
            qq = round(qq,4)
            g_temp = g.copy()
            P_GOUT = run_simulation(g_temp, m,initialization_ran,qq)
            print(qq, P_GOUT)
            out_file.write('%e\t%e\n' %(qq, P_GOUT))
            out_file.flush()

def scan_gamma(para_m,qq_list):
    out_file_path="result_simu.txt"
    for qq in qq_list:
        qq = round(qq, 4)
        initialization_ran = np.random.random(n_node)
        out_file = open(out_file_path, 'w')
        for gamma in np.arange(1.1, 3.5, 0.001):
            g = gen_direct_SF(gamma)
            g_temp = g.copy()
            P_GOUT = run_simulation(g_temp,para_m,initialization_ran,qq)
            print(qq, P_GOUT)
            out_file.write('%e\t%e\n' %(qq, P_GOUT))
            out_file.flush()

def test_gamma():
    m= 2
    gamma_list=[1.9069,1.9621,2.0379,2.0655,2.0793,2.0931]
    m = 4
    gamma_list = [1.73,1.7962,1.8872,1.9286,1.9617,2.0055]
    m = 5
    gamma_list = [1.71,1.7845,1.8176,1.8838,1.9252,1.9417]
    m = 6
    gamma_list = [1.6,1.7061,1.7551,1.8531,1.9102,1.9592]

def test_q():
    m = 2
    q_list=[0.0015,0.0029,0.0034,0.0054,0.0156]
    m = 4
    q_list = [0.0015,0.0039,0.0049,0.0088,0.0156]
    m = 5
    q_list = [0.0015,0.0044,0.0059,0.0103,0.0156]
    m = 6
    q_list = [0.0025,0.0049,0.0059,0.0112,0.0156]

def test_gamma_second_GOUT():
    param_m = 3
    gamma_list = [1.7]
    scan_gamma_1(param_m,gamma_list)

if __name__=="__main__":
    test_gamma_second_GOUT()
