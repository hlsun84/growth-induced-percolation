import itertools
import math
import statistics
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols
import statsmodels.stats.anova as anova
import networkx as nx
import numpy as np

def preprocess():
    file_name="/Users/hongliangsun/Desktop/workspace/percolation/Carrer_burst/data/dblp/dblp-v5.txt"
    year = 0
    author_papers={}
    edge_list_list = []
    year_set = set()
    author_list = []
    with open(file_name,'r') as file:
        for line in file.readlines():
            stripped_line = line.strip()
            if len(stripped_line) != 0:
                if stripped_line.startswith('#year'):
                    year = int(stripped_line.replace('#year', ''))
                    year_set.add(year)
                if stripped_line.startswith('#@') and year >= 2002:
                    author_list = stripped_line.replace('#@', '').split(',')
                    edge_list_list.append([edge for edge in list(itertools.combinations(author_list, 2))])
                if stripped_line.startswith('#citation') and year >= 2002:
                    numbers = math.fabs(int(stripped_line.replace('#citation', '')))
                    for author in author_list:
                        if author not in author_papers:
                            author_papers[author] = {}
                            author_papers[author][year-2001]=[numbers]
                        else:
                            if author_papers.get(author).get(year-2001) != None:
                                author_papers.get(author).get(year-2001).append(numbers)
                            else:
                                author_papers[author][year-2001] = [numbers]
    sorted_years = sorted(year_set)
    print(sorted_years)
    return author_papers,sorted_years,edge_list_list

def count_consecutive_items(lists,thr):
    count = 0
    count_ls = []
    for item in lists:
        if item <= thr:
            if count>=1:
                count_ls.append(count)
                count = 0
            continue
        else:
            count += 1
    if count>0:
        count_ls.append(count)
    max_len = max(count_ls) if len(count_ls) > 0 else 0
    mean_ls = sum(count_ls)/len(count_ls) if len(count_ls) > 0 else 0
    return mean_ls,max_len


def calc_citations():
    author_papers, years, edge_list_list = preprocess()
    author_dict = {}
    line_num = 1
    author_citations = {}
    for author, total_papers in author_papers.items():
        paper_all = []
        papers_1 = []
        papers_2 = []
        citations = 2 * [0]
        for year, papers in total_papers.items():
            paper_all = paper_all + papers
            if year <= 4:
                papers_1 = papers_1 + papers
            if year > 4:
                papers_2 = papers_2 + papers

        if len(paper_all) >= 1 and sum(paper_all) >= 1:
            author_dict[author] = 1
        line_num = line_num + 1
        citations[0] = statistics.mean(papers_1) if len(papers_1) >= 1 and sum(papers_1) >= 1 else []
        citations[1] = statistics.mean(papers_2) if len(papers_2) >= 1 and sum(papers_2) >= 1 else []
        author_citations[author] = citations

    print("======")
    unGraph = nx.Graph()
    for edge_list in edge_list_list:
        for edge in edge_list:
            if author_dict.get(edge[0]) != None and author_dict.get(edge[1]) != None:
                unGraph.add_edge(edge[0], edge[1])
    n = len(unGraph.nodes())
    m = len(unGraph.edges())
    print('number of nodes:' + str(n))
    print('number of edges:' + str(m))
    return author_citations, unGraph

def calc_hot_streaks():
    author_papers,years,edge_list_list = preprocess()
    author_hot_streak = {}
    author_dict = {}
    line_num = 1

    for author,total_papers in author_papers.items():
        paper_all = []
        papers_1 = []
        papers_2 = []

        hot_streaks = 2*[0]

        for year,papers in total_papers.items():
            paper_all = paper_all + papers
            if year <= 4:
                papers_1 = papers_1 + papers
            if year > 4:
                papers_2 = papers_2 + papers

        median = statistics.median(paper_all) if len(paper_all)>=1 and sum(paper_all)>=1 else 0
        median1 = statistics.median(papers_1) if len(papers_1) >= 1 and sum(papers_1) >= 1 else 0
        median2 = statistics.median(papers_2) if len(papers_2) >= 1 and sum(papers_2) >= 1 else 0
        if len(paper_all) >= 1 and sum(paper_all) >= 1 :
            author_dict[author] = 1

        line_num = line_num + 1

        hot_streaks[0] = count_consecutive_items(papers_1,median) if len(papers_1) >= 1 and sum(papers_1) >= 1 else []
        hot_streaks[1] = count_consecutive_items(papers_2,median) if len(papers_2) >= 1 and sum(papers_2) >= 1 else []
        author_hot_streak[author] = hot_streaks

    print("======")
    unGraph = nx.Graph()
    for edge_list in edge_list_list:
        for edge in edge_list:
            if author_dict.get(edge[0]) != None and author_dict.get(edge[1]) != None:
                unGraph.add_edge(edge[0],edge[1])
    n = len(unGraph.nodes())
    m = len(unGraph.edges())
    print('number of nodes:' + str(n))
    print('number of edges:' + str(m))
    return author_hot_streak,unGraph

def get_g_stat(unGraph,author_hot_streaks):
    g_stat_1 = {}
    g_stat_2 = {}
    for node in unGraph.nodes():
        g_stat_1[node] = 1 if len(author_hot_streaks[node][0])!=0 and author_hot_streaks[node][0][1]>=1 else 0
        g_stat_2[node] = 1 if len(author_hot_streaks[node][1])!=0 and author_hot_streaks[node][1][1]>=1 else 0
    return g_stat_1,g_stat_2

def m_index(g,g_stat,i):
    m = 0
    j_list = list(g.neighbors(i))
    k_list = []
    m_list = []
    for j in j_list:
        tmp_k = 0
        if g_stat[j] == 1:
            k_list = list(g.neighbors(j))
            tmp_k = 0
            for k in k_list:
                if g_stat[k]==1 and k!=i:
                    tmp_k += 1
        m_list.append(tmp_k)
    if len(m_list)> 0:
        m = max(m_list)
    else:
        m = 0
    return int(m)

def calc_mIndx():
    author_hot_streak, unGraph = calc_hot_streaks()
    g_stat_1,g_stat_2 = get_g_stat(unGraph,author_hot_streak)
    m_indx1_dict = {}
    m_indx2_dict = {}
    for i in unGraph.nodes():
        m_indx1_dict[i] = m_index(unGraph,g_stat_1,i)
        m_indx2_dict[i] = m_index(unGraph,g_stat_2,i)
    #print(statistics.mean(m_indx1_dict.values()))
    return unGraph,g_stat_1,g_stat_2,m_indx1_dict,m_indx2_dict,author_hot_streak

def calc_mIndx_citations():
    author_citations, unGraph = calc_citations()
    author_hot_streak, unGraph = calc_hot_streaks()
    g_stat_1,g_stat_2 = get_g_stat(unGraph,author_hot_streak)
    m_indx1_dict = {}
    m_indx2_dict = {}
    for i in unGraph.nodes():
        m_indx1_dict[i] = m_index(unGraph,g_stat_1,i)
        m_indx2_dict[i] = m_index(unGraph,g_stat_2,i)
    return unGraph,g_stat_1,g_stat_2,m_indx1_dict,m_indx2_dict,author_citations

def calc_degree():
    author_hot_streak, unGraph = calc_hot_streaks()
    g_stat_1,g_stat_2 = get_g_stat(unGraph,author_hot_streak)
    deg_indx1_dict = {}
    deg_indx2_dict = {}
    for i in unGraph.nodes():
        deg_indx1_dict[i] = unGraph.degree[i]
        deg_indx2_dict[i] = unGraph.degree[i]
    return unGraph,g_stat_1,g_stat_2,deg_indx1_dict,deg_indx2_dict,author_hot_streak

def calc_degree_induce_m():
    author_hot_streak, unGraph = calc_hot_streaks()
    g_stat_1,g_stat_2 = get_g_stat(unGraph,author_hot_streak)
    deg_indx1_dict = {}
    deg_indx2_dict = {}
    m_indx1_dict = {}
    m_indx2_dict = {}
    for i in unGraph.nodes():
        deg_indx1_dict[i] = unGraph.degree[i]
        deg_indx2_dict[i] = unGraph.degree[i]
        m_indx1_dict[i] = m_index(unGraph, g_stat_1, i)
        m_indx2_dict[i] = m_index(unGraph, g_stat_2, i)
    return unGraph,g_stat_1,g_stat_2,deg_indx1_dict,deg_indx2_dict,m_indx1_dict,m_indx2_dict,author_hot_streak

def display_10(Hot_streak,x_label,y_label):
    weights = [
        (Hot_streak[0]),
        (Hot_streak[1]),
        (Hot_streak[2]),
        (Hot_streak[3]),
        (Hot_streak[4]),
        (Hot_streak[5]),
        (Hot_streak[6]),
        #(Hot_streak[7]),
        #(Hot_streak[8]),
        #(Hot_streak[9]),
    ]
    labels_deg = ['[2-21]', '[22-41]', '[42-61]','[62-81]','[82-101]','[102-121]','[122-141]','[142-161]','[162-181]','[182-194]']
    labels_induced = []
    colors = ['peachpuff', 'orange', 'tomato','red','peachpuff', 'orange', 'tomato','red','peachpuff', 'orange']

    fig, ax = plt.subplots()
    ax.set_ylabel(y_label)
    ax.set_xlabel(x_label)
    bplot = ax.boxplot(weights,
                       patch_artist=True,  # fill with color
                       labels=labels_deg,
                       showmeans=True,
                       showfliers=False)  # will be used to label x-ticks
    # fill with colors
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
    plt.show()

def display_10_indirect(Hot_streak,x_label,y_label):
    weights = [
        (Hot_streak[0]),
        (Hot_streak[1]),
        (Hot_streak[2]),
        (Hot_streak[3]),
        (Hot_streak[4]),
        (Hot_streak[5]),
        (Hot_streak[6])
    ]
    labels_deg = ['[1-90]', '[91-180]', '[181-270]','[271-360]','[361-450]','[451-540]','[541-893]']
    labels_induced = []
    colors = ['peachpuff', 'orange', 'tomato','red','peachpuff', 'orange', 'tomato','red','peachpuff', 'orange']

    fig, ax = plt.subplots()
    ax.set_ylabel(y_label)
    ax.set_xlabel(x_label)

    bplot = ax.boxplot(weights,
                       patch_artist=True,  # fill with color
                       labels=labels_deg,
                       showmeans=True,
                       showfliers=False)  # will be used to label x-ticks
    # fill with colors
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
    plt.show()

def display_5(Hot_streak,x_label,y_label):
    weights = [
        (Hot_streak[0]),
        (Hot_streak[1]),
        (Hot_streak[2]),
        (Hot_streak[3]),
        (Hot_streak[4]),
    ]
    labels = ['1-25','26-50','51-75','76-100','101-127']
    colors = ['peachpuff', 'orange', 'tomato','red','peachpuff']

    fig, ax = plt.subplots()
    ax.set_ylabel(y_label)
    ax.set_xlabel(x_label)

    bplot = ax.boxplot(weights,
                       patch_artist=True,  # fill with color
                       labels=labels,
                       showmeans=True,
                       showfliers=False)  # will be used to label x-ticks
    # fill with colors
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
    plt.show()

def display_errorbar(x,y,yerr):
    plt.errorbar(x, y, yerr=yerr, fmt='o', ecolor='r', capsize=5, capthick=2, markeredgecolor='blue',
                 markerfacecolor='green')
    plt.xlabel('X-axis label')
    plt.ylabel('Y-axis label')
    plt.title('Plot with Error Bars')
    plt.show()

def display_errorbar_demo(segments):
    x = [1 for i in len(segments[0])]
    y = segments[0]
    yerr = statistics.stdev(segments[0])  # error bar values (constant in this case)

    # Plotting the data with error bars
    plt.errorbar(x, y, yerr=yerr, fmt='o', ecolor='r', capsize=5, capthick=2, markeredgecolor='blue',
                 markerfacecolor='green')
    plt.xlabel('X-axis label')
    plt.ylabel('Y-axis label')
    plt.title('Plot with Error Bars')
    plt.show()

def main_indirect_max():
    unGraph, g_stat_1, g_stat_2, m_indx1_dict, m_indx2_dict, author_hot_streaks = calc_mIndx()
    num1 = 0
    num2 = 0
    induced_vertex = []
    for v in unGraph.nodes():
        if g_stat_1[v] == 0 and g_stat_2[v] ==1 and (m_indx1_dict[v]>=2 or m_indx2_dict[v]>=2):
            num1 = num1 + 1
            induced_vertex.append(v)
        if g_stat_1[v] == 0 and g_stat_2[v] ==0 and (m_indx1_dict[v]>=2 or m_indx2_dict[v]>=2):
            num2 = num2 + 1
    print("2")
    induced_m_idx = {}
    induced_hot_streaks = {}

    Hot_streak = {}
    div_num = 5
    for i in np.arange(div_num):
        Hot_streak[i] = []

    for vertex in induced_vertex:
        if m_indx2_dict[vertex] != 0:
            induced_m_idx[vertex] = m_indx2_dict[vertex]
            induced_hot_streaks[vertex] = author_hot_streaks[vertex][1][1]

    sorted_m_idx = sorted(induced_m_idx.items(), key=lambda item: item[1])
    #sorted_hot_streaks = sorted(induced_hot_streaks.items(), key=lambda item: item[1])
    #print(sorted_hot_streaks)

    top_20_per_count = int(len(induced_m_idx) * 0.2)
    top_20_percent_keys = [item[0] for item in sorted_m_idx[:top_20_per_count]]
    top_20_percent_values = [item[1] for item in sorted_m_idx[:top_20_per_count]]
    hot_streaks = []
    for key in top_20_percent_keys:
        hot_stk = induced_hot_streaks[key]
        hot_streaks.append(hot_stk)
    x = str(min(top_20_percent_values))+"-"+str(max(top_20_percent_values))
    y = statistics.mean(hot_streaks)
    yerr = statistics.stdev(hot_streaks)
    plt.errorbar(x, y, yerr=yerr, fmt='o', ecolor='r', capsize=5, capthick=2, markeredgecolor='blue',
                 markerfacecolor='green')

    for ii in range(1,5):
        top_x_per_count = int(len(induced_m_idx) * 0.2 * ii)
        top_x1_per_count = int(len(induced_m_idx) * 0.2 * (ii+1))
        top_x1_percent_keys = [item[0] for item in sorted_m_idx[top_x_per_count+1:top_x1_per_count]]
        top_x1_percent_values = [item[1] for item in sorted_m_idx[top_x_per_count+1:top_x1_per_count]]
        hot_streaks = []
        for key in top_x1_percent_keys:
            hot_stk = induced_hot_streaks[key]
            hot_streaks.append(hot_stk)
        x = str(min(top_x1_percent_values)) + "-" + str(max(top_x1_percent_values))
        y = statistics.mean(hot_streaks)
        yerr = statistics.stdev(hot_streaks)
        plt.errorbar(x, y, yerr=yerr, fmt='o', ecolor='r', capsize=5, capthick=2, markeredgecolor='blue',
                 markerfacecolor='green')

    plt.xlabel('Induced index m(i)')
    plt.ylabel('Hot streaks')
    plt.title('Hot streaks of computer scientists from 2002-2011')
    plt.show()

def main_indirect_max2():
    unGraph, g_stat_1, g_stat_2, m_indx1_dict, m_indx2_dict, author_hot_streaks = calc_mIndx()
    num1 = 0
    num2 = 0
    induced_vertex = []

    x_list = [i for i in range(1,21)]
    y_list = []
    y_err_list = []

    for v in unGraph.nodes():
        if g_stat_1[v] == 0 and g_stat_2[v] ==1 and (m_indx1_dict[v]>=2 or m_indx2_dict[v]>=2):
            num1 = num1 + 1
            induced_vertex.append(v)
        if g_stat_1[v] == 0 and g_stat_2[v] ==0 and (m_indx1_dict[v]>=2 or m_indx2_dict[v]>=2):
            num2 = num2 + 1
    print("2")

    induced_m_idx = {}
    induced_hot_streaks = {}

    for vertex in induced_vertex:
        if m_indx2_dict[vertex] != 0:
            induced_m_idx[vertex] = m_indx2_dict[vertex]
            induced_hot_streaks[vertex] = author_hot_streaks[vertex][1][1]

    sorted_m_idx = sorted(induced_m_idx.items(), key=lambda item: item[1])
    top_20_per_count = int(len(induced_m_idx) * 0.05)
    top_20_percent_keys = [item[0] for item in sorted_m_idx[:top_20_per_count]]
    top_20_percent_values = [item[1] for item in sorted_m_idx[:top_20_per_count]]
    hot_streaks = []
    for key in top_20_percent_keys:
        hot_stk = induced_hot_streaks[key]
        if hot_stk !=0:
            hot_streaks.append(hot_stk)
    x = str(min(top_20_percent_values))+"-"+str(max(top_20_percent_values))

    y,yerr = test_bootstrap_statistics(hot_streaks)
    #y = statistics.mean(hot_streaks)
    y_list.append(y)
    #yerr = statistics.stdev(hot_streaks)
    y_err_list.append(yerr)
    plt.errorbar(x, y, yerr=yerr, fmt='o', ecolor='r', capsize=5, capthick=2, markeredgecolor='blue',
                 markerfacecolor='green')

    for ii in range(1,20):
        top_x_per_count = int(len(induced_m_idx) * 0.05 * ii)
        top_x1_per_count = int(len(induced_m_idx) * 0.05 * (ii+1))
        top_x1_percent_keys = [item[0] for item in sorted_m_idx[top_x_per_count+1:top_x1_per_count]]
        top_x1_percent_values = [item[1] for item in sorted_m_idx[top_x_per_count+1:top_x1_per_count]]
        hot_streaks = []
        for key in top_x1_percent_keys:
            hot_stk = induced_hot_streaks[key]
            hot_streaks.append(hot_stk)
        x = str(min(top_x1_percent_values)) + "-" + str(max(top_x1_percent_values))
        #y = statistics.mean(hot_streaks)
        y, yerr = test_bootstrap_statistics(hot_streaks)
        y_list.append(y)
        #yerr = statistics.stdev(hot_streaks)
        y_err_list.append(yerr)
        plt.errorbar(x, y, yerr=yerr, fmt='o', ecolor='r', capsize=5, capthick=2, markeredgecolor='blue',
                 markerfacecolor='green')

    print(y_list)
    print(y_err_list)

    plt.xlabel('Induced Index m(i)')
    plt.ylabel('Max Hot Streaks')
    plt.title('Hot streaks of computer scientists from 2002-2011')
    plt.show()

def main_direct_max():
    unGraph, g_stat_1, g_stat_2, m_deg_dict1, m_deg_dict2, author_hot_streaks = calc_degree()
    print("1")
    num1 = 0
    num2 = 0
    induced_vertex = []
    for v in unGraph.nodes():
        if g_stat_1[v] == 0 and g_stat_2[v] ==1 and (m_deg_dict1[v]>=1 and m_deg_dict2[v]>=1):
            num1 = num1 + 1
            induced_vertex.append(v)
        if g_stat_1[v] == 0 and g_stat_2[v] ==0 and (m_deg_dict1[v]>=1 and m_deg_dict2[v]>=1):
            num2 = num2 + 1
    print("2")
    induced_m_idx = {}
    induced_hot_streaks = {}

    x_list = [i for i in range(1, 11)]
    y_list = []
    y_err_list = []

    for vertex in induced_vertex:
        if m_deg_dict2[vertex] != 0:
            induced_m_idx[vertex] = m_deg_dict2[vertex]
            induced_hot_streaks[vertex] = author_hot_streaks[vertex][1][1]

    sorted_m_idx = sorted(induced_m_idx.items(), key=lambda item: item[1])
    top_10_per_count = int(len(induced_m_idx) * 0.05)
    top_10_percent_keys = [item[0] for item in sorted_m_idx[:top_10_per_count]]
    top_10_percent_values = [item[1] for item in sorted_m_idx[:top_10_per_count]]
    hot_streaks = []
    for key in top_10_percent_keys:
        hot_stk = induced_hot_streaks[key]
        if hot_stk != 0:
            hot_streaks.append(hot_stk)
    x = str(min(top_10_percent_values)) + "-" + str(max(top_10_percent_values))
    y,yerr = test_bootstrap_statistics(hot_streaks)
    #y = statistics.mean(hot_streaks)
    #yerr = statistics.stdev(hot_streaks)
    y_list.append(y)
    y_err_list.append(yerr)
    plt.errorbar(x, y, yerr=yerr, fmt='o', ecolor='r', capsize=5, capthick=2, markeredgecolor='blue',
                 markerfacecolor='green')

    for ii in range(1, 10):
        top_x_per_count = int(len(induced_m_idx) * 0.1 * ii)
        top_x1_per_count = int(len(induced_m_idx) * 0.1 * (ii + 1))
        top_x1_percent_keys = [item[0] for item in sorted_m_idx[top_x_per_count + 1:top_x1_per_count]]
        top_x1_percent_values = [item[1] for item in sorted_m_idx[top_x_per_count + 1:top_x1_per_count]]
        hot_streaks = []
        for key in top_x1_percent_keys:
            hot_stk = induced_hot_streaks[key]
            if hot_stk != 0:
                hot_streaks.append(hot_stk)
        x = str(min(top_x1_percent_values)) + "-" + str(max(top_x1_percent_values))
        y, yerr = test_bootstrap_statistics(hot_streaks)
        #y = statistics.mean(hot_streaks)
        #yerr = statistics.stdev(hot_streaks)
        y_list.append(y)
        y_err_list.append(yerr)
        plt.errorbar(x, y, yerr=yerr, fmt='o', ecolor='r', capsize=5, capthick=2, markeredgecolor='blue',
                     markerfacecolor='green')

    print(y_list)
    print(y_err_list)

    plt.xlabel('Degree d(i)')
    plt.ylabel('Hot streaks')
    plt.title('Hot streaks of computer scientists from 2002-2011')
    plt.show()

def main_direct_indirect_max2():
    unGraph, g_stat_1, g_stat_2, deg_indx1_dict, deg_indx2_dict, m_indx1_dict, m_indx2_dict, author_hot_streak = calc_degree_induce_m()
    print("1")
    num1 = 0
    num2 = 0
    induced_vertex = []
    for v in unGraph.nodes():
        if g_stat_1[v] == 0 and g_stat_2[v] ==1 and (deg_indx1_dict[v]>=1 and deg_indx2_dict[v]>=1 ):
            num1 = num1 + 1
            induced_vertex.append(v)
        if g_stat_1[v] == 0 and g_stat_2[v] ==0 and (deg_indx1_dict[v]>=1 and deg_indx2_dict[v]>=1):
            num2 = num2 + 1
    print("2")
    induced_m_idx = {}
    induced_hot_streaks = {}

    for vertex in induced_vertex:
        if deg_indx2_dict[vertex] != 0:
            induced_m_idx[vertex] = deg_indx2_dict[vertex]
            induced_hot_streaks[vertex] = author_hot_streak[vertex][1][1]

    sorted_m_idx = sorted(induced_m_idx.items(), key=lambda item: item[1])
    top_10_per_count = int(len(induced_m_idx) * 0.1)
    top_10_percent_keys = [item[0] for item in sorted_m_idx[:top_10_per_count]]
    top_10_percent_values = [item[1] for item in sorted_m_idx[:top_10_per_count]]
    hot_streaks = []
    y_list = []
    y_err_list = []
    for key in top_10_percent_keys:
        hot_stk = induced_hot_streaks[key]
        if hot_stk != 0:
            hot_streaks.append(hot_stk)
    x = str(min(top_10_percent_values)) + "-" + str(max(top_10_percent_values))
    #y = statistics.mean(hot_streaks)
    #yerr = statistics.stdev(hot_streaks)
    y,yerr = test_bootstrap_statistics(hot_streaks)
    y_list.append(y)
    y_err_list.append(yerr)
    plt.errorbar(x, y, yerr=yerr, fmt='o', ecolor='r', capsize=5, capthick=2, markeredgecolor='blue',
                 markerfacecolor='green')

    for ii in range(1, 10):
        top_x_per_count = int(len(induced_m_idx) * 0.1 * ii)
        top_x1_per_count = int(len(induced_m_idx) * 0.1 * (ii + 1))
        top_x1_percent_keys = [item[0] for item in sorted_m_idx[top_x_per_count + 1:top_x1_per_count]]
        top_x1_percent_values = [item[1] for item in sorted_m_idx[top_x_per_count + 1:top_x1_per_count]]
        hot_streaks = []
        for key in top_x1_percent_keys:
            hot_stk = induced_hot_streaks[key]
            hot_streaks.append(hot_stk)
        x = str(min(top_x1_percent_values)) + "-" + str(max(top_x1_percent_values))
        #y = statistics.mean(hot_streaks)
        #yerr = statistics.stdev(hot_streaks)
        y, yerr = test_bootstrap_statistics(hot_streaks)
        y_list.append(y)
        y_err_list.append(yerr)
        plt.errorbar(x, y, yerr=yerr, fmt='o', ecolor='r', capsize=5, capthick=2, markeredgecolor='blue',
                     markerfacecolor='green')
    print(y_list)
    print(y_err_list)
    plt.xlabel('Degree d(i)')
    plt.ylabel('Mean Hot streaks')
    plt.title('Hot streaks of computer scientists from 2002-2011')
    plt.show()

def main_indirect_mean():
    unGraph, g_stat_1, g_stat_2, deg_indx1_dict, deg_indx2_dict, m_indx1_dict, m_indx2_dict, author_hot_streak = calc_degree_induce_m()
    print("1")
    num1 = 0
    num2 = 0
    induced_vertex = []
    for v in unGraph.nodes():
        if g_stat_1[v] == 0 and g_stat_2[v] ==1 and (m_indx1_dict[v]>=1 and m_indx2_dict[v]>=1 ):
            num1 = num1 + 1
            induced_vertex.append(v)
        if g_stat_1[v] == 0 and g_stat_2[v] ==0 and (m_indx1_dict[v]>=1 and m_indx2_dict[v]>=1):
            num2 = num2 + 1
    print("2")
    induced_m_idx = {}
    induced_hot_streaks = {}

    for vertex in induced_vertex:
        if m_indx2_dict[vertex] != 0:
            induced_m_idx[vertex] = m_indx2_dict[vertex]
            induced_hot_streaks[vertex] = author_hot_streak[vertex][1][1]

    sorted_m_idx = sorted(induced_m_idx.items(), key=lambda item: item[1])
    top_10_per_count = int(len(induced_m_idx) * 0.05)
    top_10_percent_keys = [item[0] for item in sorted_m_idx[:top_10_per_count]]
    top_10_percent_values = [item[1] for item in sorted_m_idx[:top_10_per_count]]
    hot_streaks = []
    y_list = []
    y_err_list = []
    for key in top_10_percent_keys:
        hot_stk = induced_hot_streaks[key]
        if hot_stk != 0:
            hot_streaks.append(hot_stk)
    x = str(min(top_10_percent_values)) + "-" + str(max(top_10_percent_values))
    #y = statistics.mean(hot_streaks)
    #yerr = statistics.stdev(hot_streaks)
    y,yerr = test_bootstrap_statistics(hot_streaks)
    y_list.append(y)
    y_err_list.append(yerr)
    plt.errorbar(x, y, yerr=yerr, fmt='o', ecolor='r', capsize=5, capthick=2, markeredgecolor='blue',
                 markerfacecolor='green')

    for ii in range(1, 20):
        top_x_per_count = int(len(induced_m_idx) * 0.05 * ii)
        top_x1_per_count = int(len(induced_m_idx) * 0.05 * (ii + 1))
        top_x1_percent_keys = [item[0] for item in sorted_m_idx[top_x_per_count + 1:top_x1_per_count]]
        top_x1_percent_values = [item[1] for item in sorted_m_idx[top_x_per_count + 1:top_x1_per_count]]
        hot_streaks = []
        for key in top_x1_percent_keys:
            hot_stk = induced_hot_streaks[key]
            hot_streaks.append(hot_stk)
        x = str(min(top_x1_percent_values)) + "-" + str(max(top_x1_percent_values))
        #y = statistics.mean(hot_streaks)
        #yerr = statistics.stdev(hot_streaks)
        y, yerr = test_bootstrap_statistics(hot_streaks)
        y_list.append(y)
        y_err_list.append(yerr)
        plt.errorbar(x, y, yerr=yerr, fmt='o', ecolor='r', capsize=5, capthick=2, markeredgecolor='blue',
                     markerfacecolor='green')
    print(y_list)
    print(y_err_list)
    plt.xlabel('Induced Index m(i)')
    plt.ylabel('Mean Hot streaks')
    plt.title('Hot streaks of computer scientists from 2002-2011')
    plt.show()

def main_indirect_max3_citations():
        unGraph, g_stat_1, g_stat_2, m_indx1_dict, m_indx2_dict, author_citations = calc_mIndx_citations()
        num1 = 0
        num2 = 0
        induced_vertex = []
        for v in unGraph.nodes():
            if g_stat_1[v] == 0 and g_stat_2[v] == 1 and (m_indx1_dict[v] >= 1 or m_indx2_dict[v] >= 1):
                num1 = num1 + 1
                induced_vertex.append(v)
            if g_stat_1[v] == 0 and g_stat_2[v] == 0 and (m_indx1_dict[v] >= 1 or m_indx2_dict[v] >= 1):
                num2 = num2 + 1

        induced_m_idx = {}
        induced_citations = {}

        for vertex in induced_vertex:
            if m_indx2_dict[vertex] != 0:
                induced_m_idx[vertex] = m_indx2_dict[vertex]
                induced_citations[vertex] = author_citations[vertex][1]

        sorted_m_idx = sorted(induced_m_idx.items(), key=lambda item: item[1])

        top_count = int(len(induced_m_idx) * 0.95)
        top_keys = [item[0] for item in sorted_m_idx[:top_count] if item[1] <= 2]
        top_values = [item[1] for item in sorted_m_idx[:top_count] if item[1] <= 2]

        cit_ations = []
        for key in top_keys:
            citation = induced_citations[key]
            if citation != 0:
                cit_ations.append(citation)
        x = str(min(top_values)) + "-" + str(max(top_values))
        y = statistics.mean(cit_ations)
        yerr = statistics.stdev(cit_ations)
        plt.errorbar(x, y, yerr=yerr, fmt='o', ecolor='r', capsize=5, capthick=2, markeredgecolor='blue',
                     markerfacecolor='green')
        #plt.plot(x, y)
        for ii in range(2, 20):
            top_x1_percent_keys = [item[0] for item in sorted_m_idx[:top_count] if item[1] >= ii and item[1] <= 1 + ii]
            top_x1_percent_values = [item[1] for item in sorted_m_idx[:top_count] if
                                     item[1] >= ii and item[1] <= 1 + ii]
            citation_list = []
            for key in top_x1_percent_keys:
                cit_ation = induced_citations[key]
                citation_list.append(cit_ation)
            x = str((min(top_x1_percent_values) if len(top_x1_percent_values) > 0 else 0)) + "-" + str(
                max(top_x1_percent_values) if len(top_x1_percent_values) > 0 else 0)
            y = statistics.mean(citation_list)
            yerr = statistics.stdev(citation_list)
            plt.errorbar(x, y, yerr=yerr, fmt='o', ecolor='r', capsize=5, capthick=2, markeredgecolor='blue',
                         markerfacecolor='green')
            #plt.plot(x,y)
        plt.xlabel('Induced index m(i)')
        plt.ylabel('Mean citations')
        plt.title('Hot streaks of computer scientists from 2002-2011')
        plt.show()

def main_indirect_max4():
    unGraph, g_stat_1, g_stat_2, m_indx1_dict, m_indx2_dict, author_hot_streaks = calc_mIndx()
    print("1")
    num1 = 0
    num2 = 0
    induced_vertex = []
    for v in unGraph.nodes():
        if g_stat_1[v] == 0 and g_stat_2[v] == 1 and (m_indx1_dict[v] >= 1 or m_indx2_dict[v] >= 1):
            num1 = num1 + 1
            induced_vertex.append(v)
        if g_stat_1[v] == 0 and g_stat_2[v] == 0 and (m_indx1_dict[v] >= 1 or m_indx2_dict[v] >= 1):
            num2 = num2 + 1
    print("2")
    induced_m_idx = {}
    induced_hot_streaks = {}

    for vertex in induced_vertex:
        if m_indx2_dict[vertex] != 0:
            induced_m_idx[vertex] = m_indx2_dict[vertex]
            induced_hot_streaks[vertex] = author_hot_streaks[vertex][1][1]

    sorted_m_idx = sorted(induced_m_idx.items(), key=lambda item: item[1])

    top_count = int(len(induced_m_idx) * 0.95)
    top_keys = [item[0] for item in sorted_m_idx[:top_count] if item[1] <= 2]
    top_values = [item[1] for item in sorted_m_idx[:top_count] if item[1] <= 2]

    hot_streaks = []
    for key in top_keys:
        hot_stk = induced_hot_streaks[key]
        if hot_stk != 0:
            hot_streaks.append(hot_stk)
    x = str(min(top_values)) + "-" + str(max(top_values))
    y = statistics.mean(hot_streaks)
    yerr = statistics.stdev(hot_streaks)
    plt.errorbar(x, y, yerr=yerr, fmt='o', ecolor='r', capsize=5, capthick=2, markeredgecolor='blue',
                 markerfacecolor='green')

    for ii in range(2, 20):
        top_x1_percent_keys = [item[0] for item in sorted_m_idx[:top_count] if item[1] >= ii and item[1] <= 1 + ii]
        top_x1_percent_values = [item[1] for item in sorted_m_idx[:top_count] if item[1] >= ii and item[1] <= 1 + ii]
        hot_streaks = []
        for key in top_x1_percent_keys:
            hot_stk = induced_hot_streaks[key]
            hot_streaks.append(hot_stk)
        x = str((min(top_x1_percent_values) if len(top_x1_percent_values) > 0 else 0)) + "-" + str(
            max(top_x1_percent_values) if len(top_x1_percent_values) > 0 else 0)
        y = statistics.mean(hot_streaks)
        yerr = statistics.stdev(hot_streaks)
        plt.errorbar(x, y, yerr=yerr, fmt='o', ecolor='r', capsize=5, capthick=2, markeredgecolor='blue',
                     markerfacecolor='green')
    plt.xlabel('Induced index m(i)')
    plt.ylabel('Max Hot streaks')
    plt.title('Hot streaks of computer scientists from 2002-2011')
    plt.show()

def main_direct_max3():
    unGraph, g_stat_1, g_stat_2, m_deg_dict1, m_deg_dict2, author_hot_streaks = calc_degree()
    print("1")
    num1 = 0
    num2 = 0
    induced_vertex = []
    for v in unGraph.nodes():
        if g_stat_1[v] == 0 and g_stat_2[v] == 1 and (m_deg_dict1[v] >= 1 or m_deg_dict2[v] >= 1):
            num1 = num1 + 1
            induced_vertex.append(v)
        if g_stat_1[v] == 0 and g_stat_2[v] == 0 and (m_deg_dict1[v] >= 1 or m_deg_dict2[v] >= 1):
            num2 = num2 + 1
    print("2")
    induced_m_idx = {}
    induced_hot_streaks = {}

    for vertex in induced_vertex:
        if m_deg_dict2[vertex] != 0:
            induced_m_idx[vertex] = m_deg_dict2[vertex]
            induced_hot_streaks[vertex] = author_hot_streaks[vertex][1][0]

    sorted_m_idx = sorted(induced_m_idx.items(), key=lambda item: item[1])

    top_count = int(len(induced_m_idx) * 0.95)
    top_keys = [item[0] for item in sorted_m_idx[:top_count] if item[1] <= 2]
    top_values = [item[1] for item in sorted_m_idx[:top_count] if item[1] <= 2]

    hot_streaks = []
    for key in top_keys:
        hot_stk = induced_hot_streaks[key]
        if hot_stk != 0:
            hot_streaks.append(hot_stk)
    x = str(min(top_values)) + "-" + str(max(top_values))
    y = statistics.mean(hot_streaks)
    yerr = statistics.stdev(hot_streaks)
    plt.errorbar(x, y, yerr=yerr, fmt='o', ecolor='r', capsize=5, capthick=2, markeredgecolor='blue',
                 markerfacecolor='green')

    for ii in range(2, 20):
        top_x1_percent_keys = [item[0] for item in sorted_m_idx[:top_count] if item[1] >= ii and item[1] <= 1 + ii]
        top_x1_percent_values = [item[1] for item in sorted_m_idx[:top_count] if item[1] >= ii and item[1] <= 1 + ii]
        hot_streaks = []
        for key in top_x1_percent_keys:
            hot_stk = induced_hot_streaks[key]
            hot_streaks.append(hot_stk)
        x = str((min(top_x1_percent_values) if len(top_x1_percent_values) > 0 else 0)) + "-" + str(
            max(top_x1_percent_values) if len(top_x1_percent_values) > 0 else 0)
        y = statistics.mean(hot_streaks)
        yerr = statistics.stdev(hot_streaks)
        plt.errorbar(x, y, yerr=yerr, fmt='o', ecolor='r', capsize=5, capthick=2, markeredgecolor='blue',
                     markerfacecolor='green')
    plt.xlabel('Degree Index d(i)')
    plt.ylabel('Mean Hot streaks')
    plt.title('Hot streaks of computer scientists from 2002-2011')
    plt.show()

def main_direct():
    unGraph, g_stat_1, g_stat_2, m_deg_dict1, m_deg_dict2, author_hot_streaks = calc_degree()
    print("1")
    num1 = 0
    num2 = 0
    induced_vertex = []
    for v in unGraph.nodes():
        if g_stat_1[v] == 0 and g_stat_2[v] ==1 and (m_deg_dict1[v]>=2 or m_deg_dict2[v]>=2):
            num1 = num1 + 1
            induced_vertex.append(v)
        if g_stat_1[v] == 0 and g_stat_2[v] ==0 and (m_deg_dict1[v]>=2 or m_deg_dict2[v]>=2):
            num2 = num2 + 1
    print("2")
    induced_vertex_set = set()
    induced_m_idx = set()
    Hot_streak = {}
    div_num = 20
    for i in np.arange(div_num):
        Hot_streak[i] = []
    for vertex in induced_vertex:
        induced_vertex_set.add(vertex)
        induced_m_idx.add(m_deg_dict2[vertex])
    print("3")
    ##min_value, max_value = stat_data(list(m_deg_dict2.values()))
    min_value = min(induced_m_idx)
    max_value = max(induced_m_idx)
    print(min_value)
    print(max_value)
    step_ = (max_value - min_value) // div_num
    bound = []
    for i in np.arange(min_value,max_value+1,step=step_):
        if i > min_value:
            bound.append(i)
    print("4")
    for vertex in induced_vertex:
        m_v = m_deg_dict2[vertex]
        hot_streak = int(author_hot_streaks[vertex][1][1])
        if m_v <= bound[0]:
            Hot_streak[0].append(hot_streak)
        elif m_v <= bound[1]:
            Hot_streak[1].append(hot_streak)
        elif m_v <= bound[2]:
            Hot_streak[2].append(hot_streak)
        elif m_v <= bound[3]:
            Hot_streak[3].append(hot_streak)
        elif m_v <= bound[4]:
            Hot_streak[4].append(hot_streak)
        elif m_v <= bound[5]:
            Hot_streak[5].append(hot_streak)
        elif m_v <= bound[6]:
            Hot_streak[6].append(hot_streak)
        elif m_v <= bound[7]:
            Hot_streak[7].append(hot_streak)
        elif m_v <= bound[8]:
            Hot_streak[8].append(hot_streak)
        elif m_v <= bound[9]:
            Hot_streak[9].append(hot_streak)
    print("5")
    x_label = "Degree d"
    y_label = "Distributions of hot streaks"
    display_10(Hot_streak,x_label,y_label)

def display_degree_distributions(G):
    degree_counts = nx.degree_histogram(G)
    degrees = range(len(degree_counts))
    # 绘制度分布
    plt.figure(figsize=(10, 6))
    plt.bar(degrees, degree_counts, width=0.8, color='b')
    plt.title('网络的度分布')
    plt.xlabel('度数')
    plt.ylabel('节点数量')
    plt.yscale('log')  # 设置 y 轴为对数尺度（如果度分布跨度较大时）
    plt.show()

def stat_data(data):
    # 计算统计数据
    Q1 = np.percentile(data, 25)
    median = np.percentile(data, 50)
    Q3 = np.percentile(data, 75)
    IQR = Q3 - Q1
    lower_whisker = Q1 - 1.5 * IQR
    upper_whisker = Q3 + 1.5 * IQR

    # 计算最小值和最大值
    min_val = np.min([ i for i in data if i >= lower_whisker])
    max_val = np.max([ i for i in data if i <= upper_whisker])

    # 异常值
    outliers = [ i for i in data if i < lower_whisker or i > upper_whisker]

    # 打印统计数据
    print(f"Q1 (25th percentile): {Q1}")
    print(f"Median: {median}")
    print(f"Q3 (75th percentile): {Q3}")
    print(f"IQR: {IQR}")
    print(f"Lower Whisker: {lower_whisker}")
    print(f"Upper Whisker: {upper_whisker}")
    print(f"Min (within whisker): {min_val}")
    print(f"Max (within whisker): {max_val}")
    print(f"Outliers: {outliers}")
    return min_val,max_val

def error_bar():
    # Sample data
    x = np.linspace(0, 10, 10)
    y = np.sin(x)
    yerr = 0.2  # Constant error

    # Create the plot
    plt.errorbar(x, y, yerr=yerr, fmt='-o', ecolor='red', capsize=5, capthick=2, marker='s', markersize=5,
                 label='Data with Error')

    # Customize the plot
    plt.title('Plot with Error Bars')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.legend()
    plt.grid(True)

    # Display the plot
    plt.show()

def test_sort_data():
    data = {
        'apple': 10,
        'banana': 2,
        'cherry': 5,
        'date': 7
    }

    # 按照值对字典进行排序，并返回排序后的键
    sorted_keys = sorted(data, key=lambda x: data[x])
    # 输出排序后的键
    print(sorted_keys)

def test_sort_data_2():
    # 示例字典
    data = {
        'a': 10,
        'b': 20,
        'c': 5,
        'd': 30,
        'e': 15,
        'f': 25,
        'g': 8,
        'h': 12,
        'i': 35,
        'j': 22,
    }

    # 按照值排序
    sorted_data = sorted(data.items(), key=lambda item: item[1])

    # 计算前10%的数据数量
    top_10_percent_count = int(len(data) * 0.3)

    # 获取前10%的键
    top_10_percent_keys = [item[0] for item in sorted_data[:top_10_percent_count]]
    top_10_percent_values = [item[1] for item in sorted_data[:top_10_percent_count]]

    print(top_10_percent_keys)
    print(top_10_percent_values)

def test_two_way_ANOVA():
    # Sample data
    # Sample data
    data = {
        'Condition1': [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4],
        'Condition2': [10, 15, 20, 10, 15, 20, 10, 15, 20, 10, 15, 20],
        'Score': [23, 20, 27, 30, 35, 34, 25, 28, 23, 22, 19, 18]
    }

    # Create DataFrame
    df = pd.DataFrame(data)

    # Fit the model
    model = ols('Score ~ Condition1 + Condition2 + Condition1:Condition2', data=df).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)

    print(anova_table)


def main_anova_1():
    unGraph, g_stat_1, g_stat_2, deg_indx1_dict, deg_indx2_dict, m_indx1_dict, m_indx2_dict, author_hot_streaks = calc_degree_induce_m()
    print("1")
    induced_vertex = []
    for v in unGraph.nodes():
        if (g_stat_1[v] == 0 and g_stat_2[v] == 1) and (deg_indx1_dict[v] >= 1 and deg_indx2_dict[v] >= 1) and (m_indx1_dict[v]>=1 and m_indx2_dict[v]>=1):
            induced_vertex.append(v)
    print("2")
    d_idx=[]
    m_idx=[]
    y = []
    num=0
    for vertex in induced_vertex:
        d_idx.append(deg_indx2_dict[vertex])
        m_idx.append(m_indx1_dict[vertex])
        y.append(author_hot_streaks[vertex][1][0])
        num = num + 1
    print(num)

    data = {
        'Condition1': d_idx,
        'Condition2': m_idx,
        'Score': y
    }
    # Create DataFrame
    df = pd.DataFrame(data)

    # Fit the model
    model = ols('Score ~ Condition1 + Condition2 + Condition1:Condition2', data=df).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)
    print(anova_table)

def main_anova_2():
    unGraph, g_stat_1, g_stat_2, deg_indx1_dict, deg_indx2_dict, m_indx1_dict, m_indx2_dict, author_hot_streaks = calc_degree_induce_m()
    print("1")
    induced_vertex = []
    for v in unGraph.nodes():
        if (g_stat_1[v] == 0 and g_stat_2[v] == 1) and (deg_indx1_dict[v] >= 1 and deg_indx2_dict[v] >= 1) and (m_indx1_dict[v]>=1 and m_indx2_dict[v]>=1):
            induced_vertex.append(v)
    print("2")
    d_idx=[]
    m_idx=[]
    y = []
    num=0
    for vertex in induced_vertex:
        if author_hot_streaks[vertex][1][0] !=0:
            d_idx.append(deg_indx2_dict[vertex])
            m_idx.append(m_indx1_dict[vertex])
            y.append(author_hot_streaks[vertex][1][0])
            num = num + 1
    print(num)

    data = {
        'Condition1': d_idx,
        'Condition2': m_idx,
        'Score': y
    }
    # Create DataFrame
    df = pd.DataFrame(data)

    # Fit the model
    model = ols('Score ~ Condition1 + Condition2 + Condition1:Condition2', data=df).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)
    #anova_table['mean_sq'] = anova_table['sum_sq'] / anova_table['df']
    anova_table['eta_sq'] = anova_table['sum_sq'] / (anova_table['sum_sq'] + anova_table.loc['Residual', 'sum_sq'])
    print(anova_table)

# 定义一个函数来执行Bootstrap抽样并计算均值和方差
def bootstrap_statistics(data, n_iterations):
    means = []
    variances = []
    n_size = len(data)

    # 进行n_iterations次采样
    for i in range(n_iterations):
        # 从原始数据中有放回地随机采样
        sample = np.random.choice(data, size=n_size, replace=True)
        # 计算每个样本的均值和方差
        sample_mean = np.mean(sample)
        sample_variance = np.var(sample, ddof=1)  # ddof=1是样本方差
        means.append(sample_mean)
        variances.append(sample_variance)

    # 返回均值和方差的结果
    return means, variances

def test_bootstrap_statistics(data):
    # 示例数据
    #data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    # 执行Bootstrap抽样并计算均值和方差
    n_iterations = 1000  # 设置抽样次数
    bootstrap_means, bootstrap_variances = bootstrap_statistics(data, n_iterations)
    # 计算Bootstrap均值和方差的平均值
    mean_of_means = np.mean(bootstrap_means)
    mean_of_variances = np.mean(bootstrap_variances)

    # 输出结果
    print(f"Bootstrap均值的平均值: {mean_of_means}")
    print(f"Bootstrap方差的标准差: {mean_of_variances}")
    return (mean_of_means,mean_of_variances)

def gen_graph_gexf():
    # 设置随机种子以保证结果可复现
    # 示例邻接矩阵，表示无向图
    G = nx.Graph()
    G.add_node(1,color='red')
    G.add_node(2, color='blue')
    G.add_node(3, color='blue')
    for i in range(4,14):
        G.add_node(i,color='yellow')

    G.add_edge(1,2)
    G.add_edge(1,3)
    G.add_edge(2,4)
    G.add_edge(2,5)
    G.add_edge(2,6)
    G.add_edge(2,7)
    G.add_edge(2,8)
    G.add_edge(2,9)
    G.add_edge(3,10)
    G.add_edge(3,11)
    G.add_edge(3,12)
    G.add_edge(3,13)
    nx.write_gexf(G, "graph_1.gexf")
    print("Graph saved as graph_output.gexf")

    G = nx.Graph()
    G.add_node(14,color='red')
    G.add_node(15, color='blue')
    G.add_node(16, color='blue')
    for i in range(17, 28):
        G.add_node(i,color='yellow')
    G.add_edge(14, 15)
    G.add_edge(14, 16)
    G.add_edge(15, 17)
    G.add_edge(15, 18)
    G.add_edge(15, 19)
    G.add_edge(15, 20)
    G.add_edge(15, 21)
    G.add_edge(15, 22)
    G.add_edge(15, 23)
    G.add_edge(16, 24)
    G.add_edge(16, 25)
    G.add_edge(16, 26)
    G.add_edge(16, 27)
    # 将图写入到GEXF文件
    nx.write_gexf(G, "graph_2.gexf")
    print("Graph saved as graph_output.gexf")

    G = nx.Graph()
    G.add_node(29, color='red')
    G.add_node(30, color='blue')
    G.add_node(31, color='blue')
    G.add_node(32, color='blue')
    for i in range(33, 43):
        G.add_node(i,color='yellow')
    G.add_edge(29, 30)
    G.add_edge(29, 31)
    G.add_edge(29, 32)
    G.add_edge(30, 33)
    G.add_edge(30, 34)
    G.add_edge(30, 35)
    G.add_edge(30, 36)
    G.add_edge(30, 37)
    G.add_edge(30, 38)
    G.add_edge(30, 39)
    G.add_edge(31, 40)
    G.add_edge(31, 41)
    G.add_edge(32, 42)
    nx.write_gexf(G, "graph_3.gexf")
    print("Graph saved as graph_output.gexf")

    G = nx.Graph()
    G.add_node(44, color='red')
    G.add_node(45, color='blue')
    G.add_node(46, color='blue')
    for i in range(47, 62):
        G.add_node(i,color='yellow')
    G.add_edge(44, 45)
    G.add_edge(44, 46)
    G.add_edge(45, 47)
    G.add_edge(45, 48)
    G.add_edge(45, 49)
    G.add_edge(45, 50)
    G.add_edge(45, 51)
    G.add_edge(45, 52)
    G.add_edge(45, 53)
    G.add_edge(45, 54)
    G.add_edge(45, 55)
    G.add_edge(45, 56)
    G.add_edge(45, 57)
    G.add_edge(45, 58)
    G.add_edge(46, 59)
    G.add_edge(46, 60)
    G.add_edge(46, 61)
    nx.write_gexf(G, "graph_4.gexf")
    print("Graph saved as graph_output.gexf")

if __name__ == "__main__":
    gen_graph_gexf()
    #main_direct_max()
    #main_indirect_max2()
    #test_bootstrap_statistics()
    #main_indirect_max3_citations()
    #test_sort_data_2()
    #main_anova_2()
    #test_two_way_ANOVA()
    #main_indirect_max3_citations()

