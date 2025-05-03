# growth-induced percolation
This project contains source codes and data sets from the paper "Growth-induced Percolation on Complex Networks" by H.L. SUN, S.H. Chen, J.R. Xie and Y.Q.Hu.

#1. Source Codes:
The source codes implemented in this project are not strictly along with equations, but with optimization from ER and SF networks.

Generally speaking, there are six parts here, including directed ER and undirected ER networks, directed SF and undirected SF networks,  experimental study from hot streaks as well as common functions.
Each part contains theoretical study and simulation on extremely large networks.
(1) directed ER networks: direc_ER_numerical.py, direc_ER_simulation.py.
(2) undirected ER networks: un_ER_numerical.py, un_ER_simulation.py.
(3) directed SF networks: direc_SF_numerical.py, direc_SF_simulation.py.
(4) undirected SF networks: un_SF_numerical.py, un_SF_simulation.py.
(5) DBLP scientist collaboration networks: Data file: dblp-v5-1k.txt, Python: hot_streaks.py.
(6) common functions: function_general.py

#2. Python enviroments.
(1) Implementation by Python 3.9 
(2) Dependency: networkx, scipy, numpy, statistics, matplotlib, pandas and function_general.

#3. Data schema
Notes:Journal name, author list, publication year, domain, citation number,h-index number, arnet-id.
-#*Transaction Management in Multidatabase Systems.
-#@Yuri Breitbart,Hector Garcia-Molina,Abraham Silberschatz
-#year1995
-#confModern Database Systems
-#citation16
-#index1
-#arnetid3

Due to the limitation of space, we merely upload a sample data file with 1000 lines. The whole data file is nearly 1GB downloaded from DBLP website.

#4. Inquiry:
For any problems, please contact hlsun84@mail.ustc.edu.cn and we feel happy to give kind support to make your work better.

#5.  Collaboration:
We are open-minded to work with any person who really loves percolation theory and applications in Complex Networks.


