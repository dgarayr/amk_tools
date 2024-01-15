#!/usr/bin/env python3
import RXReader as arx
import networkx as nx
import numpy as np
import sys
from pathlib import Path

## Calculation functions
def get_avg_k(graph):
    # Compute the average node degree of a graph
    deg_info = dict(graph.degree())
    Nn = graph.number_of_nodes()
    k = sum(list(deg_info.values()))/Nn
    #print('avg_k',k)
    return k

def get_avg_Lr(glist):
    # Compute average shortest path length for a list of Erdos-Renyi graphs
    # as in https://arxiv.org/pdf/cond-mat/0407098.pdf
    gamma = 0.5772
    Ng = len(glist)
    Nn = glist[0].number_of_nodes()
  
    kvals = [get_avg_k(gx) for gx in glist]

    sum_term = sum([1/np.log(k) for k in kvals if (k-1.0>1e-4)])
    Lr = 0.5 + (np.log(Nn) - gamma)/Ng * sum_term
    return Lr

def get_Lr(graph):
    # Compute shortest path length for a single Erdos-Renyi graph
    # as in https://arxiv.org/pdf/cond-mat/0407098.pdf
    gamma = 0.5772
    Nn = graph.number_of_nodes()
    Lr = 0.5 + (np.log(Nn) - gamma)/np.log(get_avg_k(graph)) 
    return Lr

def spl_random(glist):
    Ngraphs = len(glist)
    avg_c = 0
    avg_t = 0
    for ii in range(Ngraphs):
        avg_c += nx.average_clustering(glist[ii])
        avg_t += nx.transitivity(glist[ii])
    avg_c = avg_c / Ngraphs
    avg_t = avg_t / Ngraphs
    return avg_c,avg_t

fol = str(sys.argv[1])
net = "RXNet.cg"
path = Path(fol+"/RXNet.barrless")


data = arx.RX_parser(fol,net,check_connection=True)
if (path.is_file()):
    data_barrless = arx.RX_parser(fol,rxnfile="RXNet.barrless")
    joined_data = [data[ii]+data_barrless[ii] for ii in range(len(data))]
    data = joined_data
G2 = arx.RX_builder(fol,data)

# Basic graph info and p estimation
Nn2 = G2.number_of_nodes()
Ne2 = G2.number_of_edges()
Nem2 = Nn2*(Nn2-1)/2
p2 = Ne2/Nem2

#Create lists of w = 1000 graphs with Erdos-Renyi 
w = 1000
#print("Complete graph, ER")
rand_er2 = [nx.erdos_renyi_graph(Nn2,p2) for ii in range(w)]

# And do it for the complete one from RXNet.recalc + RXNet.barrless
#print("Erdos-Renyi, complete")
#print(get_avg_Lr(rand_er2))

#print("Actual, complete")
#print("L = %.4f" % nx.average_shortest_path_length(G2))
#print("="*42)


fw = open('rxn_stats.txt', 'w')
fw.write('#################################################### \n')
fw.write('#                                                  # \n')
fw.write('#        Properties of the reaction network        # \n')
fw.write('#                                                  # \n')
fw.write('#################################################### \n \n ')
fw.write("  Number of nodes = {0:>7} \n".format(Nn2))
fw.write("   Number of edges = {0:>7} \n".format(Ne2))
fw.write("   Average shortest path length of the current network             = {0:>7} \n".format(nx.average_shortest_path_length(G2)))
fw.write("   Average shortest path length of the equivalent random network   = {0:>7} \n".format(get_avg_Lr(rand_er2)))
fw.write("   Average clustering coefficient of the current network           = {0:>7} \n".format(nx.average_clustering(G2)))
fw.write("   Average clustering coefficient of the equivalent random network = {0:>7} \n".format(spl_random(rand_er2)[0]))
fw.write("   Transitivity of the current network                             = {0:>7} \n".format(nx.transitivity(G2)))
fw.write("   Transitivity of the equivalent random network                   = {0:>7} \n".format(spl_random(rand_er2)[1]))
fw.write("   Density of edges (edges/possible_edges)                         = {0:>7} \n".format(Ne2/Nem2))
fw.write("   Degree assortativity coefficient                                = {0:>7} \n".format(nx.degree_assortativity_coefficient(G2)))
fw.close()




