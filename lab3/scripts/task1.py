import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

# Read the network file as an edgelist, as specified in instructions
# Note that the delimiter is ","
# https://networkx.org/documentation/stable/reference/readwrite/generated/networkx.readwrite.edgelist.read_edgelist.html
G = nx.read_edgelist("../networks/human_kidney_protein_net.txt", delimiter=",")

def get_degrees(G):
    # You can get the degree from every node as a dict.
    # We use list comprehension to pick out just the degrees, ignoring the keys, which are just the node names.
    # https://networkx.org/documentation/stable/reference/classes/generated/networkx.Graph.degree.html
    degrees = [d for _, d in G.degree]


    # We can compute the degree distribution by picking out the unique degrees and counting
    udegrees = np.unique(degrees)
    Pd = np.array([sum(degrees == ud) for ud in udegrees])
    Pd = Pd/sum(Pd)

    # Plot the results
    # plt.scatter(udegrees, Pd)
    # plt.xscale("log")
    # plt.yscale("log")
    # plt.xlabel("Degree $k$")
    # plt.ylabel("$P(k)$")
    # plt.show()


def get_clustering(G):
    clustering_coeff = nx.clustering(G).values()

    x = np.arange(0,len(clustering_coeff))
    plt.scatter(x,clustering_coeff)
    plt.show()
    # print(G.degree)
    # print(nx.clustering(G))
    # print(clustering_coeff)
    # print(len(clustering_coeff))
    # print(G)


def get_shortest_path_lengths(G):
    spl = dict(nx.all_pairs_shortest_path_length(G))
    print(spl)
    # spl = []
    # nodes = G.nodes()
    nodes = list(spl.keys())

    spl_edges = []
    spl_list = []
    for i, source_node in enumerate(nodes):
        for target_node in spl[source_node]:
            if target_node not in nodes[:i+1]:
                spl_list.append(spl[source_node][target_node])
                # spl_edges.append({source_node, target_node})
    return spl_list



if __name__=="__main__":

    get_clustering(G)
