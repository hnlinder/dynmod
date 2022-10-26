import networkx as nx
import numpy as np



def remove_fraction(G,f):
    nr_nodes = G.number_of_nodes()
    # nodes_to_remove = np.random.choice(nr_nodes, int(nr_nodes*f), replace = False)
    # G.remove_nodes_from(nodes_to_remove)
    nr_nodes_to_remove = int(nr_nodes*f)
    print(nr_nodes_to_remove)

    nodes_to_remove = np.random.choice(G.nodes(),nr_nodes_to_remove,replace=False)
    # G.remove_nodes_from(nodes_to_remove)
    for n in nodes_to_remove:
        # print(f"G:{G}")
        # print(n)
        G.remove_node(n)

def generate_random_network(nodes=300, prob=.03, nr_networks=1, filename= "../networks/my_random_network.txt", write_to_file=False):
    G_arr = []
    for i in range(nr_networks):
        G_arr.append(nx.fast_gnp_random_graph(nodes, prob))
        if write_to_file:
            with open(filename,"wb") as f:
                nx.write_edgelist(G_arr[i],f,delimiter=",")
    return G_arr

def remove_top_fraction(G,f):
    if f==0:
        return 
    else:
        nr_nodes = G.number_of_nodes()
        degrees = [d for _, d in G.degree]
        nr_nodes_to_remove = int(nr_nodes*f)
        top_nodes = sorted(range(len(degrees)), key=lambda x: degrees[x])[-nr_nodes_to_remove:]
        # nodes_to_remove = np.random.choice(G.nodes(), nr_nodes_to_remove, replace=False)
        node_names = list(G.nodes())

        print([node_names[n] for n in top_nodes])
        G.remove_nodes_from([node_names[n] for n in top_nodes])#node_names[top_nodes])

network = "../networks/yeast_gene_net.txt"
if __name__=="__main__":

    nr_nodes= 100
    # G = nx.read_edgelist(network,delimiter=",")
    G = generate_random_network(10)
    G = G[0]
    # G = nx.path_graph(nr_nodes)
    print(G)
    f = .0

    remove_top_fraction(G,f)

    print(G)
