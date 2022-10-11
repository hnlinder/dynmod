import networkx as nx
import matplotlib


# filename = "../networks/my_small_random_network.txt"
filename = "../networks/my_random_network.txt"
G_rand = nx.fast_gnp_random_graph(10, .03)
# with open(filename,"wb") as f:
    # nx.write_edgelist(G_rand,f,delimiter=",")

# print(G_rand)


# G = nx.read_edgelist(filename, delimiter=",")
G = G_rand
# print(nx.shortest_path(G)["3"]["4"])

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

print(spl_list)
print(spl_edges)
nx.draw(G)
# print(len(spl_list))
print(f"Expected length: {G}\nMy length: {len(spl_list)}")












# for i, source_node in enumerate(nodes):
    # for target_node in list(nodes)[i:]:
        # # print(spl[node1][node2])
        # spl.append(nx.shortest_path_length(G,source_node,target_node))
        # print(spl)
