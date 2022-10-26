import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import glob
import scipy


def get_degree_distribution(G):
    # You can get the degree from every node as a dict.
    # We use list comprehension to pick out just the degrees, ignoring the keys, which are just the node names.
    # https://networkx.org/documentation/stable/reference/classes/generated/networkx.Graph.degree.html
    degrees = [d for _, d in G.degree]


    # We can compute the degree distribution by picking out the unique degrees and counting
    udegrees = np.unique(degrees)
    Pd = np.array([sum(degrees == ud) for ud in udegrees])
    Pd = Pd/sum(Pd)
    return udegrees, Pd


def plot_degree_dist(udegrees, Pd, show=False, fig_nr=None, network_name=None,legend_on=True):
    # Plot the results
    if fig_nr is not None:
        if fig_nr==0:
            plt.figure()
        else:
            plt.figure(fig_nr)

    if network_name is not None:
        plt.scatter(udegrees, Pd,label=network_name)
    else:
        plt.scatter(udegrees, Pd)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Degree $k$")
    plt.ylabel("$P(k)$")
    if legend_on:
        plt.legend(loc="best")
    # if network_name is not None:
        # plt.title(network_name)
    if show:
        plt.show()


def plot_clustering_coefficient(clustering_coeff, bins=None, show=False, fig_nr=None, network_name=None, legend_on=True, normalized=False, line_type= None):
    clustering_coeff = list(clustering_coeff)
    if fig_nr is not None:
        if fig_nr==0:
            plt.figure()
        else:
            plt.figure(fig_nr)
    if bins is None:
        # bins = np.arange(np.min(clustering_coeff),np.max(clustering_coeff))
        bins=np.linspace(0,1,10)
    if line_type is None:
        line_type = "-"
    ax = plt.gca()
    color = next(ax._get_lines.prop_cycler)["color"]
    if normalized:
        # increment and get the "props" cycle (and extract the color)
        (counts,bins) = np.histogram(clustering_coeff,bins=bins) #bins=np.linspace(0,1,10)
        plt.hist(bins[:-1], bins, weights=counts/max(counts), ls="dashed", lw=3, edgecolor=color, label=network_name,fc="None")

    else:
        if network_name is not None:
            (counts,bins) = np.histogram(clustering_coeff,bins=bins) #bins=np.linspace(0,1,10)
            plt.hist(bins[:-1], bins, weights=counts, histtype="step",lw=3, edgecolor=color, label=network_name,fc=color, alpha=.3, ls=line_type)
            # plt.hist(clustering_coeff,bins, ls="dashed", lw=3,label=network_name,fc="None")
        else:
            plt.hist(clustering_coeff,bins)

    plt.xlabel("Clustering coefficient")
    plt.ylabel("Number of nodes")
    plt.yscale("log")
    if legend_on:
        plt.legend(loc="best")
    # if network_name is not None:
        # plt.title(network_name)

    if show:
        plt.show()


def get_clustering_coefficient(G):
    clustering_coeff = nx.clustering(G).values()

    return clustering_coeff
    # print(G.degree)
    # print(nx.clustering(G))
    # print(clustering_coeff)
    # print(len(clustering_coeff))
    # print(G)


def get_shortest_path_lengths(G):
    spl = dict(nx.all_pairs_shortest_path_length(G))
    nodes = list(spl.keys())
    # spl_edges = []
    spl_list = []
    for i, source_node in enumerate(nodes):
        for target_node in spl[source_node]:
            if target_node not in nodes[:i+1]:
                spl_list.append(spl[source_node][target_node])
                # spl_edges.append({source_node, target_node})
    return spl_list



def plot_spl_hist(spl_list, bins=None, show=False, fig_nr=None, network_name=None, legend_on=True, line_type=None):
    if fig_nr is not None:
        if fig_nr==0:
            plt.figure()
        else:
            plt.figure(fig_nr)
    ax = plt.gca()
    color = next(ax._get_lines.prop_cycler)["color"]
    if bins is None:
        bins = np.linspace(0,10,10)
        # bins = np.arange(np.min(spl_list),np.max(spl_list)+1)
    if line_type is None: 
        line_type= "-"
    if network_name is not None:
        (counts,bins) = np.histogram(spl_list,bins=bins) #bins=np.linspace(0,1,10)
        plt.hist(bins[:-1], bins, weights=counts, histtype="step",ls=line_type,  lw=3, edgecolor=color, label=network_name,fc=color, alpha=.3)
        # plt.hist(spl_list,bins,label=network_name)
    else:
        plt.hist(spl_list,bins)
    if legend_on:
        plt.legend(loc="best")
    plt.xlabel("Shortest path length between nodes")
    plt.ylabel("Nr of node pairs")
    plt.yscale("log")
    # if network_name is not None:
        # plt.title(network_name)

    if show:
        plt.show()


def generate_random_network(nodes=300, prob=.03, filename= "../networks/my_random_network.txt", write_to_file=False):
    G = nx.fast_gnp_random_graph(nodes, prob)
    if write_to_file:
        with open(filename,"wb") as f:
            nx.write_edgelist(G_rand,f,delimiter=",")
    return G

def save_plot(filename, path="",extension=".eps", fig_nr=None, increment_filename=True):
    if fig_nr is not None:
        plt.figure(fig_nr)
    if increment_filename:
        num = []
        new_num = "0"
        matching_files = glob.glob(f"{path}{filename}*")
        print(matching_files)
        if (len(matching_files) > 0):
            max_num = -1
            for file in matching_files:
                fname = file.replace(path,"").split(".")[-2]
                digits = ''.join(x for x in fname[-3:] if x.isdigit())
                try:
                    num.append(int(digits))
                except KeyboardInterrupt:
                    pass
                except:
                    pass
            max_num = max(num)
            if ((max_num)>(-1)):
                new_num = str((max_num)+1)
            else:
                new_num = "0"
        plt.savefig(f"{path}{filename}{new_num}{extension}")


def test_plot_increment():
    plt.plot([1,2,3], [1,2,3])
    plt.figure()
    plt.plot([1,2,3], [-1,-2,-3])
    for i in range(5):
        save_plot("figtest", "../report/figs/", ".eps",fig_nr=i%2+1)

# folder_path =
network_paths = ["../networks/yeast_gene_net.txt", "../networks/human_kidney_protein_net.txt", "../networks/ecoli_metabolic_net.txt" ]
network_names = [ "Yeast gene network", "Human kidney protein", "E.Coli metabolic net"]
save_plots =True
draw_networks = False#True
line_types = ["-", "--", ":"]

plot_names = ["degree_dist","clustering_coeff", "spl"]

if __name__=="__main__":
    # Read the network file as an edgelist, as specified in instructions
    # Note that the delimiter is ","
    # https://networkx.org/documentation/stable/reference/readwrite/generated/networkx.readwrite.edgelist.read_edgelist.html
    # test_plot_increment()

    for i,network in enumerate(  network_paths ):
        G = nx.read_edgelist(network, delimiter=",")
        # G_rand = nx.fast_gnp_random_graph(300, .03)
        # network_names = 3*["Random Network"]
        # G = G_rand
        if draw_networks:
            plt.figure(i+10)
            nx.draw(G)
        print(G)
        print(network_names[i])

        unique_degrees,degree_dist=get_degree_distribution(G)
        plot_degree_dist(unique_degrees,degree_dist,fig_nr=1 , network_name=network_names[i])
        # if save_plots:
            # save_plot(network_names[i].replace(" ", "_").lower()+"_degree_dist_v","../report/figs/task1/")

        clustering_coeff=get_clustering_coefficient(G)
        plot_clustering_coefficient(clustering_coeff, fig_nr=2, network_name=network_names[i],line_type=line_types[i])
        # if save_plots:
            # save_plot(network_names[i].replace(" ", "_").lower()+"_clustering_coeff_v","../report/figs/task1/")
        spl=get_shortest_path_lengths(G)
        plot_spl_hist(spl,fig_nr=3, network_name=network_names[i], line_type=line_types[i])
        # if save_plots:
            # save_plot(network_names[i].replace(" ", "_").lower()+"_spl_v","../report/figs/task1/")
        print(np.unique(spl))

    if save_plots:
        for i in range(3):

            save_plot(plot_names[i].replace(" ", "_").lower()+"_v","../report/figs/task1/",fig_nr=i+1)


    plt.show()

