import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("tkagg")
import glob
import scipy


def get_degree_distribution(G):
    # You can get the degree from every node as a dict.
    # We use list comprehension to pick out just the degrees, ignoring the keys, which are just the node names.
    # https://networkx.org/documentation/stable/reference/classes/generated/networkx.Graph.degree.html
    degrees = [d for _, d in G.degree]
    av_degree = np.mean(degrees)


    # We can compute the degree distribution by picking out the unique degrees and counting
    udegrees = np.unique(degrees)
    Pd = np.array([sum(degrees == ud) for ud in udegrees])
    Pd = Pd/sum(Pd)
    return udegrees, Pd, av_degree


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


def plot_clustering_coefficient(clustering_coeff, bins=None, show=False, fig_nr=None, network_name=None, legend_on=True, normalized=False, line_type=None, see_through=False):
    clustering_coeff = list(clustering_coeff)
    if fig_nr is not None:
        if fig_nr==0:
            plt.figure()
        else:
            plt.figure(fig_nr)
    if line_type is None:
        line_type="-"
    alpha=1
    if see_through:
        alpha = .3

    ax = plt.gca()
    color = next(ax._get_lines.prop_cycler)["color"]
    if bins is None:
        # bins = np.arange(np.min(clustering_coeff),np.max(clustering_coeff))
        bins=np.linspace(0,1,10)
    if normalized:
        # increment and get the "props" cycle (and extract the color)
        (counts,bins) = np.histogram(clustering_coeff,bins=np.linspace(0,1,10)) #bins=np.linspace(0,1,10)
        plt.hist(bins[:-1], bins, weights=counts/max(counts), ls=line_type, lw=3, edgecolor=color, label=network_name,fc=color,alpha=alpha)

    else:
        plt.yscale("log")
        if network_name is not None:
            plt.hist(clustering_coeff,bins,histtype="step",ls=line_type, lw=3,label=network_name,edgecolor=color,fc=color,alpha=alpha)
            # plt.hist(clustering_coeff,bins, ls="dashed", lw=3,label=network_name,fc="None")
        else:
            plt.hist(clustering_coeff,bins)

    plt.xlabel("Clustering coefficient")
    plt.ylabel("Number of nodes")

    if legend_on:
        plt.legend(loc="best")
    # if network_name is not None:
        # plt.title(network_name)

    if show:
        plt.show()


def get_clustering_coefficient(G):
    clustering_coeff = nx.clustering(G).values()
    clust_arr = list(clustering_coeff)

    return clust_arr
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



def plot_spl_hist(spl_list, bins=None, show=False, fig_nr=None, network_name=None, legend_on=True, normalized=False, line_type=None, see_through=False):
    if fig_nr is not None:
        if fig_nr==0:
            plt.figure()
        else:
            plt.figure(fig_nr)
    # if bins is None:
        # bins = np.arange(np.min(spl_list),np.max(spl_list)+1)
    alpha=1
    if see_through:
        alpha = .3
    if line_type is None:
        line_type="-"
    ax = plt.gca()
    color = next(ax._get_lines.prop_cycler)["color"]
    bins=np.linspace(0,10,10)
    if normalized:
        # increment and get the "props" cycle (and extract the color)
        (counts,bins) = np.histogram(spl_list,bins=np.linspace(0,10,10)) #bins=np.linspace(0,1,10)
        plt.hist(bins[:-1], bins, weights=counts/nr_random_networks, ls=line_type, lw=3, edgecolor=color, label=network_name,fc=color,alpha=alpha)

    else:
        plt.yscale("log")
        if network_name is not None:
            plt.hist(spl_list,bins,histtype="step",ls=line_type, lw=3,label=network_name,edgecolor=color,fc=color,alpha=alpha)
        else:
            plt.hist(spl_list,bins)
    # if network_name is not None:
        # plt.hist(spl_list,bins,label=network_name)
    # else:
        # plt.hist(spl_list,bins)
    if legend_on:
        plt.legend(loc="best")
    plt.xlabel("Shortest path length between nodes")
    plt.ylabel("Nr of node pairs")
    # if network_name is not None:
        # plt.title(network_name)

    if show:
        plt.show()

def plot_hist(hist_list, bins=None,x_label=None, y_label=None, show=False, fig_nr=None, network_name=None, legend_on=True, normalized=False, line_type=None, see_through=True):
    if fig_nr is not None:
        if fig_nr==0:
            plt.figure()
        else:
            plt.figure(fig_nr)
    # if bins is None:
        # bins = np.arange(np.min(spl_list),np.max(spl_list)+1)
    alpha=1
    if see_through:
        alpha = .3
    if line_type is None:
        line_type="-"
    ax = plt.gca()
    color = next(ax._get_lines.prop_cycler)["color"]
    bins=np.linspace(0,10,10)
    if normalized:
        # increment and get the "props" cycle (and extract the color)
        (counts,bins) = np.histogram(hist_list,bins=np.linspace(0,max(hist_list),max(hist_list))) #bins=np.linspace(0,1,10)
        plt.hist(bins[:-1], bins, weights=counts/max(counts), ls=line_type, lw=3, edgecolor=color, label=network_name,fc=color,alpha=alpha)

    else:
        plt.yscale("log")
        if network_name is not None:
            plt.hist(hist_list,bins,ls=line_type, lw=3,label=network_name,edgecolor=color,fc=color,alpha=alpha)
        else:
            plt.hist(hist_list,bins)
    # if network_name is not None:
        # plt.hist(spl_list,bins,label=network_name)
    # else:
        # plt.hist(spl_list,bins)
    if legend_on:
        plt.legend(loc="best")
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    # if network_name is not None:
        # plt.title(network_name)

    if show:
        plt.show()



def save_plot(filename, path="",extension=".eps", fig_nr=None, increment_filename=True):
    if fig_nr is not None:
        plt.figure(fig_nr)
    if increment_filename:
        num = []
        new_num = "0"
        matching_files = glob.glob(f"{path}{filename}*")
        print(f"Saving plot: matching files: {matching_files}")
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

def generate_random_network(nodes=300, prob=.03, nr_networks=1, filename= "../networks/my_random_network.txt", write_to_file=False):
    G_arr = []
    for i in range(nr_networks):
        G_arr.append(nx.fast_gnp_random_graph(nodes, prob))
        if write_to_file:
            with open(filename,"wb") as f:
                nx.write_edgelist(G_arr[i],f,delimiter=",")
    return G_arr

def get_network_characteristics(G):
    nr_nodes = G.number_of_nodes()
    nr_edges = G.number_of_edges()
    return nr_nodes, nr_edges

def get_uniques(list_of_arrs):
    all_elements = []
    for arr in list_of_arrs:
        all_elements.append(list(arr))
    return np.unique(all_elements)

def mean_degree_dist(unique_degree_arr, degree_dist_arr):
    b = []
    for a in unique_degree_arr:
        b.append(list(a))
    uds = np.unique(b)
    ud_dict = {}
    deg_dict = {}
    av_degree_dist = []
    for ud in uds:
        ud_dict[ud] = 0
        deg_dict[ud] = 0

    for i,(ud_arr,deg_arr) in enumerate(zip(unique_degree_arr,degree_dist_arr)):
        for ud,deg_dist in zip(ud_arr, deg_arr):
            ud_dict[ud] = ud_dict[ud] + 1
            deg_dict[ud] = deg_dict[ud] + deg_dist

    for ud in ud_dict.keys() :
        av_degree_dist.append(deg_dict[ud]/ud_dict[ud])
    return av_degree_dist



network_paths = ["../networks/yeast_gene_net.txt", "../networks/human_kidney_protein_net.txt", "../networks/ecoli_metabolic_net.txt" ]
network_names = [ "Yeast gene network", "Human kidney protein", "E.Coli metabolic net"]

save_plots = True
draw_networks = False
nr_random_networks = 30
if __name__=="__main__":


    # For the three original networks
    for i,network in enumerate(  network_paths ):
        G = nx.read_edgelist(network, delimiter=",")
        nr_nodes, nr_edges = get_network_characteristics(G)
        # G_rand = nx.fast_gnp_random_graph(300, .03)
        # network_names = 3*["Random Network"]
        # G = G_rand
        if draw_networks:
            print(f"At plot, i={i}, type={type(i)}")
            plt.figure(int(i+10))
            nx.draw(G)
        print(G)
        print(network_names[i])

        unique_degrees,degree_dist,av_degree=get_degree_distribution(G)
        # plot_degree_dist(unique_degrees,degree_dist,fig_nr=1 , network_name=network_names[i])
        # if save_plots:
            # save_plot(network_names[i].replace(" ", "_").lower()+"_degree_dist_v","../report/figs/task1/")

        clustering_coeff=get_clustering_coefficient(G)
        plot_clustering_coefficient(clustering_coeff, fig_nr=i+1, network_name=network_names[i])
        # if save_plots:
            # save_plot(network_names[i].replace(" ", "_").lower()+"_clustering_coeff_v","../report/figs/task2/")
        spl=get_shortest_path_lengths(G)
        plot_spl_hist(spl,fig_nr=i+4, network_name=network_names[i], normalized=False)
        # if save_plots:
            # save_plot(network_names[i].replace(" ", "_").lower()+"_spl_v","../report/figs/task2/")
        print(np.unique(spl))

        print("Random networks\n--------------------------------------------")
        random_networks = generate_random_network(nodes=nr_nodes,prob=(av_degree/(nr_nodes-1)),nr_networks=nr_random_networks)
        clustering_coeff_arr = np.zeros([nr_random_networks,nr_nodes])
        unique_degrees_arr=[]
        degree_dist_arr = []
        spl_arr = []
        for j,G in enumerate(random_networks):
            print(G)
            print(network_names[i])

            unique_degrees,degree_dist,av_degree=get_degree_distribution(G)
            unique_degrees_arr.append(unique_degrees)
            degree_dist_arr.append(degree_dist)
            clustering_coeff_arr[j,:]=get_clustering_coefficient(G)
            spl_arr.append(get_shortest_path_lengths(G))
            print(f"Unique spl: {np.unique(spl)}")

        clustering_coeff = np.zeros([nr_nodes,1])
        unique_degrees = []
        degree_dist = []
        spl = []
        clustering_coeff = np.mean(clustering_coeff_arr, axis=0)
        # plt.figure(99)
        # plt.plot(clustering_coeff,label=network_names[i])
        # for k in range(nr_nodes):
            # clustering_coeff[k]= np.mean(clustering_coeff_arr[:,k])
            # unique_degrees = get_uniques(unique_degrees_arr)
            # print(f"degree dist_arr : {degree_dist_arr}")
            # degree_dist.append(np.mean(degree_dist_arr[k]))
        # for k in range(len(random_networks)):
            # spl.append(np.mean(spl_arr[k]))
            # spl.append(spl_arr[k])
        for arr in spl_arr:
            for val in arr:
                spl.append(val)
        # spl = spl_arr


        # unique_degrees = np.mean(unique_degrees_arr,axis=0)
        # degree_dist = np.mean(degree_dist_arr,axis=0)
        # plots
        # plot_degree_dist(unique_degrees,degree_dist,fig_nr=1 , network_name=f"{network_names[i]} random")
        # if save_plots:
            # save_plot(network_names[i].replace(" ", "_").lower()+"_degree_dist_v","../report/figs/task1/")
        plot_clustering_coefficient(clustering_coeff, fig_nr=i+1, network_name=f"{network_names[i]} random", line_type="--")
        if save_plots:
            save_plot(network_names[i].replace(" ", "_").lower()+f"_clustering_coeff_w_{nr_random_networks}_rand_v","../report/figs/task2/")
        plot_spl_hist(spl,fig_nr=i+4, network_name=f"{network_names[i]} random",normalized=False,line_type="--")
        if save_plots:
            save_plot(network_names[i].replace(" ", "_").lower()+f"_spl_w_{nr_random_networks}_rand_v","../report/figs/task2/")


    plt.show()

