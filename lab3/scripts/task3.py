import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import glob
# import scipy


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


def plot_clustering_coefficient(clustering_coeff, bins=None, show=False, fig_nr=None, network_name=None, legend_on=True, normalized=True, line_type=None):
    clustering_coeff = list(clustering_coeff)
    if fig_nr is not None:
        if fig_nr==0:
            plt.figure()
        else:
            plt.figure(fig_nr)
    if line_type is None:
        line_type="-"

    ax = plt.gca()
    color = next(ax._get_lines.prop_cycler)["color"]
    # if bins is None:
        # bins = np.arange(np.min(clustering_coeff),np.max(clustering_coeff))
    if normalized:
        # increment and get the "props" cycle (and extract the color)
        (counts,bins) = np.histogram(clustering_coeff,bins=np.linspace(0,1,10)) #bins=np.linspace(0,1,10)
        plt.hist(bins[:-1], bins, weights=counts/max(counts), ls=line_type, lw=3, edgecolor=color, label=network_name,fc=color,alpha=0.3)

    else:
        if network_name is not None:
            plt.hist(clustering_coeff,bins, ls="dashed", lw=3,label=network_name,fc="None")
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



def plot_spl_hist(spl_list, bins=None, show=False, fig_nr=None, network_name=None, legend_on=True, normalized=False, line_type=None, see_through=True):
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
        plt.hist(bins[:-1], bins, weights=counts/max(counts), ls=line_type, lw=3, edgecolor=color, label=network_name,fc=color,alpha=alpha)

    else:
        plt.yscale("log")
        if network_name is not None:
            plt.hist(spl_list,bins,ls=line_type, lw=3,label=network_name,edgecolor=color,fc=color,alpha=alpha)
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
    new_num = ""
    if increment_filename:
        num = []
        matching_files = glob.glob(f"{path}{filename}*")
        # print(matching_files)
        if (len(matching_files) > 0):
            max_num = -1
            for file in matching_files:
                fname = file.replace(path,"").split(".")[-2]
                digits = ''.join(x for x in fname[-3:] if x.isdigit()) #crashes if len(fname)<3
                try: #I reeeeaaally really
                    num.append(int(digits))
                except KeyboardInterrupt: #Fucking hate 
                    pass
                except: # This Try-Except
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


def remove_rand_fraction(G,f):
    nr_nodes = G.number_of_nodes()
    nr_nodes_to_remove = int(nr_nodes*f)
    nodes_to_remove = np.random.choice(G.nodes(), nr_nodes_to_remove, replace=False)
    G.remove_nodes_from(nodes_to_remove)


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
        # G.remove_nodes_from(top_nodes)
        G.remove_nodes_from([node_names[n] for n in top_nodes])


def get_network_diameter(spl_list):
    return max(spl_list)






network_paths = ["../networks/yeast_gene_net.txt", "../networks/human_kidney_protein_net.txt", "../networks/ecoli_metabolic_net.txt" ]
network_names = [ "Yeast gene network", "Human kidney protein", "E.Coli metabolic net"]



save_plots = True
draw_networks =False
nr_runs_per_f = 20
nr_fs = 6
fs = np.linspace(0,.25,nr_fs)

if __name__=="__main__":
    # For the three original networks
    for i,network in enumerate(  network_paths ):
        # G_og = nx.read_edgelist(network, delimiter=",")
        nd = np.zeros([nr_fs,1])
        for j,f in enumerate(fs):
            for k in range(nr_runs_per_f):
                print(f"f:{f}, k:{k}")
                G = nx.read_edgelist(network, delimiter=",")
                # f = .1
                # remove_rand_fraction(G,f)
                # print(f"f:{f}")
                # print(G)
                remove_top_fraction(G,f)
                print(f"After:  {G}")
                spl=get_shortest_path_lengths(G)
                nd_tmp =  get_network_diameter(spl)
                print(f"Network diameter:{ nd_tmp}")
                nd[j] += nd_tmp
                print(f"Network diameter: {nd}")
            nd[j] = nd[j]/nr_runs_per_f
        plt.plot(fs, nd, label=network_names[i])
    plt.xlabel("Fraction removed")
    plt.ylabel("Network diameter")
    plt.legend(loc="best")
    if save_plots:
        save_plot(network_names[i].replace(" ", "_").lower()+f"_top_rm_{nr_runs_per_f}_runs_per_f_v","../report/figs/task3/")
    

        ## G_rand = nx.fast_gnp_random_graph(300, .03)
        ## network_names = 3*["Random Network"]
        ## G = G_rand
        # if draw_networks:
            # plt.figure(i+10)
            # nx.draw(G)
        # print(G)
        # print(network_names[i])

        # unique_degrees,degree_dist,av_degree=get_degree_distribution(G)
        # plot_degree_dist(unique_degrees,degree_dist,fig_nr=1 , network_name=network_names[i])
        # if save_plots:
            # save_plot(network_names[i].replace(" ", "_").lower()+"_degree_dist_v","../report/figs/task1/")

        # clustering_coeff=get_clustering_coefficient(G)
        # plot_clustering_coefficient(clustering_coeff, fig_nr=2, network_name=network_names[i])
        # if save_plots:
            # save_plot(network_names[i].replace(" ", "_").lower()+"_clustering_coeff_v","../report/figs/task1/")
        # spl=get_shortest_path_lengths(G)
        # plot_spl_hist(spl,fig_nr=3, network_name=network_names[i])
        # if save_plots:
            # save_plot(network_names[i].replace(" ", "_").lower()+"_spl_v","../report/figs/task1/")
        # print(np.unique(spl))

        # print("Random networks")
        # nr_random_networks = 2
        # random_networks = generate_random_network(nodes=nr_nodes,prob=(av_degree/(nr_nodes-1)),nr_networks=nr_random_networks)
        # clustering_coeff_arr = np.zeros([nr_random_networks,nr_nodes])
        # unique_degrees_arr=[]
        # degree_dist_arr = []
        # spl_arr = []
        # for j,G in enumerate(random_networks):
            # print(G)
            # print(network_names[i])

            # unique_degrees,degree_dist,av_degree=get_degree_distribution(G)
            # unique_degrees_arr.append(unique_degrees)
            # degree_dist_arr.append(degree_dist)
            # clustering_coeff_arr[j,:]=get_clustering_coefficient(G)
            # spl_arr.append(get_shortest_path_lengths(G))
            # print(np.unique(spl))

        # clustering_coeff = np.zeros([nr_nodes,1])
        # unique_degrees = []
        # degree_dist = []
        # spl = []
        # for k in range(nr_nodes):
            # clustering_coeff[k]= np.mean(clustering_coeff_arr[:,k])
        # for item in unique_degrees_arr:
            # unique_degrees.append(np.mean(item))
        # for item in degree_dist_arr:
            # degree_dist.append(np.mean(item))
        # for item in spl_arr:
            # spl.append(np.mean(item))


        # # unique_degrees = np.mean(unique_degrees_arr,axis=0)
        # # degree_dist = np.mean(degree_dist_arr,axis=0)
        # # plots
        # plot_degree_dist(unique_degrees,degree_dist,fig_nr=1 , network_name=f"{network_names[i]} random")
        # if save_plots:
            # save_plot(network_names[i].replace(" ", "_").lower()+"_degree_dist_v","../report/figs/task1/")
        # plot_clustering_coefficient(clustering_coeff, fig_nr=2, network_name=f"{network_names[i]} random")
        # if save_plots:
            # save_plot(network_names[i].replace(" ", "_").lower()+"_clustering_coeff_v","../report/figs/task1/")
        # plot_spl_hist(spl,fig_nr=3, network_name=f"{network_names[i]} random")
        # if save_plots:
            # save_plot(network_names[i].replace(" ", "_").lower()+"_spl_v","../report/figs/task1/")


    plt.show()

