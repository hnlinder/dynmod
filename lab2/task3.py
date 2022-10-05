import numpy as np
import matplotlib.pyplot as plt
from numpy.random import rand
import time
import random

plt.rcParams.update({'font.size':12})
class Nucleosome():
    def __init__(self, methylation = 0):
        self.methylation = methylation

    # Here I define methylated = 1, unmethylated = 0, acethylated = -1
    def update_methylation(self, reference):
        if reference != self.methylation:
            self.methylation += reference

    def __neighbor_states(self):
        if abs(self.methylation):
            return 0
        else:
            return [-1,1][rand()<.5]

    def noisy_change(self):
        # if (rand() < (1-alpha)):
        self.methylation = self.__neighbor_states()

def init_nucleosome_chain(L):
    nucleosome_chain = []
    for i in range(L):
        # random distribution of methylation
        nucleosome_chain.append(Nucleosome(np.random.randint(-1, 2)))
    return nucleosome_chain


def get_random_index(nr_indices):
    # inds = []
    # for i in range(nr_indices):
    inds = random.sample(range(60),nr_indices)
    return inds



def recruited_conversion(inds, model2):
    # if (rand() < alpha):
    n1 = nucleosome_chain[inds[0]]
    n2 = nucleosome_chain[inds[1]]
    n3 = nucleosome_chain[inds[2]]
    if model2:
        if (n2.methylation==n3.methylation):
            n1.update_methylation(n2.methylation)
    else:
        n1.update_methylation(n2.methylation)



def update_nr_methylated(nr_methylated,t,inds, nucleosome_chain):
    for i in range(-1,2):
        # nr_methylated[i+1, t] = nucleosome_chain.count(i)
        # return nr_methylated
        for n in nucleosome_chain:
            if (n.methylation == i):
                nr_methylated[i+1, t] += 1
    return nr_methylated

def get_distances(ind, L):
    return np.abs(np.arange(0,L))


def get_index_of_n2(ind, gamma, L=60):
    d = np.abs(np.arange(0,L)-ind)
    d[ind] = 1
    P = np.float_power(d,gamma)
    d[ind] = 0
    Pc = np.cumsum(P)
    r = np.random.uniform(0,1)*Pc[-1]
    n2_ind = np.argmax(Pc>=r)
    return n2_ind



def plot_methylated(nr_methylated, t, legend_str):
    plt.figure()
    x = np.arange(0,t+1)/60
    plt.plot(x, nr_methylated[-1,:], label=legend_str)
    plt.fill_between(x, nr_methylated[-1,:])
    plt.xlabel("Time [updates per nucleosome]")
    plt.ylabel("Number of methylated nucleosomes [1]")
    plt.legend(loc="best")
# def save_array_to_file
def mk_hist(arr, boxwidth, label_string):
    # histstop = .4
    plt.figure()
    histstop = 60
    boxwidth = 1
    nrboxes = histstop/boxwidth
    bins = np.linspace(-histstop, histstop, 2*int(nrboxes))
    # bins = [0, .02, .04,.06]
    plt.hist(arr, bins, label=label_string)
    plt.xlabel("M-A")
    plt.ylabel("Nr of nucleosomes in state")
    plt.legend(loc="best")


# def save_array_to_file

L = 60 #nr of nucleosomes
tmax = 10000
filename = "data/dummy.txt"
model2 = False
gamma = -2
gammalist = np.linspace(-1,-3,6)
write_to_file = False


if  write_to_file:
    with open(filename, "a+") as f:
        np.savetxt(f,gammalist)
for gamma in gammalist:
    for F in [6]:#[2, 4, 6]:
        # F = 4
        alpha =  F/(1 + F)
        nr_methylated = np.zeros([3,L*tmax])
        nucleosome_chain = init_nucleosome_chain(L)

        start_time = time.time()
        # main loop
        for t in range(L*tmax):
            inds = get_random_index(3)
            inds[1] = get_index_of_n2(inds[0], gamma, L)
            # if not abs(nucleosome_chain[inds[0]].methylation): #only if this nucleosome is unmethylated, ie n==0
            if (rand() < alpha):
                recruited_conversion(inds, model2)
            if (rand() < (1-alpha)):
                nucleosome_chain[inds[-1]].noisy_change()
            nr_methylated = update_nr_methylated(nr_methylated,t,inds, nucleosome_chain)
            # else:
                # nr_methylated[:,t] = nr_methylated[:,t-1]



        # plt.figure()
        # plt.plot(nr_methylated[0,:])
        # plt.title(f"Acethylated, gamma = {gamma}")
        # plt.figure()
        # plt.plot(nr_methylated[2,:])
        # plt.title(f"Methylated, gamma = {gamma}")
        legend_str = fr"$\gamma = {gamma}$"
        plot_methylated(nr_methylated, t, legend_str)
        mk_hist(nr_methylated[-1,:]-nr_methylated[0,:],1,legend_str)

        print("\n--- %s seconds ---" % (time.time() - start_time))

        if write_to_file:
            with open(filename, "a+") as f:
                np.savetxt(f, nr_methylated)

# loaded = np.loadtxt(filename)
# plt.figure()
# plt.plot(loaded[0,:])
# plt.figure()
# plt.plot(loaded[2,:])
# print(nr_methylated ==loaded[-3:, :])

plt.show()



