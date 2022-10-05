import numpy as np
import matplotlib.pyplot as plt
from numpy.random import rand
import time

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
    inds = []
    for i in range(nr_indices):
        inds.append(np.random.randint(0,60))
    return inds


def recruited_conversion(inds):
    # if (rand() < alpha):
    n1 = nucleosome_chain[inds[0]]
    n2 = nucleosome_chain[inds[1]]
    n1.update_methylation(n2.methylation)

def update_nr_methylated(nr_methylated,t,inds, nucleosome_chain):
    for i in range(-1,2):
        # nr_methylated[i+1, t] = nucleosome_chain.count(i)
        # return nr_methylated
        for n in nucleosome_chain:
            if n.methylation == i:
                nr_methylated[i+1, t] += 1
    return nr_methylated

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
    # plt.legend(loc="best")
    # plt.show()

def plot_methylated(nr_methylated, t):
    plt.figure()
    x = np.arange(0,t+1)/60
    plt.plot(x, nr_methylated[-1,:])
    plt.fill_between(x, nr_methylated[-1,:])
    plt.xlabel("Time [updates per nucleosome]")
    plt.ylabel("Number of methylated nucleosomes [1]")




L = 60 #nr of nucleosomes
tmax = 10000
filename = "data/task1_nr_methylated_v6_f6_t3e6.txt"
write_to_file =False

for F in [2, 4, 6]:
    # F = 4
    alpha =  F/(1 + F)
    nr_methylated = np.zeros([3,L*tmax])
    nucleosome_chain = init_nucleosome_chain(L)

    start_time = time.time()
    # main loop
    for t in range(L*tmax):
        inds = get_random_index(3)
        if (rand() < alpha):
            recruited_conversion(inds)
        if (rand() < (1-alpha)):
            nucleosome_chain[inds[2]].noisy_change()
        nr_methylated = update_nr_methylated(nr_methylated,t,inds, nucleosome_chain)

    plot_methylated(nr_methylated, t)
    mk_hist(nr_methylated[-1,:]-nr_methylated[0,:],1,f"F={F}")
    print("\n--- %s seconds ---" % (time.time() - start_time))

    if write_to_file:
        with open(filename, "a+") as f:
            np.savetxt(f, nr_methylated)


# mk_hist(nr_methylated[-1,:], boxwidth=1, label_string="\# methylated nucleosomes")

# loaded = np.loadtxt(filename)
# plt.figure()
# plt.plot(loaded[0,:])
# plt.figure()
# plt.plot(loaded[2,:])
# print(nr_methylated ==loaded[-3:, :])

plt.show()



