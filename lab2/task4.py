import numpy as np
import matplotlib.pyplot as plt
from numpy.random import rand
import time
import random

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
def zero_nr_methylated(nr_methylated, t):
    nr_methylated[:,t] = 0
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

def cell_division(nucleosome_chain):
    print("Cell division")
    for n in nucleosome_chain:
        n.methylation = n.methylation*(rand() < .5)
def plot_methylated(nr_methylated, t, cell_division_index):
    plt.figure()
    x = np.arange(0,t+1)/60
    plt.vlines(cell_division_index/60, 0,60,"r", label="Cell division")
    plt.plot(x, nr_methylated[-1,:])
    plt.fill_between(x, nr_methylated[-1,:])
    plt.xlabel("Time [updates per nucleosome]")
    plt.ylabel("Number of methylated nucleosomes [1]")
    plt.legend(loc="best")
    # plt.title(f"Methylated, gamma = {gamma}")




# def save_array_to_file

L = 60 #nr of nucleosomes
tmax = 10000
cell_div_counter = 0
filename = "data/dummy.txt"
# filename = "data/dummy.txt"
model2 = False
write_to_file = False
# gamma = -2
gammalist = [-1] #np.linspace(-1,-3,6)
cell_division_index = np.arange(int(L*tmax/11), int(L*tmax*(1-1/11)+1),int(L*tmax*(1-2/11)/10+1)) #10 evenly spaced cell divisions along t

print(cell_division_index)

if write_to_file:
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
            if (cell_div_counter < 10) and ( t == cell_division_index[cell_div_counter] ) :
                cell_division(nucleosome_chain)
                nr_methylated = zero_nr_methylated(nr_methylated,t)
                nr_methylated = update_nr_methylated(nr_methylated,t,inds, nucleosome_chain)
                print(cell_div_counter)
                cell_div_counter += 1



        plot_methylated(nr_methylated, t, cell_division_index)

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



