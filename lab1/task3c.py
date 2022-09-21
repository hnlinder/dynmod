import numpy as np
import matplotlib.pyplot as plt
from numpy.random import rand
from time import sleep

L  = 30
# alpha = .9
# beta = 2
q = 1

Nr = 4
tmax = 3600 * 100
dt = 1/2

kp = 1/600
kd = 1/1800
Np0 = 0
kd_rib = 1/1800
kd_mrna = 1/300
kp_mrna = 1/600

meanlength = int(tmax/2)
mean = np.empty(10)
tarr = np.linspace(0, tmax-1, int((tmax-1)/dt))


class mRNA:
    def __init__(self, lattice = np.zeros(L), ribs = 0):
        self.lattice = lattice
        self.ribs = ribs

    # one time step of the ribosomes
    def rib_step(self, proteins_produced, Nr, alpha, beta ):
        rib_index = np.where(self.lattice>0)[0]
        np.random.shuffle(rib_index)
        for i in rib_index:
            if i != L-1:
                if ((rand() < dt * q) and valid_move(self.lattice,i)):
                    self.lattice[i+1] += 1
                    self.lattice[i] -= 1
        if ((self.lattice[0]==0) and (rand() < dt * alpha * Nr)):
            self.lattice[0] += 1
            Nr -=1
        if ((self.lattice[L-1]>0) and (rand() < dt * beta)):
            self.lattice[L-1] -= 1
            proteins_produced += 1
            # print("MADE PROTEIN")
            Nr +=1
        if (rand() < 1*dt * kd * proteins_produced):
            # print("REMOVED PROTEIN")
            proteins_produced -= 1

        return proteins_produced, Nr


# class LM:
    # def __init__(self, mRNA = mRNA(), matr = []):
        # self.mRNA = mRNA
        # self.matr = matr

    # def create_mRNA(self):
        # self.matr.append(mRNA(np.zeros(L)))

    # def remove_mRNA(self, index):
        # ribs = np.sum(self.matr[index].lattice)
        # del(self.matr[index])
        # return ribs

class LM:
    def __init__(self, mRNA = mRNA(), matr = [], nr_mRNAs = np.zeros(len(tarr))):
        self.mRNA = mRNA
        self.matr = matr
        self.nr_mRNAs = nr_mRNAs

    def create_mRNA(self, t):
        self.matr.append(mRNA(np.zeros(L)))
        self.nr_mRNAs[t:] += 1

    def remove_mRNA(self, index, t):
        ribs = np.sum(self.matr[index].lattice)
        del(self.matr[index])
        self.nr_mRNAs[t:] -= 1
        return ribs

def valid_move(lattice, i):
    return (lattice[i+1] == 0) and (lattice[i] > 0)


# lm.create_mRNA(0)
# main loop
graph_type = ["-", "-."]
prot_prod_decay = 0
alpha_arr = [.9]#np.linspace(0,1, 10):
beta_arr = [2]# [.25, .5]
for prot_prod_decay in [0,1]:
    for alpha in alpha_arr:
        for beta in beta_arr:
            lm = LM()
            lm.create_mRNA(0)
            ind = 0
            proteins_produced = np.zeros(len(tarr))
            for t in range(len(tarr)):
                proteins_produced[ind] = proteins_produced[ind - 1]
                # make new mrna
                if (prot_prod_decay*kp_mrna * dt > rand())  :
                    lm.create_mRNA(t)
                # Do a time step with ribosomes for each mRNA
                for index, mrna in enumerate(lm.matr):
                    proteins_produced[ind], Nr = mrna.rib_step(proteins_produced[ind], Nr, alpha, beta)
                    # remove mrna
                    if (prot_prod_decay*kd_mrna * dt * lm.nr_mRNAs[t] > rand()):
                        # print(lm.nr_mRNAs[-1])
                        Nr += lm.remove_mRNA(index, t)
                    # print(mrna.lattice)
                # sleep(.1)

                # print(Nr)
                ind +=1
            print(f"mean nr MRNAs: {np.mean(lm.nr_mRNAs)}")
            # sleep(.5)
            # plt.figure()
            # plt.title(f"alpha : {alpha}, beta = {beta}")
            # plt.ylim([0, 130])
            if prot_prod_decay:
                legend_lable = "Multiple mRNAs"
            else :
                legend_lable = "One mRNA"
            plt.plot(tarr, proteins_produced,graph_type[prot_prod_decay], label=legend_lable)

    fano = np.var(proteins_produced)/np.mean(proteins_produced)
    print(f"Fano factor : {fano}")
plt.xlabel("Time [s]")
plt.ylabel("Proteins produced [1]")
plt.legend(loc="best")
plt.show()
