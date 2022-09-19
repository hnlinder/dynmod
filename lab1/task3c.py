import numpy as np
import matplotlib.pyplot as plt
from numpy.random import rand
from time import sleep

L  = 30
# alpha = .9
# beta = 2
q = .3

Nr = 14
tmax = 3600 * 100
dt = 1/2

kp = 1/60
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


class LM:
    def __init__(self, mRNA = mRNA(), matr = []):
        self.mRNA = mRNA
        self.matr = matr

    def create_mRNA(self):
        self.matr.append(mRNA(np.zeros(L)))

    def remove_mRNA(self, index):
        ribs = np.sum(self.matr[index].lattice)
        del(self.matr[index])
        return ribs


def valid_move(lattice, i):
    return (lattice[i+1] == 0) and (lattice[i] > 0)


lm = LM()
lm.create_mRNA()
# main loop
for alpha in np.linspace(.2, 1, 5):
    for beta in [.25, .5]:
        ind = 0
        proteins_produced = np.zeros(len(tarr))
        for t in range(len(tarr)):
            proteins_produced[ind] = proteins_produced[ind - 1]
            # make new mrna
            if (1*kp_mrna * dt > rand())  :
                lm.create_mRNA()
            # Do a time step with ribosomes for each mRNA
            for index, mrna in enumerate(lm.matr):
                proteins_produced[ind], Nr = mrna.rib_step(proteins_produced[ind], Nr, alpha, beta)
                # remove mrna
                if (1*kd_mrna * dt > rand()):
                    Nr += lm.remove_mRNA(index)
                # print(mrna.lattice)
            # sleep(.1)

            # print(Nr)
            ind +=1
        plt.figure()
        plt.title(f"alpha : {alpha}, beta = {beta}")
        plt.ylim([0, 130])
        plt.plot(tarr, proteins_produced)

plt.show()
