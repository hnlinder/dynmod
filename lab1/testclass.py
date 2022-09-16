import numpy as np
import matplotlib.pyplot as plt
from numpy.random import rand
from time import sleep
L = 30
class mRNA:
    def __init__(self, lattice = np.zeros(L), ribs = 0):
        self.lattice = lattice
        self.ribs = ribs

    # one time step of the ribosomes
    def rib_step(self, proteins_produced, ind, alpha, beta, dt, lattice, Nr):
        proteins_produced[ind] = proteins_produced[ind - 1]
        rib_index = np.where(lattice>0)[0]
        # print(rib_index)
        np.random.shuffle(rib_index)
        for i in rib_index:
            if i != L-1:
                if ((rand() < dt * q) and valid_move(lattice,i)):
                    lattice[i+1] += 1
                    lattice[i] -= 1
        if ((lattice[0]==0) and (rand() < dt * alpha * Nr)):
            lattice[0] += 1
            Nr -=1
        if ((lattice[L-1]>0) and (rand() < dt * beta)):
            lattice[L-1] -= 1
            proteins_produced[ind] += 1
            Nr +=1
        if (rand() < dt * kd * proteins_produced[ind]):
            proteins_produced[ind] -= 1

        return proteins_produced, lattice, Nr





class LM:
    def __init__(self, mRNA = mRNA(), matr = []):
        self.mRNA = mRNA
        self.matr = matr

    def create_mRNA(self):
        self.matr.append(mRNA(np.zeros(L)))

    def remove_mRNA(self, index):
        ribs = self.matr[index].ribs
        del(self.matr[index])
        return ribs


lat1 =  mRNA(np.zeros(3))
# print(lat1.lattice)
lat2 =  mRNA(np.zeros(3))
lm1 = LM()
lm1.create_mRNA()
# print(lm1.matr[0].lattice)
lm1.create_mRNA()
print(f"length : {len(lm1.matr)}")
# print(lm1.matr[1].lattice)
lm1.matr[1].ribs += 1
# nr_ribs = lm1.remove_mRNA(1)
# print(nr_ribs)
# print(lm1.matr[1].lattice)
lm1.matr[1].lattice = np.ones(L)
lm1.remove_mRNA(1)
print("after remove:")
print(f"length : {len(lm1.matr)}")

for index, mrna in enumerate(lm1.matr):
    print(mrna.lattice)
    print(index)
    print(mrna.ribs)

