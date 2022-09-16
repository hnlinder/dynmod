import numpy as np
import matplotlib.pyplot as plt
from numpy.random import rand
from time import sleep
# print("hej")
# class mRNA:
    # def __init__(self, lattice):
        # self.lattice = lattice


# class LM:
    # def __init__(self, mRNA):
        # self.mRNA = mRNA

# lat1 =  mRNA(np.zeros(3))
# lm1 = LM(lat1)


# # print(lat1.lattice)
# print(lm1.mRNA.lattice)



L  = 30

alpha = .9
beta = 2
q = .3



Nr = 14
tmax = 3600 * 100
dt = 1/2

kp = 1/60
kd_rib = 1/1800
kd_mrna = 1/300
kp_mrna = 1/600
Np0 = 0
lattice  = np.zeros(L)

meanlength = int(tmax/2)
mean = np.empty(10)
tarr = np.linspace(0, tmax-1, int((tmax-1)/dt))
# print(tarr)


# def move_rib(lattice, i):
    # lattice[i] -=1
    # lattice[i+1] +=1
    # return lattice

def valid_move(lattice, i):
    return (lattice[i+1] == 0) and (lattice[i] > 0)
# rib_index = np.where(lattice>0)[0]
# print(rib_index)
# proteins_produced = 0
ind = 0
proteins_produced = np.zeros(len(tarr))
LM = np.empty([1,L])
for t in range(len(tarr)):
    if (kp_mrna * dt > rand())  :
        np.append(LM,np.zeros(L),0)
        print(LM)
    for lattice in LM:
        # task 3a, ribflow
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
        if (rand() < dt * kd_mrna * proteins_produced[ind]):
            proteins_produced[ind] -= 1
        # print(f"{lattice}\n")
        # print(proteins_produced)
        ind +=1
        # sleep(.1)
        if (kd_mrna * dt > rand()):
            Nr += np.sum(lattice)
            np.delete(LM,lattice)
    # print(proteins_produced)

# plt.plot(tarr, proteins_produced)
# plt.show()
