import numpy as np
import matplotlib.pyplot as plt
from numpy.random import rand
from time import sleep
# print("hej")

L  = 30

alpha = 1
beta = .5
q = 1

tmax = 600
dt = 10

kp = 1/60
kd = 1/1800
Np0 = 0
lattice  = np.zeros(L)

meanlength = int(tmax/2)
mean = np.empty(10)
tarr = np.linspace(0, tmax-1, int((tmax-1)/dt))
# print(tarr)


def move_rib(lattice, i):
    lattice[i] -=1
    lattice[i+1] +=1
    return lattice

def valid_move(lattice, i):
    return (lattice[i+1] == 0) and (lattice[i] > 0)
# rib_index = np.where(lattice>0)[0]
# print(rib_index)
proteins_produced = 0
alpha_array = np.linspace(.2, 1, 5)
beta_array = [.25, .5]
occupied = np.zeros(alpha_array * beta_array)
for alpha in alpha_array:
    for beta in beta_array:
        for t in range(len(tarr)):
            rib_index = np.where(lattice>0)[0]
            # print(rib_index)
            np.random.shuffle(rib_index)
            for i in rib_index:
                if i != L-1:
                    if ((rand() < dt * q) and valid_move(lattice,i)):
                        lattice[i+1] += 1
                        lattice[i] -= 1


                    # lattice[i+1] += ((rand() < q) and (lattice[i+1] == 0) and (lattice[i] > 0))
                    # lattice[i] -= ((rand() < q) and (lattice[i+1] == 0) and (lattice[i] > 0))

            lattice[0] += ((lattice[0]==0) and (rand() < dt * alpha))
            if ((lattice[L-1]>0)):
                occupied +=1
                if (rand() < dt * beta):
                    lattice[L-1] -= 1
                    proteins_produced +=1
            # print(f"{lattice}\n")
            # sleep(.1)
        print(proteins_produced)

        print(f"Occupational probability : {occupied/len(tarr)}")

