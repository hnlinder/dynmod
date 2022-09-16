import numpy as np
import matplotlib.pyplot as plt
from numpy.random import rand
from time import sleep


L  = 30
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



        # print(lattice)
        # sleep(.2)

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
# class mRNA:
    # def __init__(self, lattice = np.zeros(L), ribs = 0):
        # self.lattice = lattice
        # self.ribs = ribs

    # # one time step of the ribosomes
    # def rib_step(self, proteins_produced, ind, alpha, beta, dt, lattice, Nr):
        # proteins_produced[ind] = proteins_produced[ind - 1]
        # rib_index = np.where(lattice>0)[0]
        # # print(rib_index)
        # np.random.shuffle(rib_index)
        # for i in rib_index:
            # if i != L-1:
                # if ((rand() < dt * q) and valid_move(lattice,i)):
                    # lattice[i+1] += 1
                    # lattice[i] -= 1
        # if ((lattice[0]==0) and (rand() < dt * alpha * Nr)):
            # lattice[0] += 1
            # Nr -=1
        # if ((lattice[L-1]>0) and (rand() < dt * beta)):
            # lattice[L-1] -= 1
            # proteins_produced[ind] += 1
            # Nr +=1
        # if (rand() < dt * kd * proteins_produced[ind]):
            # proteins_produced[ind] -= 1

        # return proteins_produced, lattice, Nr





# class LM:
    # def __init__(self, mRNA, matr):
        # self.mRNA = mRNA
        # self.mrna_matrix = matr

    # def create_mRNA():
        # np.append(lm, mRNA)

    # def remove_mRNa(index):
        # np.delete(mRNA)
        # return ribs


# lat1 =  mRNA(np.zeros(3))
# lm1 = LM(lat1)


# print(lat1.lattice)
# print(lm1.mRNA.lattice)


alpha = .9
beta = 2
q = .3



Nr = 14
tmax = 3600 
dt = 1/2

kp = 1/60
kd = 1/1800
Np0 = 0
kd_rib = 1/1800
kd_mrna = 1/300
kp_mrna = 1/600
# lattice  = np.zeros(L)

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


# def create_mRNA():
    # lat = 

ind = 0
proteins_produced = np.zeros(len(tarr))
# lat = mRNA(np.zeros(L))
lm = LM()
# print(lat.lattice)
# print(proteins_produced, ind, alpha, beta, dt, lat.lattice, Nr)
for t in range(len(tarr)):
    # make new mrna
    if (kp_mrna * dt > rand())  :
        # print("NEW RNA")
        lm.create_mRNA()
        # print(f"nr of RNAs : {len(lm.matr)}")
    # Do a time step with ribosomes for each mRNA
    for index, mrna in enumerate(lm.matr):
        proteins_produced, mrna.lattice, Nr = mrna.rib_step(proteins_produced, ind, alpha, beta, dt, mrna.lattice, Nr)
        # remove mrna
        if (kd_mrna * dt > rand()):
            # print("REMOVING RNA")
            # print(f"RNAs before remove : {len(lm.matr)}")
            Nr += lm.remove_mRNA(index)
            # print(f"RNAs after remove : {len(lm.matr)}")

    # if proteins_produced[ind] ==0:
        # try:
            # print(lm.matr[0].lattice)
            # sleep(.5)
        # except KeyboardInterrupt: 
            # pass
        # else: pass


    # print(f"{lattice}\n")
    # print(proteins_produced[ind])
    ind +=1
    # sleep(.1)
print(proteins_produced)

plt.plot(tarr, proteins_produced)
plt.show()
