import numpy as np
import matplotlib.pyplot as plt
from numpy.random import rand
from time import sleep
# print("hej")

L  = 30

alpha = .1
beta = .5
q = 1

tmax = 60 * 1000
dt = 1

kp = 1/600
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
alpha_array = np.linspace(0, 1,20)

beta_array = [.25, .5]

# J = np.zeros([len(alpha_array), len(beta_array)])

proteins_produced = np.zeros([len(alpha_array), len(beta_array), len(tarr)])
occupied = np.zeros([len(alpha_array) ,len(beta_array)])
for a, alpha in enumerate(alpha_array):
    for b, beta in enumerate(beta_array):
        lattice  = np.zeros(L)
        for t in range(len(tarr)):
            proteins_produced[a,b,t] = proteins_produced[a,b,t-1]

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
                if tarr[t] > meanlength:
                    occupied[a,b] +=1
                if (rand() < dt * beta):
                    lattice[L-1] -= 1
                    proteins_produced[a, b, t] +=1
                    # print(proteins_produced[a,b,t])
            # print(f"{lattice}\n")
            # sleep(.1)
            if (rand() < kd * proteins_produced[a,b,t]):
                proteins_produced[a,b,t] -= 1

        
        if (a == 0 and b == 0) or (a == 0 and b ==1) or (alpha == 1 and beta == .25) or (alpha == 1 and beta == .5):
            plt.figure(10)
            plt.plot(tarr, proteins_produced[a, b, :], label=fr"$\alpha = {alpha}, \beta = {beta}$")
        # print(f"Occupational probability : {occupied/len(tarr)}")

def plot_mean_field():
    alpha_under_half = np.linspace(0,.5,100)
    # print(alpha_under_half)
    alpha_over_half = np.linspace(.5,1,100)
    y1 =  alpha_under_half*(1-alpha_under_half)
    # print(y)
    plt.plot(alpha_under_half,y1, "r","--")
    plt.hlines(.25*(1-.25),.2445, 1,"k","--" ,label="Mean-field solution, "+r"$\alpha\leq 0.5,  \beta=0.25$")
    plt.hlines(.5*(1-.5),.5, 1, "r","-.",label="Mean-field solution, "+r"$\alpha>0.5, \beta=0.5$")
# print(occupied/len(tarr)*[1, 0])
plt.figure(10)
plt.legend(loc="best")
plt.xlabel("Time [s]")
plt.ylabel("Proteins produced [1]")

J = occupied/len(tarr)*beta_array
# print(J)
plt.figure(1)
plt.plot(alpha_array, J[:,0], label=r"$\beta = 0.25$")
plt.plot(alpha_array, J[:,1], label=r"$\beta = 0.5$")
plt.xlabel(r"$\alpha [s^{-1}]$")
plt.ylabel(r"$J [s^{-1}]$")
plot_mean_field()
plt.legend(loc="best")
plt.show()
