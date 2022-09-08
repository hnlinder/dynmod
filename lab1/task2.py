import numpy as np
import matplotlib.pyplot as plt
# print("hej")

L  = 10
lattice  = np.zeros(L)

alpha = .1
beta = .1
q = .1
Np0 = 0
tmax = 10
kp = 1/60
kd = 1/1800

meanlength = int(tmax/2)
mean = np.empty(10)
tarr = np.linspace(0,tmax-1,tmax)



for t in range(tmax):
    for i,nr_ribosomes in enumerate(lattice):
        rand1 = np.random.rand()
        rand2 = np.random.rand()
        rand3 = np.random.rand()
        if i == 0:
            lattice[i] += (rand1 < alpha)
        elif i == L-1:
            lattice[i] += (rand1 < beta)
        else:
            print(i)
            lattice[i+1] += (rand1 < q)

    while 2 not in lattice:
        for i,nr_ribosomes in enumerate(lattice):
            if i == 0:
                if nr_ribosomes >1:
                    lattice[i] -= 1
            else:
                if nr_ribosomes > 1:
                    lattice[i-1] += 1
        print(lattice)

    print(lattice)


# print(mean)
# print(np.mean(mean))
# plt.plot(tarr, Np)
# plt.show()
