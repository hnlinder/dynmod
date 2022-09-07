import numpy as np
import matplotlib.pyplot as plt
# print("hej")

L  = 10
lattice  = np.zeros(L)

alpha = 1
beta = 1
Np0 = 0
tmax = 50000
kp = 1/60
kd = 1/1800

meanlength = int(tmax/2)
mean = np.empty(10)
tarr = np.linspace(0,tmax-1,tmax)
for i in range(0,10):
    Np0 = i*10

    Np = np.empty(tmax)
    Np[0] = Np0

    for t in range(1,tmax):
        randi1 = np.random.randint(0,600)

        if Np[t-1]!=0:
            randi2 = np.random.randint(0,1800/Np[t-1])
        else: randi2 = 1
        dNp =(randi1==0) -(randi2==0)
        Np[t] = Np[t-1]+dNp

    # plt.figure()
    plt.plot(tarr, Np)
    mean[i] = np.mean(Np[meanlength:])

    # plt.hold()
    # plt.show()
# plt.legend([])


print(mean)
print(np.mean(mean))
plt.show()
