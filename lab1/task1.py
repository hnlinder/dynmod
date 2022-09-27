import numpy as np
import matplotlib.pyplot as plt
# print("hej")

L  = 10
lattice  = np.zeros(L)

alpha = 1
beta = 1
Np0 = 0
tmax = 1800 * 100
kp = 1/60
kd = 1/1800

# proteins_per_mRNA = 1
nr_runs = 1
meanlength = int(tmax/2)
mean = np.empty([nr_runs, 2])
tarr = np.linspace(0,tmax-1,tmax)
nparr = np.empty([nr_runs, 2, tmax])
twomean = np.zeros([2, tmax])
Np0 = 0
Np = np.empty([nr_runs, 2, tmax])
for n in range(nr_runs):

    for i, proteins_per_mRNA in enumerate([1,10]):

        Np[0,0,0] = Np0

        for t in range(1,tmax):
            randi1 = np.random.randint(0,600)

            if Np[n,i,t-1]!=0:
                randi2 = np.random.randint(0,1800/Np[n,i,t-1])
            else: randi2 = 1
            dNp =proteins_per_mRNA*(randi1==0) -(randi2==0)
            Np[n, i, t] = Np[n,i,t-1]+dNp

        # plt.figure()
        # plt.plot(tarr, Np)
        twomean[i, :] += Np[n, i, :]
        mean[n, i] = np.mean(Np[n,i,meanlength:])

    # plt.hold()
    # plt.show()
# plt.legend([])

twomean = twomean/nr_runs

print(f"mean3 {(Np[0,0,meanlength:])}")
# print(np.mean(mean))

fano3= np.var(Np[0,0,meanlength:])/np.mean(Np[0,0,meanlength:])
fano30= np.var(Np[0,1,meanlength:])/np.mean(Np[0,1,meanlength:])

print(f"Fano factor 3: {fano3}\nFano factor 30: {fano30}")
# np.mean()

# mean1 = np.mean(Np[:,1,:], axis=0)
# print(f"meanmean : {np.mean(mean1)}")
plt.plot(tarr, Np[0,0,:] ,label="1 protein per mRNA")
plt.plot(tarr,twomean[1,:], "--", label="10 proteins per mRNA")
plt.xlabel("Time [s]")
plt.ylabel("Proteins produced [1]")
# plt.plot(tarr, np.ones(tmax)*3,"-.", label =  r"$N_p = 3$")
# plt.plot(tarr, 30*np.ones(tmax),".", label= r"$N_p = 30$")
plt.hlines(3,tarr[0], tarr[-1], "b", label= r"$N_p = 3$")
plt.hlines(30,tarr[0], tarr[-1], "r", label= r"$N_p = 30$")
plt.legend()
plt.show()
