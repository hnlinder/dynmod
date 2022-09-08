import numpy as np
import matplotlib.pyplot as plt
from random import random as rand

# Make sure you choose a time step so that all probabilities aren't >= 1. In
# the first task, the rates are 1/600 and 1/1800 so we can use a big timestep to
# reduce times where nothing happens.
t_end = 10*3600  # end time of simulation (ten hours)
t_0 = 0  # start time
dt = 10  # time step

# create a list of all steps in time [0, dt, 2*dt, ..., t_end-dt, t_end]
times = range(t_0, t_end, dt)

# Define the rates
kp = 1/600
kd = 1/1800

# We store the number of proteins each time, starting at 0
N_prot = np.zeros(len(times))

# Loop through all times (skip the first one since that's the initial)
for i in range(1, len(times)):
    # The old protein amount carries over
    N_prot[i] = N_prot[i-1]

    # For each timestep, the probability of an event with rate k to happen is
    # k*dt. To check this each timestep, we check if a random uniform number
    # between [0, 1) is smaller than this probability.

    if rand() < kp*dt:  # a protein is created
        N_prot[i] += 1

    if rand() < kd*dt*N_prot[i-1]:  # a protein decays
        N_prot[i] -= 1

# Plot the results
plt.plot(times, N_prot)
plt.show()
