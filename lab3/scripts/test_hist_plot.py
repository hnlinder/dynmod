# import matplotlib as mpl
# import matplotlib.pyplot as plt
# import numpy as np

# x = np.random.random(20)
# y = np.random.random(20)
# z = np.random.random(20)

# fig = plt.figure()
# ax = fig.add_subplot(111)

# ax.hist(x, bins=np.arange(0, 1, 0.1), edgecolor='None', alpha = 0.5, color= 'b')
# ax.hist(y, bins=np.arange(0, 1, 0.1), edgecolor='None', alpha = 0.5, color= 'r')
# ax.hist(z, bins=np.arange(0, 1, 0.1), edgecolor="None", alpha = 0.5, color= 'k')


# ax.hist(x, bins=np.arange(0, 1, 0.1), ls='dashed', lw=3, facecolor="None")
# ax.hist(y, bins=np.arange(0, 1, 0.1), ls='dotted', lw=3, facecolor="None")
# ax.hist(z, bins=np.arange(0, 1, 0.1), lw=3, facecolor="None")

# plt.show()


# ORIGINAL Q
import matplotlib.pyplot as plt
import numpy as np

x = np.random.random(20)
y =np.random.random(20)
z= np.random.random(20)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(x, bins=np.arange(0, 1, 0.1), ls='dashed', edgecolor="k", alpha = 0.5, lw=3, facecolor= 'None')
ax.hist(y, bins=np.arange(0, 1, 0.1), ls='dotted', alpha = 0.5, lw=3, color= 'r')
ax.hist(z, bins=np.arange(0, 1, 0.1), alpha = 0.5, lw=3, color= 'k')
ax.set_xlim(-0.5, 1.5)
ax.set_ylim(0, 7)
plt.show()
