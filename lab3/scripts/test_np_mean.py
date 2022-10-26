import numpy as np



len1 = 3
len2 = 5
myarr = np.zeros([len1,len2])

for i in range(len1):
    for j in range(len2):
        myarr[i,j] = j
print(myarr)
ax0 = np.mean(myarr,axis=0)
ax1 = np.mean(myarr,axis=1)
print(ax0)
print(ax1)

