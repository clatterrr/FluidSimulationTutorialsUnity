import numpy as np
nmax = 7
tmax = 1000
v = np.zeros((nmax-1))
v2 = np.zeros((nmax-1))
div = np.zeros((nmax))
div2 = np.zeros((nmax))
dx = 1
for i in range(1,nmax-2):
    v[i] = i#这步假设了一个初始有散度的速度场
for i in range(1,nmax-1):
    div[i] = (v[i]-v[i-1])/dx
p = np.zeros((tmax,nmax))
left = np.ones((nmax))
right = np.ones((nmax))
left[4] = 2
right[3] = 2
for t in range(0,tmax-2):
    for i in range(1,nmax-1):
        term = (left[i]*p[t,i-1] + right[i]*p[t,i+1] - dx*dx*div[i])/(left[i] + right[i])
        p[t+1,i] = p[t,i] + (term - p[t,i])
for i in range(0,nmax-1):
    v2[i] = v[i] - right[i]*(p[t,i+1] - p[t,i])
for i in range(1,nmax-1):#重新计算散度，应当接近零
    div2[i] = v2[i]-v2[i-1]
    div2[i] = div[i] - (p[t,i-1] - p[t,i])*left[i]*dx - (p[t,i+1] - p[t,i])*right[i]*dx