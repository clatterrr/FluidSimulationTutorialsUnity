import numpy as np
nmax = 7
tmax = 1000
v = np.zeros((nmax-1))
v2 = np.zeros((nmax-1))
div = np.zeros((nmax))
div2 = np.zeros((nmax))
dx = 0.1
for i in range(1,nmax-2):
    v[i] = i+1#这步假设了一个初始有散度的速度场
for i in range(1,nmax-1):
    div[i] = (v[i]-v[i-1])/dx
p = np.zeros((tmax,nmax))
for t in range(0,tmax-2):
    for i in range(1,nmax-1):
        term = (p[t,i-1] + p[t,i+1] - dx*dx*div[i])/2
        p[t+1,i] = p[t,i] + (term - p[t,i])
for i in range(0,nmax-1):
    v2[i] = v[i] + (p[t,i] - p[t,i+1])/dx
for i in range(1,nmax-1):#重新计算散度，应当是零
    div2[i] = (v2[i]-v2[i-1])/dx