import numpy as np
#hll scheme for 1d Burgers equation
#https://zhuanlan.zhihu.com/p/331771977
nmax = 8
tmax = 100
dx = 1.0 / nmax
dt = 1.0 / tmax
U = np.zeros((tmax,nmax))
F = np.zeros((tmax,nmax))
f = np.zeros((tmax,nmax + 1))
for i in range(0,nmax):
    if(i < 2):
        U[0,i] = 2
    else:
        U[0,i] = 1
F[0,:] = 0.5 * U[0,:] * U[0,:]
cfl = 0.8
dx = 0.1
dt = 0.01

for k in range(0,tmax-1):
    for i in range(0,nmax-1):
        sl = U[k,i]
        sr = U[k,i+1]
        if(sl >= 0):
            f[k,i+1] = F[k,i]
        elif(sr <= 0):
            f[k,i+1] = F[k,i+1]
        else:
            f[k,i+1] = (sr*F[k,i] - sl*F[k,i+1] + sr*sl*(U[k,i+1] - U[k,i]))/(sr - sl)
    f[k,0] = F[k,1]
    f[k,-1] = F[k,-1]
    for i in range(0,nmax):
        U[k+1,i] = U[k,i] - dt / dx * (f[k,i+1] - f[k,i])
    F[k+1,:] = 0.5 * U[k+1,:] * U[k+1,:]