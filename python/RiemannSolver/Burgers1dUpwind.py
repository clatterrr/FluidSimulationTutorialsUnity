import numpy as np
#Upwind Scheme for 1d Burgers equation
#介绍：https://zhuanlan.zhihu.com/p/331771977
nmax = 8
tmax = 100
dx = 1.0 / nmax
dt = 1.0 / tmax
U = np.zeros((tmax,nmax))
F = np.zeros((tmax,nmax))
for i in range(0,nmax):
    if(i < nmax/2):
        U[0,i] = 0.1*i
    else:
        U[0,i] = 0.2*i
F[0,:] = 0.5 * U[0,:] * U[0,:]

for k in range(0,tmax-1):
    for i in range(1,nmax):
        U[k+1,i] = U[k,i] - (dt/dx)*(F[k,i] - F[k,i-1])
    U[k+1,0] = U[k+1,1]
    F[k+1,:] = 0.5 * U[k+1,:] * U[k+1,:]
    