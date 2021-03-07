import numpy as np
#lax-Friedchs Scheme for 1d Burgers equation
#介绍：https://zhuanlan.zhihu.com/p/331771977
nmax = 8
tmax = 100
dx = 1.0 / nmax
dt = 1.0 / tmax
U = np.zeros((tmax,nmax))
F = np.zeros((tmax,nmax))
f = np.zeros((tmax,nmax + 1))
for i in range(0,nmax):
    if(i < nmax/2):
        U[0,i] = 2
    else:
        U[0,i] = 1
F[0,:] = 0.5 * U[0,:] * U[0,:]

cfl = 0.8
dx = 0.1
dt = 0.01

for k in range(0,tmax - 1):
    for i in range(0,nmax+1):
        if(i == 0):
            f[k,i] = F[k,i]
        elif(i == nmax):
            f[k,i] = F[k,i-1]
        else:
            f[k,i] = 0.5*(F[k,i] + F[k,i-1]) - 0.5*cfl*dx/dt*(U[k,i] - U[k,i-1])
            
    for i in range(0,nmax):
        U[k+1,i] = U[k,i] - dt/dx*(f[k,i+1] - f[k,i])
    
    F[k+1,:] = 0.5 * U[k+1,:] * U[k+1,:] 