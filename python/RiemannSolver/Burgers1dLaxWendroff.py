import numpy as np
#lax-Wendroff Scheme for 1d Burgers equation
#版本1为详细版，版本2为简化版
#介绍：https://zhuanlan.zhihu.com/p/331771977
nmax = 8
tmax = 100
dx = 1.0 / nmax
dt = 1.0 / tmax
U = np.zeros((tmax,nmax))
F = np.zeros((tmax,nmax))
upredp = np.zeros((tmax,nmax))#plus 1/2
upredm = np.zeros((tmax,nmax))#minus 1/2
fpredp = np.zeros((tmax,nmax))#plus 1/2
fpredm = np.zeros((tmax,nmax))#minus 1/2
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
    for i in range(0,nmax-1):
        upredp[k,i] = 0.5 * (U[k,i+1] + U[k,i] - dt / dx * (F[k,i+1] - F[k,i]))
    for i in range(1,nmax):
        upredm[k,i] = 0.5 * (U[k,i] + U[k,i-1] - dt / dx * (F[k,i] - F[k,i-1]))
        
    upredp[k,-1] = upredp[k,-2]
    upredm[k,0] = upredm[k,1]
    fpredp[k,:] = 0.5 * upredp[k,:] * upredp[k,:]
    fpredm[k,:] = 0.5 * upredm[k,:] * upredm[k,:]
            
    for i in range(1,nmax-1):
        U[k+1,i] = U[k,i] - dt/dx*(fpredp[k,i] - fpredm[k,i])
    
    U[k+1,0] = U[k+1,1]
    U[k+1,-1] = U[k+1,-2]
    F[k+1,:] = 0.5 * U[k+1,:] * U[k+1,:]