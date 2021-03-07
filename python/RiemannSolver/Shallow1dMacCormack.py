import numpy as np
#MacCormack for 1d Shallow Water Equation
#介绍：https://zhuanlan.zhihu.com/p/331771977
nmax = 15
tmax = 100
h = np.ones((tmax,nmax+2))
for i in range(0,4):
    h[0,i] = 2
    
for i in range(13,17):
    h[0,i] = 2
u = np.zeros((tmax,nmax+2))
U1 = np.zeros((tmax,nmax+2))
U2 = np.zeros((tmax,nmax+2))
F1 = np.zeros((tmax,nmax+2))
F2 = np.zeros((tmax,nmax+2))

U1pred = np.zeros((tmax,nmax+2))
U2pred = np.zeros((tmax,nmax+2))
U1[0,:] = h[0,:]
U2[0,:] = F1[0,:] = h[0,:]*u[0,:]
F2[0,:] = h[0,:]*u[0,:]*u[0,:] + 0.5*10*h[0,:]*h[0,:]

for k in range(0,tmax-2):
    dx = 10
    cfl = 0.8
    dt = cfl*dx/max(abs(u[k,:]) + np.sqrt(abs(10*h[k,:])))
    c = dt / dx
    for i in range(1,nmax+1):
        U1pred[k,i] = U1[k,i] - c*(F1[k,i+1]-F1[k,i])
        U2pred[k,i] = U2[k,i] - c*(F2[k,i+1]-F2[k,i])
        h[k+1,i] = U1pred[k,i]
        u[k+1,i] = U2pred[k,i]/U1pred[k,i]
        
    h[k+1,0] = h[k+1,1]
    u[k+1,0] = u[k+1,1]
    h[k+1,-1] = h[k+1,-2]
    u[k+1,-1] = u[k+1,-2] 
    
    F1[k+1,:] = h[k+1,:]*u[k+1,:]
    F2[k+1,:] = h[k+1,:]*u[k+1,:]*u[k+1,:] + 0.5*10*h[k+1,:]*h[k+1,:]
    
    for i in range(1,nmax+1):
        U1[k+1,i] = 0.5 * (U1[k,i] + U1pred[k,i] - c * (F1[k+1,i] - F1[k+1,i-1]))
        U2[k+1,i] = 0.5 * (U2[k,i] + U2pred[k,i] - c * (F2[k+1,i] - F2[k+1,i-1]))
        h[k+1,i] = U1[k+1,i]
        u[k+1,i] = U2[k+1,i]/U1[k+1,i]
    
    h[k+1,0] = h[k+1,1]
    u[k+1,0] = u[k+1,1]
    h[k+1,-1] = h[k+1,-2]
    u[k+1,-1] = u[k+1,-2] 
    
    F1[k+1,:] = h[k+1,:]*u[k+1,:]
    F2[k+1,:] = h[k+1,:]*u[k+1,:]*u[k+1,:] + 0.5*10*h[k+1,:]*h[k+1,:]