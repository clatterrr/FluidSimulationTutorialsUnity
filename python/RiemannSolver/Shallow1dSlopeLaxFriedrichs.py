import numpy as np
nmax = 16
tmax = 100
h = np.ones((tmax,nmax))
for i in range(0,4):
    h[0,i] = 2

u = np.zeros((tmax,nmax))
U1 = np.zeros((tmax,nmax))
U2 = np.zeros((tmax,nmax))
F1 = np.zeros((tmax,nmax))
F2 = np.zeros((tmax,nmax))
z = np.zeros(nmax + 1)
H = np.zeros((tmax,nmax))

z[5] = 1
z[6] = 2
z[7] = 3
z[8] = 2
z[9] = 1
U1[0,:] = h[0,:]
U2[0,:] = F1[0,:] = h[0,:]*u[0,:]
F2[0,:] = h[0,:]*u[0,:]*u[0,:] + 0.5*10*h[0,:]*h[0,:]

f1 = np.zeros((tmax,nmax+1))
f2 = np.zeros((tmax,nmax+1))

cfl = 0.8
dx = 0.25

for k in range(0,tmax-1):
    dt = cfl*dx/max(abs(u[k,:]) + np.sqrt(abs(10*h[k,:])))
    for i in range(0,nmax+1):
        if(i == 0):
            f1[k,i] = F1[k,i]
            f2[k,i] = F2[k,i]
        elif(i == nmax):
            f1[k,i] = F1[k,i-1]
            f2[k,i] = F2[k,i-1]
        else:
            f1[k,i] = 0.5*(F1[k,i]+F1[k,i-1]) - 0.5*cfl*dx/dt * (U1[k,i] - U1[k,i-1])
            f2[k,i] = 0.5*(F2[k,i]+F2[k,i-1]) - 0.5*cfl*dx/dt * (U2[k,i] - U2[k,i-1])
            
    for i in range(0,nmax):
        U1[k+1,i] = U1[k,i] - dt/dx * (f1[k,i+1] - f1[k,i])
        U2[k+1,i] = U2[k,i] - dt/dx * (f2[k,i+1] - f2[k,i]) - dt * 10 * h[k,i] * (z[i+1] - z[i])
    
    h[k+1,:] = U1[k+1,:]
    u[k+1,:] = U2[k+1,:]/U1[k+1,:]
    F1[k+1,:] = U2[k+1,:]
    F2[k+1,:] = h[k+1,:]*u[k+1,:]*u[k+1,:] + 0.5*10*h[k+1,:]*h[k+1,:]
    
    for i in range(0,nmax):
        H[k,i] = h[k,i] + z[i]