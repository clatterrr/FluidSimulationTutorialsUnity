import numpy as np
#Lax Friedrichs for 2d Shallow Water Equation
#介绍：https://zhuanlan.zhihu.com/p/331771977
nmax = 8
tmax = 100
h = np.ones((nmax,nmax))
for i in range(0,nmax):
    h[i,:] = 1
for i in range(0,4):
    h[i,:] = 2

u = np.zeros((nmax,nmax))
v = np.zeros((nmax,nmax))
U1 = np.zeros((nmax,nmax))
U2 = np.zeros((nmax,nmax))
U3 = np.zeros((nmax,nmax))
F1 = np.zeros((nmax,nmax))
F2 = np.zeros((nmax,nmax))
F3 = np.zeros((nmax,nmax))
G1 = np.zeros((nmax,nmax))
G2 = np.zeros((nmax,nmax))
G3 = np.zeros((nmax,nmax))

f1 = np.zeros((nmax+1,nmax))
f2 = np.zeros((nmax+1,nmax))
f3 = np.zeros((nmax+1,nmax))
g1 = np.zeros((nmax,nmax+1))
g2 = np.zeros((nmax,nmax+1))
g3 = np.zeros((nmax,nmax+1))

dx = dy = 2.0/nmax
cfl = 0.8

for i in range(0,nmax):
    for j in range(0,nmax):
        U1[i,j] = h[i,j]
        U2[i,j] = h[i,j] * u[i,j]
        U3[i,j] = h[i,j] * v[i,j]
        F1[i,j] = h[i,j] * u[i,j]
        F2[i,j] = h[i,j] * u[i,j] * u[i,j] + 0.5 * 10 * h[i,j] * h[i,j]
        F3[i,j] = h[i,j] * u[i,j] * v[i,j]
        G1[i,j] = h[i,j] * v[i,j]
        G2[i,j] = h[i,j] * u[i,j] * v[i,j]
        G3[i,j] = h[i,j] * v[i,j] * v[i,j] + 0.5 * 10 * h[i,j] * h[i,j]



for k in range(0,tmax-1):
    dxx = dx/np.max(abs(u)+np.sqrt(10*h))
    dyy = dy/np.max(abs(v)+np.sqrt(10*h))
    dt = 0.5 * cfl * min(dxx,dyy)
    dx = 0.1
    dt = 0.01
    for i in range(0,nmax+1):
        for j in range(0,nmax):
            if(i == 0):
                f1[i,j] = F1[i,j]
                f2[i,j] = F2[i,j]
                f3[i,j] = F3[i,j]
            elif(i == nmax):
                f1[i,j] = F1[i-1,j]
                f2[i,j] = F2[i-1,j]
                f3[i,j] = F3[i-1,j]
            else:
                f1[i,j] = 0.5*(F1[i,j]+F1[i-1,j]) - dx/dt*0.25*cfl*(U1[i,j]-U1[i-1,j])
                f2[i,j] = 0.5*(F2[i,j]+F2[i-1,j]) - dx/dt*0.25*cfl*(U2[i,j]-U2[i-1,j])
                f3[i,j] = 0.5*(F3[i,j]+F3[i-1,j]) - dx/dt*0.25*cfl*(U3[i,j]-U3[i-1,j])
    
    t = 1    
    for j in range(0,nmax+1):
        for i in range(0,nmax):
            if(j == 0):
                g1[i,j] = G1[i,j]
                g2[i,j] = G2[i,j]
                g3[i,j] = G3[i,j]
            elif(j == nmax):
                g1[i,j] = G1[i,j-1]
                g2[i,j] = G2[i,j-1]
                g3[i,j] = G3[i,j-1]
            else:
                g1[i,j] = 0.5*(G1[i,j]+G1[i,j-1]) - dy/dt*0.25*cfl*(U1[i,j]-U1[i,j-1])
                g2[i,j] = 0.5*(G2[i,j]+G2[i,j-1]) - dy/dt*0.25*cfl*(U2[i,j]-U2[i,j-1])
                g3[i,j] = 0.5*(G3[i,j]+G3[i,j-1]) - dy/dt*0.25*cfl*(U3[i,j]-U3[i,j-1])
                
    t = 1    
    for i in range(0,nmax):
        for j in range(0,nmax):
            U1[i,j] = U1[i,j] - dt/dx*(f1[i+1,j]-f1[i,j]) - dt/dy*(g1[i,j+1]-g1[i,j])
            U2[i,j] = U2[i,j] - dt/dx*(f2[i+1,j]-f2[i,j]) - dt/dy*(g2[i,j+1]-g2[i,j])
            U3[i,j] = U3[i,j] - dt/dx*(f3[i+1,j]-f3[i,j]) - dt/dy*(g3[i,j+1]-g3[i,j])
    
    t = 1        
    h = U1.copy()
    u = U2.copy() / U1.copy()
    v = U3.copy() / U1.copy()
    for i in range(0,nmax):
        for j in range(0,nmax):
           F1[i,j] = h[i,j] * u[i,j]
           F2[i,j] = h[i,j] * u[i,j] * u[i,j] + 0.5 * 10 * h[i,j] * h[i,j]
           F3[i,j] = h[i,j] * u[i,j] * v[i,j]
           G1[i,j] = h[i,j] * v[i,j]
           G2[i,j] = h[i,j] * u[i,j] * v[i,j]
           G3[i,j] = h[i,j] * v[i,j] * v[i,j] + 0.5 * 10 * h[i,j] * h[i,j]




