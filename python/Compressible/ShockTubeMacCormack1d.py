import numpy as np
#Lax-Friedriches For 1d Shock Tube
tmax = 1000
nmax = 12



rho = np.zeros((tmax,nmax))#密度
u = np.zeros((tmax,nmax))#速度
E = np.zeros((tmax,nmax))#能量
p = np.zeros((tmax,nmax))#压力
e = np.zeros((tmax,nmax))
a = np.zeros((tmax,nmax))#声速
H = np.zeros((tmax,nmax))#焓变

U1 = np.zeros((tmax,nmax))#rho
U2 = np.zeros((tmax,nmax))#rho*u
U3 = np.zeros((tmax,nmax))#rho*E

F1 = np.zeros((tmax,nmax))#rho
F2 = np.zeros((tmax,nmax))#rho*u
F3 = np.zeros((tmax,nmax))#rho*E

u1pred = np.zeros((tmax,nmax))
u2pred = np.zeros((tmax,nmax))
u3pred = np.zeros((tmax,nmax))

gamma = 1.4
dt = 1
dx = 0.2
cfl = 0.9

for i in range(0,nmax):
    if(i < nmax/2):
        rho[0,i] = 1
        p[0,i] = 1
    else:
        rho[0,i] = 0.125
        p[0,i] = 0.1

for i in range(0,nmax):
    U1[0,i] = rho[0,i]
    U2[0,i] = rho[0,i]*u[0,i]
    #E = p ./ ((gamma-1) * rho) + u.^2/2;
    E[0,i] = p[0,i]/((gamma-1)*rho[0,i]) + u[0,i]*u[0,i]/2
    U3[0,i] = rho[0,i]*E[0,i]
    F1[0,i] = rho[0,i]*u[0,i]
    F2[0,i] = rho[0,i]*u[0,i]*u[0,i] + p[0,i]
    F3[0,i] = rho[0,i]*u[0,i]*H[0,i]

for t in range(0,tmax-1):
    
    smax = max(abs(u[t,:] + np.sqrt(gamma * p[t,:]/rho[t,:])))
    dt = dx * cfl/smax
    
    for i in range(0,nmax-1):
        u1pred[t,i] = U1[t,i] - dt / dx * (F1[t,i+1] - F1[t,i])
        u2pred[t,i] = U2[t,i] - dt / dx * (F2[t,i+1] - F2[t,i])
        u3pred[t,i] = U3[t,i] - dt / dx * (F3[t,i+1] - F3[t,i])
    
    for i in range(0,nmax-1):
        rho[t+1,i] = u1pred[t,i]
        u[t+1,i] = u2pred[t,i]/rho[t+1,i]
        E[t+1,i] = u3pred[t,i]/rho[t+1,i]
        p[t+1,i] = (gamma-1)*(E[t+1,i] - u[t+1,i]*u[t+1,i]/2)*rho[t+1,i]
        e[t+1,i] = p[t+1,i]/(gamma-1)/rho[t+1,i]
        a[t+1,i] = np.sqrt(gamma * p[t+1,i]/rho[t+1,i])
        H[t+1,i] = (p[t+1,i] * gamma)/((gamma-1) * rho[t+1,i]) + u[t+1,i]*u[t+1,i]/2
        
        F1[t+1,i] = rho[t+1,i]*u[t+1,i]
        F2[t+1,i] = rho[t+1,i]*u[t+1,i]*u[t+1,i]+p[t+1,i]
        F3[t+1,i] = rho[t+1,i]*u[t+1,i]*H[t+1,i]
    
    
    for i in range(1,nmax-1):
        U1[t+1,i] = 0.5*(U1[t,i] + u1pred[t,i] - dt/dx * (F1[t+1,i] - F1[t+1,i-1]))
        U2[t+1,i] = 0.5*(U2[t,i] + u2pred[t,i] - dt/dx * (F2[t+1,i] - F2[t+1,i-1]))
        U3[t+1,i] = 0.5*(U3[t,i] + u3pred[t,i] - dt/dx * (F3[t+1,i] - F3[t+1,i-1]))
        
    U1[t+1,0] = U1[t+1,1]
    U1[t+1,-1] = U1[t+1,-2]
    U2[t+1,0] = U2[t+1,1]
    U2[t+1,-1] = U2[t+1,-2]
    U3[t+1,0] = U3[t+1,1]
    U3[t+1,-1] = U3[t+1,-2]
    
    
    for i in range(0,nmax):
        rho[t+1,i] = U1[t+1,i]
        u[t+1,i] = U2[t+1,i]/rho[t+1,i]
        E[t+1,i] = U3[t+1,i]/rho[t+1,i]
        p[t+1,i] = (gamma-1)*(E[t+1,i] - u[t+1,i]*u[t+1,i]/2)*rho[t+1,i]
        e[t+1,i] = p[t+1,i]/(gamma-1)/rho[t+1,i]
        a[t+1,i] = np.sqrt(gamma * p[t+1,i]/rho[t+1,i])
        H[t+1,i] = (p[t+1,i] * gamma)/((gamma-1) * rho[t+1,i]) + u[t+1,i]*u[t+1,i]/2
        
        F1[t+1,i] = rho[t+1,i]*u[t+1,i]
        F2[t+1,i] = rho[t+1,i]*u[t+1,i]*u[t+1,i]+p[t+1,i]
        F3[t+1,i] = rho[t+1,i]*u[t+1,i]*H[t+1,i]
    