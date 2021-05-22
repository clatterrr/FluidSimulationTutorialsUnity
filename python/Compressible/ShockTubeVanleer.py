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

f1 = np.zeros((tmax,nmax+1))#rho
f2 = np.zeros((tmax,nmax+1))#rho*u
f3 = np.zeros((tmax,nmax+1))#rho*E

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

for t in range(0,tmax-2):
    
    
    for i in range(0,nmax):
        rho[t,i] = U1[t,i]
        u[t,i] = U2[t,i]/rho[t,i]
        E[t,i] = U3[t,i]/rho[t,i]
        p[t,i] = (gamma-1)*(E[t,i] - u[t,i]*u[t,i]/2)*rho[t,i]
        e[t,i] = p[t,i]/(gamma-1)/rho[t,i]
        a[t,i] = np.sqrt(gamma * p[t,i]/rho[t,i])
        H[t,i] = (p[t,i] * gamma)/((gamma-1) * rho[t,i]) + u[t,i]*u[t,i]/2
        
        F1[t,i] = rho[t,i]*u[t,i]
        F2[t,i] = rho[t,i]*u[t,i]*u[t,i]+p[t,i]
        F3[t,i] = rho[t,i]*u[t,i]*H[t,i]
        
    smax = max(abs(u[t,:] + np.sqrt(gamma * p[t,:]/rho[t,:])))
    dt = dx * cfl/smax