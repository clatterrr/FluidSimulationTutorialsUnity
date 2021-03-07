import numpy as np
#初始条件为一次函数h(x) = -x + 10，速度为常数的有限体积法
nmax = 10
tmax = 50
U = np.zeros((tmax,nmax))
Uhalf = np.zeros((tmax,nmax+1))
F = np.zeros((tmax,nmax+1))
a = 1#速度
dx = 1
for i in range(0,nmax):
    U[0,i] = 9.5 - i#piecewise 

for t in range(0,tmax-1):
    for i in range(0,nmax-1):
        Uhalf[t,i+1] = U[t,i] + (U[t,i+1] - U[t,i])/2/dx
        F[t,i+1] = a*Uhalf[t,i+1]
        
    #边界条件，外插法
    Uhalf[t,0] = U[t,0] - (U[t,1] - U[t,0])/2/dx
    F[t,0] = a*Uhalf[t,0]
    Uhalf[t,nmax] = U[t,nmax-1] + (U[t,nmax-1] - U[t,nmax-2])/2/dx
    F[t,nmax] = a*Uhalf[t,nmax]
    
    for i in range(0,nmax):
        U[t+1,i] = U[t,i] + (F[t,i] - F[t,i+1])