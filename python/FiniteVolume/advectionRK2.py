import numpy as np
#初始条件为二次函数，使用二阶Runge-Kutta的有限体积法
nmax = 20
tmax = 50
U = np.zeros((tmax,nmax))
U2 = np.zeros((tmax,nmax))
Uhalf = np.zeros((tmax,nmax+1))
F = np.zeros((tmax,nmax+1))
a = 1#速度
dx = 1
b = np.zeros((2,nmax+1))
for i in range(0,nmax+1):
    b[0,i] = -i*i/2 + i*10
    b[1,i] = -i*i*i/6 + i*i*5
for i in range(0,nmax):
    U[0,i] = b[1,i+1] - b[1,i]
    
def rk2():
    for i in range(0,nmax-1):
        Uhalf[t,i+1] = U[t,i] + (U[t,i+1] - U[t,i])/2/dx
        F[t,i+1] = a*Uhalf[t,i+1]    
    #边界条件，外插法
    Uhalf[t,0] = U[t,0] - (U[t,1] - U[t,0])/2/dx
    F[t,0] = a*Uhalf[t,0]
    Uhalf[t,nmax] = U[t,nmax-1] + (U[t,nmax-1] - U[t,nmax-2])/2/dx
    F[t,nmax] = a*Uhalf[t,nmax]
    

for t in range(0,tmax-1):

    U2[t,:] = U[t,:]
    rk2()
    for i in range(0,nmax):
        U[t,i] = U[t,i] + (F[t,i] - F[t,i+1])
    
    rk2()
    for i in range(0,nmax):
        U[t+1,i] = (U2[t,i] + U[t,i] + (F[t,i] - F[t,i+1]))/2
    U[t,:] = U2[t,:]