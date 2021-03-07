import numpy as np
#LaxFriedrichs算法，二阶精度。只有当a = 1时才完全准确。其它时候数值耗散很大
nmax = 510
tmax = 1005
U = np.zeros((tmax,nmax))
f = np.zeros((tmax,nmax+1))
F = np.zeros((tmax,nmax))

a = 0.5#速度
b = np.zeros((2,nmax+1))
dt = 1
dx = 1

for i in range(0,nmax+1):
    b[0,i] = -i*i/2 + i*10
    b[1,i] = -i*i*i*i/6 + i*i*5
    j = i + 0.5
    #b[1,i] = -j*j*j/6 + j*j*5
for i in range(0,nmax):
    if i <= 8:
        U[0,i] = -i*i/16 + i
    elif i >= 16:
        U[0,i] = -(i-8)*(i-8)/16 + (i-8)
    else:
        U[0,i] = 4
    # U[0,i] = -i*i*i/16 + i #原来的真正的二阶函数
    
F[0,:] = a*U[0,:]
for t in range(0,tmax-1):
    for i in range(0,nmax+1):
        if(i == 0):
            f[t,i] = F[t,i]
        elif(i == nmax):
            f[t,i] = F[t,i-1]
        else:
            f[t,i] = 0.5*(F[t,i] + F[t,i-1]) - 0.5*dx/dt*(U[t,i] - U[t,i-1])
    for i in range(0,nmax):
        U[t+1,i] = U[t,i] - dt/dx*(f[t,i+1] - f[t,i])
    F[t+1,:] = a*U[t+1,:]
    
    