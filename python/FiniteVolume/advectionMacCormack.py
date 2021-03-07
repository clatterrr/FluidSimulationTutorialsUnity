import numpy as np
#MacCormack算法。数值易爆炸
nmax = 510
tmax = 1005
U = np.zeros((tmax,nmax))
F = np.zeros((tmax,nmax))
upred = np.zeros((tmax,nmax))#plus 1/2

a = 0.5 #速度
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
    #U[0,i] = -i*i/16 + i #原来的真正的二阶函数
    
F[0,:] = a*U[0,:]
for t in range(0,tmax-1):
    for i in range(1,nmax-1):
        upred[t,i] = U[t,i] - dt / dx * (F[t,i+1] - F[t,i])
        
    upred[t,0] = upred[t,1]
    upred[t,-1] = upred[t,-2]#边界情况
    F[t+1,:] = a * upred[t,:]
            
    for i in range(1,nmax-1):
        U[t+1,i] = 0.5*(U[t,i] + upred[t,i]) - dt/dx*(F[t+1,i] - F[t+1,i-1])
    U[t+1,0] = U[t+1,1]
    U[t+1,-1] = U[t+1,-2]
    F[t+1,:] = a*U[0,:]
    
    