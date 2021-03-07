import numpy as np
#Muscl minmod 算术平均通量
nmax = 50
tmax = 100
U = np.zeros((tmax,nmax))
U2 = np.zeros((tmax,nmax))
U3 = np.zeros((tmax,nmax))
ul = np.zeros((tmax,nmax+1))
ur = np.zeros((tmax,nmax+1))
F = np.zeros((tmax,nmax+1))
k = np.zeros((tmax,nmax))
phi = np.zeros((tmax,nmax+1))
a = 1#速度
dx = 1
b = np.zeros((2,nmax+1))

for i in range(0,nmax+1):
    b[0,i] = -i*i/2 + i*10
    b[1,i] = -i*i*i/6 + i*i*5
    j = i + 0.5
    #b[1,i] = -j*j*j/6 + j*j*5
for i in range(0,nmax):
    if i <= 8:
        U[0,i] = -i*i/16 + i
    elif i >= 16:
        U[0,i] = -(i-8)*(i-8)/16 + (i-8)
    else:
        U[0,i] = 4
    #U[0,i] = b[1,i+1] - b[1,i]
    

def rk2():
    for i in range(1,nmax-1):
        k[t,i] = (U[t,i+1] - U[t,i-1])/2#二阶精度的中心差分即可
    k[t,0] = k[t,1] + (k[t,1] - k[t,2])
    k[t,nmax-1] = k[t,nmax-2] + (k[t,nmax-2] - k[t,nmax-3])
    
    smax = 0.5
    for i in range(0,nmax-1):
        ul[t,i+1] = U[t,i] + k[t,i]/2
        ur[t,i+1] = U[t,i+1] - k[t,i+1]/2
        phi[t,i+1] = a*(ul[t,i+1] + ur[t,i+1])/2
        F[t,i+1] = a*(ul[t,i+1] + ur[t,i+1] + smax*(ul[t,i+1] - ur[t,i+1]))/2#算术平均
        
    #边界条件，外插法
    ul[t,0] = U[t,0] - k[t,0]/2
    F[t,0] = a*ul[t,0]
    ur[t,nmax] = U[t,nmax-1] + k[t,nmax-1]/2
    F[t,nmax] = a*ur[t,nmax]

for t in range(0,tmax-1):

    U2[t,:] = U[t,:]
    rk2()
    for i in range(0,nmax):
        U[t,i] = U[t,i] + (F[t,i] - F[t,i+1])
    U3[t,:] = U[t,:]
    rk2()
    for i in range(0,nmax):
        U[t+1,i] = (U2[t,i] + U[t,i] + (F[t,i] - F[t,i+1]))/2
    U[t,:] = U2[t,:]
    
    