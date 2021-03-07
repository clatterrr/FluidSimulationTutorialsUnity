import numpy as np
#Godunov Scheme for 1d Burgers equation
#介绍：https://zhuanlan.zhihu.com/p/331771977
nmax = 8
tmax = 100
dx = 1.0 / nmax
dt = 1.0 / tmax
U = np.zeros((tmax,nmax))
F = np.zeros((tmax,nmax + 1))
for i in range(0,nmax):
    if(i < 2):
        U[0,i] = 2
    else:
        U[0,i] = 1

cfl = 0.8
dx = 0.1
dt = 0.01

def Flux(uL,uR):# Godunov
    FL = 0.5 * uL * uL
    FR = 0.5 * uR * uR
    s = 0.5*(uL + uR)
    if (uL < uR):
        if (uL > 0.0):   #对应第一种情况
            return FL
        elif (uR < 0.0): #对应第一种情况
            return FR
        else:            #对应第三种情况
            return 0.0
    else:
        if (s > 0.0):    #对应第四种情况
            return FL
        else:            #对应第五种情况
            return FR 

for k in range(0,tmax-1):
    for i in range(1,nmax):
        uL = U[k,i-1]
        uR = U[k,i]
        F[k,i] = Flux(uL,uR)
        
    if(U[k,0] < 0.0):
        uL = 2.0 * U[k,0] - U[k,1]
    else:
        uL = U[k,0]
    uR = U[k,0]
    F[k,0] = Flux(uL,uR)
    
    if(U[k,nmax-1] > 0.0):
       uR = 2.0 * U[k,nmax-1] - U[k,nmax-2]
    else:
       uR = U[k,nmax-1]
    uL = U[k,nmax-1]
    F[k,nmax] = Flux(uL,uR)
    
    for i in range(0,nmax):
        U[k+1,i] = U[k,i] - dt/dx * (F[k,i+1] - F[k,i])