import numpy as np
#涡方法模拟顶盖驱动流
#https://zhuanlan.zhihu.com/p/345332340
nmax = 8
psi = np.zeros((nmax,nmax))
v = np.zeros((nmax,nmax))
u = np.zeros((nmax,nmax))
omega = np.zeros((nmax,nmax))
w = np.zeros((nmax,nmax))
dx = dy = dt = 1
uwall = 1
Re = 10
for k1 in range(0,200):
    for i in range(1,nmax-1):#用中心差分计算涡量运输方程
        for j in range(1,nmax-1):
            w[i,j] = -(psi[i,j+1] - psi[i,j-1])/(2*dy)*(omega[i+1,j] - omega[i-1,j])/(2*dx)
            w[i,j] += (psi[i+1,j] - psi[i-1,j])/(2*dx)*(omega[i,j+1] - omega[i,j-1])/(2*dy)
            w[i,j] += 1/Re*(omega[i+1,j] -2*omega[i,j] + omega[i-1,j])/(dx*dx)
            w[i,j] += 1/Re*(omega[i,j+1] -2*omega[i,j] + omega[i,j-1])/(dy*dy)
    for i in range(1,nmax-1):
        for j in range(1,nmax-1):
            omega[i,j] += dt*w[i,j]
    omega[0,:] = -2*psi[1,:]/(dx*dx)#左面
    omega[nmax-1,:] = -2*psi[nmax-2,:]/(dx*dx)#右面
    omega[:,nmax-1] = -2*(psi[:,nmax-2])/(dy*dy) - uwall*2/dy#上面
    omega[:,0] = -2*psi[:,1]/(dy*dy)#下面  
    for k2 in range(0,100):#泊松方程求解流函数
        for i in range(1,nmax-1):
            for j in range(1,nmax-1):
                term = (psi[i+1,j] + psi[i-1,j])*dy*dy + (psi[i,j+1] + psi[i,j-1])*dx*dx
                psi[i,j] = (term + omega[i,j])/(2*(dx*dx + dy*dy))
    for i in range(1,nmax-1):
        for j in range(1,nmax-1):
            u[i,j] = (psi[i,j+1] - psi[i,j-1])/(2*dy)
            v[i,j] = (psi[i-1,j] - psi[i+1,j])/(2*dx)
    u[:,nmax-1] = uwall