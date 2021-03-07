import numpy as np
#涡方法中的流函数涡方程
#https://zhuanlan.zhihu.com/p/345332340
nmax = 8
psi = np.zeros((nmax,nmax))
v = np.zeros((nmax,nmax))
u = np.zeros((nmax,nmax))
omega = np.zeros((nmax,nmax))
dx = dy = 1
omega[3,3] = 1#涡量为正则逆时针旋转
for k in range(0,100):
    for i in range(1,nmax-1):
        for j in range(1,nmax-1):
            term = (psi[i+1,j] + psi[i-1,j])*dy*dy + (psi[i,j+1] + psi[i,j-1])*dx*dx
            psi[i,j] = (term + omega[i,j])/(2*(dx*dx + dy*dy))
for i in range(1,nmax-1):
    for j in range(1,nmax-1):
        u[i,j] = (psi[i,j+1] - psi[i,j-1])/(2*dy)
        v[i,j] = (psi[i-1,j] - psi[i+1,j])/(2*dx)