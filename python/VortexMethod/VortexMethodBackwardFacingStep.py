import numpy as np
#涡方法模拟后台阶流动(Backward Facing Step Flow)
#https://zhuanlan.zhihu.com/p/345332340
#由于使用了中心差分这种不精确的格式
#所以如果调大时间步或调多迭代次数会让数值爆炸
#注意看i = 1那一列的u和v，u是负的而v是正的，说明有漩涡产生
nmax = 16
psi = np.zeros((nmax,nmax))
psis = np.zeros((nmax,nmax))
v = np.zeros((nmax,nmax))
u = np.zeros((nmax,nmax))
omega = np.zeros((nmax,nmax))
w = np.zeros((nmax,nmax))
block = np.zeros((nmax,nmax))
dx = dy = 1
dt = 0.4
uwall = 1
Re = 30
psi[:,0] = 0
height = int(nmax/2)
for j in range(height,nmax):
    psi[0,j] = uwall*(j-height)
psi[1:nmax-1,nmax-1] = psi[0,nmax-1]
    
for k1 in range(0,100):
    for k2 in range(0,200):#泊松方程求解流函数
        for i in range(1,nmax-1):
            for j in range(1,nmax-1):
                term = (psi[i+1,j] + psi[i-1,j])*dy*dy + (psi[i,j+1] + psi[i,j-1])*dx*dx
                psis[i,j] = (term + omega[i,j])/(2*(dx*dx + dy*dy))
        psi[1:nmax-1,1:nmax-1] = psis[1:nmax-1,1:nmax-1]
        psi[nmax-1,:] = psi[nmax-2,:]
    for i in range(1,nmax-1):
        for j in range(1,nmax-1):
            u[i,j] = (psi[i,j+1] - psi[i,j-1])/(2*dy)
            v[i,j] = (psi[i-1,j] - psi[i+1,j])/(2*dx)
    for i in range(1,nmax-1):#用中心差分计算涡量运输方程
        for j in range(1,nmax-1):
            w[i,j] = -(psi[i,j+1] - psi[i,j-1])/(2*dy)*(omega[i+1,j] - omega[i-1,j])/(2*dx)
            w[i,j] += (psi[i+1,j] - psi[i-1,j])/(2*dx)*(omega[i,j+1] - omega[i,j-1])/(2*dy)
            w[i,j] += 1/Re*(omega[i+1,j] -2*omega[i,j] + omega[i-1,j])/(dx*dx)
            w[i,j] += 1/Re*(omega[i,j+1] -2*omega[i,j] + omega[i,j-1])/(dy*dy)
    for i in range(1,nmax-1):
        for j in range(1,nmax-1):
            omega[i,j] += dt*w[i,j]
    for j in range(0,nmax):
        if(j < height):
            omega[0,j] = -2*(psi[1,j] - psi[0,j])/(dx*dx)
        else:
            omega[0,j] = 0#流函数在墙壁上才是常数，但入口并不是墙壁
    omega[nmax-1,:] = omega[nmax-2,:]#出口
    omega[:,nmax-1] = 2*(psi[:,nmax-1] - psi[:,nmax-2])/(dy*dy)#上边界
    omega[:,0] = 2*(psi[:,0] - psi[:,1])/(dy*dy)#下边界     
    if(k1 % 10 == 0):
        print(k1)