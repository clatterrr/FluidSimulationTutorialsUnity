import numpy as np
'''
二阶精度的线性对流方程，有限差分法
使用二阶龙格库塔，也就是二阶泰勒公式
'''
nmax = 8
mmax = 8
f = np.zeros((nmax,mmax))
for i in range(nmax):
    for j in range(mmax):
        f[i,j] = i * i + j * j * 3 - 2 * i * j - 5 * i + 7 * j + 5 
        
dt = 0.5
tmax = 10
u = 1
v = 1
dx = np.zeros((nmax,mmax))
dy = np.zeros((nmax,mmax))
ft = np.zeros((tmax+1,nmax,mmax))
ft[0,:,:] = f[:,:]

def RungeKutta():
    for i in range(1,nmax-1):
        for j in range(0,mmax):
            dx[i,j] = (f[i+1,j] - f[i-1,j])/2
        
    dx[0,:] = dx[1,:] + (dx[1,:] - dx[2,:])
    dx[nmax-1,:] = dx[nmax-2,:] + (dx[nmax-2,:] - dx[nmax-3,:])
    for i in range(0,nmax):
        for j in range(1,mmax-1):
            dy[i,j] = (f[i,j+1] - f[i,j-1])/2
    
    dy[:,0] = dy[:,1] + (dy[:,1] - dy[:,2])
    dy[:,nmax-1] = dy[:,nmax-2] + (dy[:,nmax-2] - dy[:,nmax-3])

for t in range(tmax):
    rk0 = f.copy()
    RungeKutta()
    f = f - dt * (dx*u + dy*v)
    RungeKutta()
    f = f - dt * (dx*u + dy*v)
    rk2 = f.copy()
    f = (rk0 + rk2)/2
    ft[t+1,:,:] = f[:,:]
            