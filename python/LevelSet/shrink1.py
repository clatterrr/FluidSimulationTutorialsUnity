import numpy as np
# D:\FluidSim\FluidSim\LevelSet\Level-Set-Method-master
# 参考 https://github.com/BenStyles95/Level-Set-Method
# 水平集，效果不错
nmax = 30
tmax = 100
L = 10
mx = L / (nmax - 1)#网格的长度
p = np.zeros((nmax,nmax))
dx = np.zeros((nmax,nmax))
dy = np.zeros((nmax,nmax))
dxx = np.zeros((nmax,nmax))
dyy = np.zeros((nmax,nmax))
dxy = np.zeros((nmax,nmax))

dxp = np.zeros((nmax,nmax))#Delta x plus forward
dxm = np.zeros((nmax,nmax))#Delta x minus backward
dyp = np.zeros((nmax,nmax))
dym = np.zeros((nmax,nmax))
gradNeg = np.zeros((nmax,nmax))
gradPos = np.zeros((nmax,nmax))
tp = np.zeros((tmax,nmax,nmax))

F = np.zeros((nmax,nmax))#表面张力的大小，如果大于零则向外扩张，小于零则向内收缩
for i in range(0,nmax):
    for j in range(0,nmax):
        p[i,j] = np.sqrt((i*mx-4)*(i*mx-4) + (j*mx-4)*(j*mx-4))-2
        p[i,j] = min(p[i,j],np.sqrt((i*mx-6)*(i*mx-6) + (j*mx-6)*(j*mx-6))-2)
        
for t in range(0,tmax):
    for i in range(1,nmax-1):
        for j in range(0,nmax):
            dx[i,j] = (p[i+1,j] - p[i-1,j])/2/mx
            dxx[i,j] = (p[i+1,j] - 2*p[i,j] + p[i-1,j])/mx/mx
    dx[0,:] = dx[1,:]
    dx[nmax-1,:] = dx[nmax-2,:]
    dxx[0,:] = dxx[0,:]
    dxx[nmax-1,:] = dxx[nmax-2,:]
    for i in range(0,nmax):
        for j in range(1,nmax-1):
            dy[i,j] = (p[i,j+1] - p[i,j-1])/2/mx
            dyy[i,j] = (p[i,j+1] - 2*p[i,j] + p[i,j-1])/mx/mx
    dy[:,0] = dy[:,1]
    dy[:,nmax-1] = dy[:,nmax-1]
    dyy[0,:] = dyy[1,:]
    dyy[nmax-1,:] = dyy[nmax-2,:]
    for i in range(0,nmax):
        for j in range(1,nmax-1):
            dxy[i,j] = (dx[i,j+1] - dx[i,j-1])/2/mx
    dxy[:,0] = dxy[:,1]
    dxy[:,nmax-1] = dxy[:,nmax-2]
    #3-3 对 4-3
    F[:,:] = (dxx[:,:]*dy[:,:]*dy[:,:] -2*dx[:,:]*dy[:,:]*dxy[:,:] + 
             dyy[:,:]*dx[:,:]*dx[:,:])/pow((dx[:,:]*dx[:,:] + dy[:,:]*dy[:,:]),1.5)
    for i in range(0,nmax):
        for j in range(0,nmax):
            F[i,j] = -min(max(F[i,j],-1/mx),1/mx)/100
    
    for i in range(0,nmax-1):
        dxp[i,:] = (p[i+1,:] -  p[i,:])/mx
        dxm[i+1,:] = dxp[i,:]
    dxp[nmax-1,:] = dxp[nmax-2,:]
    dxm[0,:] = dxm[1,:]
    for j in range(0,nmax-1):
        dyp[:,j] = (p[:,j+1] - p[:,j])/mx
        dym[:,j+1] = dyp[:,j]
    dyp[:,nmax-1] = dyp[nmax-2,:]
    dym[:,0] = dym[:,1]
    tp[t,:,:] = p[:,:]
    for i in range(0,nmax):
        for j in range(0,nmax):
            gradPos[i,j] = np.sqrt(max(dxm[i,j],0)**2 + min(dxp[i,j],0)**2
                                   + max(dym[i,j],0)**2 + min(dyp[i,j],0)**2)
            gradNeg[i,j] = np.sqrt(min(dxm[i,j],0)**2 + max(dxp[i,j],0)**2
                                   + min(dym[i,j],0)**2 + max(dyp[i,j],0)**2)
            p[i,j] = p[i,j] - 0.8*(max(F[i,j],0)*gradPos[i,j] + min(F[i,j],0)*gradNeg[i,j])