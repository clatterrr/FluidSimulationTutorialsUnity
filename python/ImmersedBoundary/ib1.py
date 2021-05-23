import numpy as np
'''
浸入边界法Immersed Boundary
参考自特别好的库
https://github.com/nickabattista/IB2d

'''
nmax = 64
mmax = 64
Lx = 1
Ly = 1
dx = Lx / nmax
dy = Ly / mmax
supp = 4 #  Delta-function support

mu = 1.8e-5
rho = 1.225
dt = 1e-4
TFinal = 0.4750
Tcurrent = 0

# Lagrangian Spacing
ds = min(Lx/(2*nmax),Ly / (2*mmax))
x = np.zeros((nmax))
y = np.zeros((mmax))
u = np.zeros((nmax,mmax))
v = np.zeros((nmax,mmax))
for i in range(0,nmax):
    x[i] = i * dx
for j in range(0,mmax):
    y[j] = j * dy
pnum = 126
posx = np.zeros((pnum))
posy = np.zeros((pnum))
for i in range(0,pnum):
    if i < 25:
        posx[i] = 0.1428
        posy[i] = 0.78 - 0.0078*i
    else:
        posx[i] = 0.1428 + (i - 24)*0.0072
        posy[i] = 0.05

Xind = np.zeros((pnum,16))
Yind = np.zeros((pnum,16))
for t in range(0,1):
    posxPrev = posx.copy()
    posyPrev = posy.copy()
    uPrev = u.copy()
    vPrev = v.copy()
    
    Xind[:,1] = np.floor(posx / dx) + 1
    Xind[:,5] = Xind[:,9] = Xind[:,13] = Xind[:,1]
    Xind[:,0] = Xind[:,4] = Xind[:,8] = Xind[:,12] = Xind[:,1] - 1
    Xind[:,2] = Xind[:,6] = Xind[:,10] = Xind[:,14] = Xind[:,1] + 1
    Xind[:,3] = Xind[:,7] = Xind[:,11] = Xind[:,15] = Xind[:,1] + 2
    
    Yind[:,1] = np.floor(posy / dy) + 1
    Yind[:,5] = Yind[:,9] = Yind[:,13] = Yind[:,1]
    Yind[:,0] = Yind[:,4] = Yind[:,8] = Yind[:,12] = Yind[:,1] - 1
    Yind[:,2] = Yind[:,6] = Yind[:,10] = Yind[:,14] = Yind[:,1] + 1
    Yind[:,3] = Yind[:,7] = Yind[:,11] = Yind[:,15] = Yind[:,1] + 2
    
    xResize = np.zeros((pnum,supp**2))
    yResize = np.zeros((pnum,supp**2))
    distX = np.zeros((pnum,supp**2))
    distY = np.zeros((pnum,supp**2))
    for i in range(0,supp**2):
        xResize[:,i] = posx[:]
        yResize[:,i] = posy[:]
    for i in range(0,pnum):
        for j in range(0,supp**2):
            dis = abs(x[int(Xind[i,j])] - xResize[i,j])
            distX[i,j] = min(dis,Lx - dis)
            dis = abs(y[int(Yind[i,j])] - yResize[i,j])
            distY[i,j] = min(dis,Ly - dis)
            
    deltaX = np.zeros((pnum,supp**2))
    deltaY = np.zeros((pnum,supp**2))
    for i in range(0,pnum):
        for j in range(0,supp**2):
            r = abs(distX[i,j])/dx
            if r < 1:
                deltaX[i,j] = ((3 - 2*r + np.sqrt(1 + 4*r - 4*r*r))/(8*dx))
            elif ((r < 2) & (r >= 1)):
                deltaX[i,j] = ((5 - 2*r + np.sqrt(-7 + 12*r - 4*r*r))/(8*dx))
            r = abs(distY[i,j])/dx
            if r < 1:
                deltaY[i,j] = ((3 - 2*r + np.sqrt(1 + 4*r - 4*r*r))/(8*dx))
            elif ((r < 2) & (r >= 1)):
                deltaY[i,j] = ((5 - 2*r + np.sqrt(-7 + 12*r - 4*r*r))/(8*dx))
    
    moveX = np.zeros((pnum,supp**2))
    moveY = np.zeros((pnum,supp**2))
    moveXsum = 0
    moveYsum = 0
    for i in range(0,pnum):
        for j in range(0,supp**2):
            idx = int(Xind[i,j])
            idy = int(Yind[i,j])
            
            moveX[i,j] = u[idx,idy]*deltaX[i,j]*deltaY[i,j]
            moveY[i,j] = v[idx,idy]*deltaX[i,j]*deltaY[i,j]
        moveXsum += moveX[i,1] / dx / dy
        moveYsum += moveY[i,1] / dx / dy
        
    posxNow = np.zeros((pnum))
    posyNow = np.zeros((pnum))
    posxNow = posx + dt * moveXsum
    posyNow = posx + dt * moveYsum