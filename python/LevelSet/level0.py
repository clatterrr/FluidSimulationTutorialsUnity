import numpy as np
import scipy.io as scio
'''
参考自
Distance Regularized Level Set Evolution and Its Application to Image Segmentation
https://github.com/balcilar/DRLSE-Image-Segmentation
最后更新时间
2021-5-23
'''
dataFile = 'image0.mat'
data = scio.loadmat(dataFile)
img = data['Img_smooth']
nmax = img.shape[0]
mmax = img.shape[1]

gradx = np.zeros((nmax,mmax))#x方向上的梯度
grady = np.zeros((nmax,mmax))#y方向上的梯度
f = np.zeros((nmax,mmax))# 梯度的长度
edge = np.zeros((nmax,mmax))# 边缘指示器

# 计算梯度

def gradient(fun):
    gx = np.zeros((nmax,mmax))
    gy = np.zeros((nmax,mmax))
    for i in range(1,nmax-1):
        gx[i,:] = (fun[i+1,:] - fun[i-1,:])/2
    for j in range(1,mmax-1):
        gy[:,j] = (fun[:,j+1] - fun[:,j-1])/2
    gx[0,:] = fun[1,:] - fun[0,:]
    gx[nmax-1,:] = fun[nmax-1,:] - fun[nmax-2,:]
    gy[:,0] = fun[:,1] - fun[:,0]
    gy[:,nmax-1] = fun[:,nmax-1] - fun[:,nmax-2] 
    return gx,gy

def del2(fun):
    res = np.zeros((nmax,mmax))
    for i in range(1,nmax-1):
        for j in range(1,mmax-1):
            res[i,j] = fun[i+1,j] + fun[i-1,j] + fun[i,j+1] + fun[i,j-1] - 4*fun[i,j]
    res[0,:] = res[1,:]
    res[nmax-1,:] = res[nmax-2,:]
    res[:,0] = res[:,1]
    res[:,nmax-1] = res[:,nmax-2]
    return res
gradx,grady = gradient(img)
for i in range(0,nmax):
    for j in range(0,mmax):
        f[i,j] = gradx[i,j]**2 + grady[i,j]**2
        edge[i,j]= 1 / (1 + f[i,j])
        
c0 = 2 # 水平集方向初始值
phi = np.zeros((nmax,mmax))
phi[:,:] = c0
phi[24:35,19:25] = -c0
phi[24:35,39:50] = -c0

dt = 1
mu = 0.2 / dt
lamb = 5
alpha = -3

tmax = 100
for t in range(0,tmax):
    vx,vy = gradient(edge)
    kmax = 10
    for k in range(0,kmax):
        # 边界条件
        phi[0,:] = phi[1,:]
        phi[nmax-1,:] = phi[nmax-2,:]
        phi[:,0] = phi[:,1]
        phi[:,nmax-1] = phi[:,nmax-2]
        phix,phiy = gradient(phi)
        s = np.sqrt(phix**2 + phiy**2)
        small = 1e-10
        normalx = phix / (s + small)
        normaly = phiy / (s + small)
        
        nxx,junk = gradient(normalx)
        junk,nyy = gradient(normaly)
        curvature = nxx + nyy
        
        # Double wall
        mata = np.zeros((nmax,mmax))
        matb = np.zeros((nmax,mmax))
        ps = np.zeros((nmax,mmax))
        dps = np.zeros((nmax,mmax))
        for i in range(0,nmax):
            for j in range(0,mmax):
                if (s[i,j] >= 0) & (s[i,j] <= 1):
                    mata[i,j] = 1
                if s[i,j] > 1:
                    matb[i,j] = 1
                ps[i,j] = mata[i,j]*np.sin(2*np.pi*s[i,j]
                    )/(2*np.pi) + matb[i,j]*(s[i,j] - 1)
                ps0 = (ps[i,j] == 0)
                s0 = (s[i,j] == 0)
                dps[i,j] = ((1-ps0)*ps[i,j]+ps0)/((1-s0)*s[i,j]+s0)
        # Single-Wall
        # distRegTerm = del2(phi) - curvature        
        # Double-Wall
        nxx,junk = gradient(dps*phix - phix)
        junk,nyy = gradient(dps*phiy - phiy)
        d4 =  del2(phi)
        distRegTerm = nxx + nyy 
        distRegTerm += d4
        
        #Dirac
        eps = 1.5
        diracPhi = np.zeros((nmax,mmax))
        for i in range(0,nmax):
            for j in range(0,mmax):
                if (phi[i,j] <= eps) & (phi[i,j] >= -eps):
                    diracPhi[i,j] = (1/2/eps)*(1 + np.cos(np.pi*phi[i,j]/eps))
        areaTerm = diracPhi * edge
        edgeTerm = diracPhi * (vx * normalx + vy * normaly) + diracPhi * edge * curvature
        phi = phi + dt * (mu * distRegTerm + lamb * edgeTerm + alpha * areaTerm)
        debug = 1
