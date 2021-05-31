import numpy as np
import matplotlib.pyplot as plt
'''
levelsetOldDemo系列均是参考自如下
https://github.com/JesseLu/level-set
写完了，效果未知
'''
nmax = 80
mmax = 80
phi = np.zeros((nmax,mmax))
phi[:,:] = 1000

def CreateCircle(xcenter,ycenter,radius):
    for i in range(0,nmax):
        for j in range(0,nmax):
            phi[i,j] = min(np.sqrt((i - xcenter)**2 + (j - ycenter)**2) - radius,phi[i,j])
            
CreateCircle(50,40,3)
CreateCircle(30,40,3)

space = 1
# Compute the (smoothed) signed function S.
S = np.zeros((nmax,mmax))
dxforward = np.zeros((nmax,mmax))
dxbackward = np.zeros((nmax,mmax))
dxcentral = np.zeros((nmax,mmax))
dyforward = np.zeros((nmax,mmax))
dybackward = np.zeros((nmax,mmax))
dycentral = np.zeros((nmax,mmax))

dxx = np.zeros((nmax,mmax))
dyy = np.zeros((nmax,mmax))
dxy = np.zeros((nmax,mmax))

phix2 = np.zeros((nmax,mmax))
phiy2 = np.zeros((nmax,mmax))
dt = 0.5
g = np.zeros((nmax,mmax))

def SignedDistance():
    
    for i in range(0,nmax):
        for j in range(0,mmax):
            S[i,j] = phi[i,j]/np.sqrt(phi[i,j]*phi[i,j] + space*space)

    for k in range(0,10):
        for i in range(0,nmax):
            for j in range(0,mmax):
                ip = int((i + 1)%nmax)
                im = int((i - 1 + nmax)%nmax)
                dxforward[i,j] = phi[ip,j] - phi[i,j]
                dxcentral[i,j] = (phi[ip,j] - phi[im,j])/2
                dxbackward[i,j] = phi[i,j] - phi[im,j]
            
                jp = int((j + 1)%nmax)
                jm = int((j - 1 + nmax)%nmax)
                dyforward[i,j] = phi[i,jp] - phi[i,j]
                dycentral[i,j] = (phi[i,jp] - phi[i,jm])/2
                dybackward[i,j] = phi[i,j] - phi[i,jm]
            
                term1 = dxbackward[i,j]**2*(dxbackward[i,j] > 0)
                term2 = dxforward[i,j]**2*(dxforward[i,j] < 0)
                term3 = dxbackward[i,j]**2*(dxbackward[i,j] < 0)
                term4 = dxforward[i,j]**2*(dxforward[i,j] > 0)
                phix2[i,j] = (S[i,j] >= 0)*max(term1,term2) + (S[i,j] < 0)*max(term3,term4)
            
                term1 = dybackward[i,j]**2*(dybackward[i,j] > 0)
                term2 = dyforward[i,j]**2*(dyforward[i,j] < 0)
                term3 = dybackward[i,j]**2*(dybackward[i,j] < 0)
                term4 = dyforward[i,j]**2*(dyforward[i,j] > 0)
                phiy2[i,j] = (S[i,j] >= 0)*max(term1,term2) + (S[i,j] < 0)*max(term3,term4)
            
                g[i,j] = np.sqrt(phix2[i,j] + phiy2[i,j])
            
        # reinitialization equation
        # phi_t + S(phi)(|nabla phi| - 1) = 
        err = 0
        for i in range(0,nmax):
            for j in range(0,mmax):
                phi[i,j] = phi[i,j] - dt * S[i,j] * (g[i,j] - 1)
                err += (g[i,j] - 1)**2
        errnorm = np.sqrt(err)
        if k % 10 == 0:
            print('%d which is :%f' % (k, errnorm))
        if errnorm < 1e-5:
            return
                
SignedDistance()
# 速度
Vx = np.zeros((nmax,mmax))
Vy = np.zeros((nmax,mmax))
Vx[0:nmax//2,:] = -0.1
Vx[nmax//2:nmax,:] = 0.1

Hext = np.zeros((nmax,mmax))
kdenom = np.zeros((nmax,mmax))
curvature = np.zeros((nmax,mmax))
Vy = np.zeros((nmax,mmax))

# 开始迭代
for t in range(0,100):
    # Update Interface
    for i in range(1,nmax-1):
        for j in range(1,mmax-1):
            dxforward[i,j] = phi[i+1,j] - phi[i,j]
            dxcentral[i,j] = (phi[i+1,j] - phi[i-1,j])/2
            dxbackward[i,j] = phi[i,j] - phi[i-1,j]
            
            dyforward[i,j] = phi[i,j+1] - phi[i,j]
            dycentral[i,j] = (phi[i,j+1] - phi[i,j-1])/2
            dybackward[i,j] = phi[i,j] - phi[i,j-1]
            
            dxx[i,j] = dxforward[i,j] - 2*phi[i,j] + dxbackward[i,j]
            dyy[i,j] = dyforward[i,j] - 2*phi[i,j] + dybackward[i,j]
            dxy[i,j] = (phi[i+1,j+1] - phi[i+1,j-1] - phi[i-1,j+1] + phi[i-1,j-1])/4
            
            Hext[i,j] = Vx[i,j]*((Vx[i,j] >= 0)*dxbackward[i,j]
                     + (Vx[i,j] < 0)*dxforward[i,j] 
            ) + Vy[i,j]*((Vy[i,j] >= 0)*dybackward[i,j] + (Vy[i,j] < 0)*dyforward[i,j])
            
            term1 = dxbackward[i,j]**2*(dxbackward[i,j] > 0)
            term2 = dxforward[i,j]**2*(dxforward[i,j] < 0)
            term3 = dxbackward[i,j]**2*(dxbackward[i,j] < 0)
            term4 = dxforward[i,j]**2*(dxforward[i,j] > 0)
            phix2[i,j] = (S[i,j] >= 0)*max(term1,term2) + (S[i,j] < 0)*max(term3,term4)
            
            term1 = dybackward[i,j]**2*(dybackward[i,j] > 0)
            term2 = dyforward[i,j]**2*(dyforward[i,j] < 0)
            term3 = dybackward[i,j]**2*(dybackward[i,j] < 0)
            term4 = dyforward[i,j]**2*(dyforward[i,j] > 0)
            phiy2[i,j] = (S[i,j] >= 0)*max(term1,term2) + (S[i,j] < 0)*max(term3,term4)
            
            g[i,j] = np.sqrt(phix2[i,j] + phiy2[i,j])
            
            kdenom[i,j] = (dxcentral[i,j]**2 + dycentral[i,j]**2)
            if kdenom[i,j] < 1e-8:
                kdenom[i,j] = 1
                
            curvature[i,j] = (dxcentral[i,j]**2*dyy[i,j] 
                - 2*dxcentral[i,j]*dycentral[i,j]*dxy[i,j] 
                + dycentral[i,j]**2*dxx[i,j])/kdenom[i,j]
    a = 0
    b = 0
    Hnormal = a * g
    alpha = 0.9
    maxH = 0
    for i in range(0,nmax):
        for j in range(0,mmax):
            good = (phi[i,j] >= -3) & (phi[i,j] <= 3)
            maxH = max(maxH,abs(Vx[i,j]*good)+abs(Vy[i,j]*good))
    dt = alpha / maxH 
    phi = phi - dt * ((Hext + Hnormal) - b * curvature)
    SignedDistance()
    plt.imshow(phi, interpolation='bilinear')
    plt.show()
    a = 0