import numpy as np
"""
参考
A Front-tracking/Finite-Volume Navier-Stokes
Solver for Direct Numerical Simulations of
Multiphase Flows
Gr´etar Tryggvason
October 19, 2012
此代码最后更新日期
2021-5-23
"""
Lx = 1
Ly = 1
gx = 0
gy = -100
rho1 = 1
rho2 = 2
m0 = 0
rro = rho1
nmax = 32
mmax = 32
dt = 0.00125
tmax = 100
beta = 1.2

u = np.zeros((nmax+1,mmax+2))
ut = np.zeros((nmax+1,mmax+2))
uu = np.zeros((nmax+1,mmax+2))
v = np.zeros((nmax+2,mmax+1))
vt = np.zeros((nmax+2,mmax+1))
vv = np.zeros((nmax+2,mmax+1))
p = np.zeros((nmax+2,mmax+2))
tmp1 = np.zeros((nmax+2,mmax+2))
tmp2 = np.zeros((nmax+2,mmax+2))

unorth = 0
usouth = 0
veast = 0
vwest = 0

dx = Lx / nmax
dy = Ly / mmax
r = np.zeros((nmax+2,mmax+2))
r[:,:] = rho1
x = np.zeros((nmax + 2))
y = np.zeros((mmax + 2))
for i in range(0,nmax+2):
    x[i] = dx * (i - 0.5)
for i in range(0,mmax+2):
    y[i] = dy * (i - 0.5)
xcenter = 0.5
ycenter = 0.8
rad = 0.15
for i in range(0,nmax+2):
    for j in range(0,mmax+2):
        if ((x[i] - xcenter)**2 + (y[j] - ycenter)**2) < (rad**2):
            r[i,j] = rho2
            
# 设置前线
Nf = 100
xf = np.zeros((Nf + 2))
yf = np.zeros((Nf + 2))
uf = np.zeros((Nf + 2))
vf = np.zeros((Nf + 2))
tx = np.zeros((Nf + 2))
ty = np.zeros((Nf + 2))
            
for i in range(0,Nf + 2):
    xf[i] = xcenter - rad * np.sin(2*np.pi*i / Nf)
    yf[i] = ycenter + rad * np.cos(2*np.pi*i / Nf)

for t in range(0,tmax):
    u[:,0] = 2 * usouth - u[:,1]
    u[:,mmax+1] = 2 * unorth - u[:,mmax]
    v[0,:] = 2 * vwest - v[1,:]
    v[nmax+1,:] = 2 * veast - v[nmax,:]
    
    for i in range(1,nmax):
        for j in range(1,mmax+1):
            term1 = ((u[i+1,j] + u[i,j])**2 - (u[i,j] + u[i-1,j])**2)/dx
            term2 = (u[i,j+1] + u[i,j])*(v[i+1,j] + v[i,j])
            term3 = (u[i,j] + u[i,j-1])*(v[i+1,j-1] + v[i,j-1])
            term4 = (term2 - term3)/dy
            
            term5 = m0 / (0.5 * (r[i+1,j] + r[i,j]))
            term6 = (u[i+1,j] - 2*u[i,j] + u[i-1,j])/dx/dx
            term7 = (u[i,j+1] - 2*u[i,j] + u[i,j-1])/dy/dy
            
            ut[i,j] = u[i,j] + dt * (-0.25*(term1 + term4) + term5*(term6 + term7) + gx)
    
    for i in range(1,nmax+1):
        for j in range(1,mmax):
            term1 = (u[i,j+1] + u[i,j])*(v[i+1,j] + v[i,j])
            term2 = (u[i-1,j+1] + u[i-1,j])*(v[i,j] + v[i-1,j])
            term3 = (term1 - term2)/dx
            term4 = ((v[i,j+1] + v[i,j])**2 - (v[i,j] + v[i,j-1])**2)/dy
            
            term5 = m0 / (0.5 * (r[i,j+1] + r[i,j]))
            term6 = (v[i+1,j] - 2*v[i,j] + v[i-1,j])/dx/dx
            term7 = (v[i,j+1] - 2*v[i,j] + v[i,j-1])/dy/dy
            
            vt[i,j] = v[i,j] + dt*(-0.25*(term3 + term4) + term5*(term6 + term7) + gy)
            
    rt = r.copy()
    rlarge = 1000
    rt[0,:] = rt[nmax+1,:] = rt[:,0] = rt[:,mmax+1] = rlarge
    
    for i in range(1,nmax+1):
        for j in range(1,mmax+1):
            tmp1[i,j] = (0.5/dt)*((ut[i,j] - ut[i-1,j])/dx + (vt[i,j] - vt[i,j-1])/dy)
            
            term1 = 1 / (dx * (rt[i+1,j] + rt[i,j]))
            term2 = 1 / (dx * (rt[i-1,j] + rt[i,j]))
            term3 = 1 / (dy * (rt[i,j+1] + rt[i,j]))
            term4 = 1 / (dy * (rt[i,j-1] + rt[i,j]))
            
            tmp2[i,j] = 1 / ( 1 / dx * (term1 + term2) + 1 / dy * (term3 + term4))
            
    kmax = 200
    for k in range(0,kmax):
        oldp = p.copy()
        
        for i in range(1,nmax+1):
            for j in range(1,mmax+1):
                term1 = p[i+1,j] / (dx * (rt[i+1,j] + rt[i,j]))
                term2 = p[i-1,j] / (dx * (rt[i-1,j] + rt[i,j]))
                term3 = p[i,j+1] / (dy * (rt[i,j+1] + rt[i,j]))
                term4 = p[i,j-1] / (dy * (rt[i,j-1] + rt[i,j]))
                
                term5 = 1 / dx * (term1 + term2) + 1 / dy * (term3 + term4)
                
                p[i,j] = (1 - beta)*p[i,j] + beta*tmp2[i,j]*(term5 - tmp1[i,j])
                
    for i in range(1,nmax):
        for j in range(1,mmax+1):
            u[i,j] = ut[i,j] - dt * (2.0 / dx)*(p[i+1,j] - p[i,j])/(r[i+1,j] + r[i,j])
            
    for i in range(1,nmax+1):
        for j in range(1,mmax):
            v[i,j] = vt[i,j] - dt * (2.0 / dy)*(p[i,j+1] - p[i,j])/(r[i,j+1] + r[i,j])
            
    test = 1
    # Advect Front
    for i in range(1,Nf + 1):
        ip = int(np.floor(xf[i] / dx))
        jp = int(np.floor((yf[i] + 0.5 * dy)/dy))
        ax = xf[i] / dx - ip 
        ay = (yf[i] + 0.5*dy)/dy - jp
        uf[i] = (1 - ax)*(1 - ay)*u[ip,jp] + ax*(1 - ay)*u[ip+1,jp]
        uf[i] += (1 - ax)*ay*u[ip,jp+1] + ax*ay*u[ip+1,jp+1]
        
        ip = int(np.floor((xf[i] + 0.5*dx)/dx))
        jp = int(np.floor(yf[i] / dy))
        ax = (xf[i] + 0.5 * dx)/dx - ip
        ay = yf[i] / dy - jp
        vf[i] = (1 - ax)*(1 - ay)*v[ip,jp] + ax*(1 - ay)*v[ip+1,jp]
        vf[i] += (1 - ax)*ay*v[ip,jp+1] + ax*ay*v[ip+1,jp+1] 
        
    for i in range(1,Nf + 1):
        xf[i] = xf[i] + dt * uf[i]
        yf[i] = yf[i] + dt * vf[i]
    xf[0] = xf[Nf]
    yf[0] = yf[Nf]
    xf[Nf + 1] = xf[1]
    yf[Nf + 1] = yf[1]
    
    xfold = xf.copy()
    yfold = yf.copy()
    j = 0
    for i in range(1,Nf + 1):
        ds = np.sqrt(((xfold[i] - xf[j])/dx)**2 + ((yfold[i] - yf[j])/dy)**2)
        if ds > 0.5:
            j = j + 1
            xf[j] = 0.5 * (xfold[i] + xf[j - 1])
            yf[j] = 0.5 * (yfold[i] + yf[j - 1])
            j = j + 1
            xf[j] = xfold[i]
            yf[j] = yfold[j]
        elif ds < 0.25:
            # 什么都不做
            ds = ds
        else:
            j = j + 1
            xf[j] = xfold[i]
            yf[j] = yfold[i]
    Nf = j
    xf[0] = xf[Nf]
    yf[0] = yf[Nf]
    xf[Nf + 1] = xf[1]
    yf[Nf + 1] = yf[1]
            
    fx = np.zeros((nmax + 2,mmax + 2))
    fy = np.zeros((nmax + 2,mmax + 2))
            
    for i in range(1,Nf + 1):
        nfx = -0.5 * (yf[i + 1] - yf[i - 1])*(rho2 - rho1)
        nfy = 0.5*(xf[i + 1] - xf[i - 1])*(rho2 - rho1)
        ip = int(np.floor(xf[i] / dx)) 
        jp = int(np.floor((yf[i] + 0.5 * dy)/dy))
        ax = xf[i] / dx - ip
        ay = (yf[i] + 0.5 * dy)/dy - jp
        fx[ip,jp] = fx[ip,jp] + (1-ax)*(1-ay)*nfx/dx/dy
        fx[ip+1,jp] = fx[ip+1,jp] + ax*(1-ay)*nfx/dx/dy
        fx[ip,jp+1] = fx[ip,jp+1] + (1-ax)*ay*nfx/dx/dy
        fx[ip+1,jp+1] = fx[ip+1,jp+1] + ax*ay*nfx/dx/dy
        
        ip = int(np.floor((xf[i] + 0.5*dx)/dx))
        jp = int(np.floor(yf[i] / dy))
        ax = (xf[i] + 0.5*dx)/dx - ip
        ay = yf[i] / dy - jp
        fy[ip,jp] = fy[ip,jp] + (1-ax)*(1-ay)*nfy/dx/dy
        fy[ip+1,jp] = fy[ip+1,jp] + ax*(1-ay)*nfy/dx/dy
        fy[ip,jp+1] = fy[ip,jp+1] + (1-ax)*ay*nfy/dx/dy
        fy[ip+1,jp+1] = fy[ip+1,jp+1] + ax*ay*nfy/dx/dy
        
    rt = r.copy()
    
    for k in range(0,200):
        for i in range(1,nmax+1):
            for j in range(1,mmax+1):
                r[i,j] = 0.25*(r[i+1,j] + r[i-1,j] + r[i,j+1] + r[i,j-1])
                r[i,j] += 0.25*(dx*fx[i-1,j] - dx*fx[i,j] + dy*fy[i,j-1] - dy*fy[i,j])
            
    test = 1