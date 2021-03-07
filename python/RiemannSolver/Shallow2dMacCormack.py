import numpy as np
#MacCormack for 2d Shallow Water Equation
#介绍：https://zhuanlan.zhihu.com/p/331771977
nmax = 9
tmax = 100
h = np.ones((nmax+2,nmax+2))
for i in range(0,4):
    h[i,:] = 2

u = np.zeros((nmax+2,nmax+2))
v = np.zeros((nmax+2,nmax+2))
U1 = np.zeros((nmax+2,nmax+2))
U2 = np.zeros((nmax+2,nmax+2))
U3 = np.zeros((nmax+2,nmax+2))
F1 = np.zeros((nmax+2,nmax+2))
F2 = np.zeros((nmax+2,nmax+2))
F3 = np.zeros((nmax+2,nmax+2))
G1 = np.zeros((nmax+2,nmax+2))
G2 = np.zeros((nmax+2,nmax+2))
G3 = np.zeros((nmax+2,nmax+2))

u1pred = np.zeros((nmax+2,nmax+2))
u2pred = np.zeros((nmax+2,nmax+2))
u3pred = np.zeros((nmax+2,nmax+2))

dx = dy = 100.0/(nmax-1)
cfl = 0.8

for i in range(0,nmax+2):
    for j in range(0,nmax+2):
        U1[i,j] = h[i,j]
        U2[i,j] = h[i,j] * u[i,j]
        U3[i,j] = h[i,j] * v[i,j]
        F1[i,j] = h[i,j] * u[i,j]
        F2[i,j] = h[i,j] * u[i,j] * u[i,j] + 0.5 * 10 * h[i,j] * h[i,j]
        F3[i,j] = h[i,j] * u[i,j] * v[i,j]
        G1[i,j] = h[i,j] * v[i,j]
        G2[i,j] = h[i,j] * u[i,j] * v[i,j]
        G3[i,j] = h[i,j] * v[i,j] * v[i,j] + 0.5 * 10 * h[i,j] * h[i,j]



for k in range(0,tmax-1):
    dxx = np.max(abs(u)+np.sqrt(10*h))
    dyy = np.max(abs(v)+np.sqrt(10*h))
    dt = cfl * dx /  min(dxx,dyy)
    c = dt / dx
    for i in range(1,nmax+1):
        for j in range(1,nmax+1):
            u1pred[i,j] = U1[i,j] - c*(F1[i+1,j] - F1[i,j]) - c*(G1[i,j+1] - G1[i,j])
            u2pred[i,j] = U2[i,j] - c*(F2[i+1,j] - F2[i,j]) - c*(G2[i,j+1] - G2[i,j])
            u3pred[i,j] = U3[i,j] - c*(F3[i+1,j] - F3[i,j]) - c*(G3[i,j+1] - G3[i,j])
            h[i,j] = u1pred[i,j]
            u[i,j] = u2pred[i,j] / u1pred[i,j]
            v[i,j] = u3pred[i,j] / u1pred[i,j]
            
    h[0,:] = h[1,:]
    h[-1,:] = h[-2,:]
    h[:,0] = h[:,1]
    h[:,-1] = h[:,-2]
    u[0,:] = u[1,:]
    u[-1,:] = u[-2,:]
    u[:,0] = u[:,1]
    u[:,-1] = u[:,-2]
    v[0,:] = v[1,:]
    v[-1,:] = v[-2,:]
    v[:,0] = v[:,1]
    v[:,-1] = v[:,-2]
    
    for i in range(0,nmax+2):
        for j in range(0,nmax+2):
            F1[i,j] = h[i,j] * u[i,j]
            F2[i,j] = h[i,j] * u[i,j] * u[i,j] + 0.5 * 10 * h[i,j] * h[i,j]
            F3[i,j] = h[i,j] * u[i,j] * v[i,j]
            G1[i,j] = h[i,j] * v[i,j]
            G2[i,j] = h[i,j] * u[i,j] * v[i,j]
            G3[i,j] = h[i,j] * v[i,j] * v[i,j] + 0.5 * 10 * h[i,j] * h[i,j]
            
    for i in range(1,nmax+1):
        for j in range(1,nmax+1):
            U1[i,j] = 0.5*(U1[i,j] + u1pred[i,j] - c * (F1[i,j] - F1[i-1,j]) - c * (G1[i,j] - G1[i,j-1]))
            U2[i,j] = 0.5*(U2[i,j] + u2pred[i,j] - c * (F2[i,j] - F2[i-1,j]) - c * (G2[i,j] - G2[i,j-1]))
            U3[i,j] = 0.5*(U3[i,j] + u3pred[i,j] - c * (F3[i,j] - F3[i-1,j]) - c * (G3[i,j] - G3[i,j-1]))
            h[i,j] = U1[i,j]
            u[i,j] = U2[i,j]/U1[i,j]
            v[i,j] = U3[i,j]/U1[i,j]
            
    h[0,:] = h[1,:]
    h[-1,:] = h[-2,:]
    h[:,0] = h[:,1]
    h[:,-1] = h[:,-2]
    u[0,:] = u[1,:]
    u[-1,:] = u[-2,:]
    u[:,0] = u[:,1]
    u[:,-1] = u[:,-2]
    v[0,:] = v[1,:]
    v[-1,:] = v[-2,:]
    v[:,0] = v[:,1]
    v[:,-1] = v[:,-2]
    
    for i in range(0,nmax+2):
        for j in range(0,nmax+2):
            F1[i,j] = h[i,j] * u[i,j]
            F2[i,j] = h[i,j] * u[i,j] * u[i,j] + 0.5 * 10 * h[i,j] * h[i,j]
            F3[i,j] = h[i,j] * u[i,j] * v[i,j]
            G1[i,j] = h[i,j] * v[i,j]
            G2[i,j] = h[i,j] * u[i,j] * v[i,j]
            G3[i,j] = h[i,j] * v[i,j] * v[i,j] + 0.5 * 10 * h[i,j] * h[i,j]
    