import numpy as np
import matplotlib.pyplot as plt
nmax = 16
f = np.zeros((nmax))
u = np.zeros((nmax))# 0 无中断，1中断左边界，2中断中间，3中断右边界
k = np.zeros((nmax))
n = np.zeros((nmax))
b = np.zeros((nmax))
qL = np.zeros((nmax))
qR = np.zeros((nmax))
for i in range(nmax):
    if i < 8:
        f[i] = (i+1)**2 - i**2
    else:
        f[i] = -3*((i+1)**2 - i**2)
n[7] = 1
b[7] = 1
tmax = 10
u[:] = 1
ft = np.zeros((nmax,tmax+1))
flux = np.zeros((nmax))
res = np.zeros((nmax))
ft[:,0] = f[:]
dt = 0.2
cx = np.zeros((nmax))
for i in range(nmax):
    cx[i] = i
for t in range(tmax):
    res[:] = 0
    for i in range(3,nmax):
        if ((n[i-1] == 0) | (b[i-1] == 1)) & (n[i-2] == 0) & ((n[i-3] == 0) | (b[i-3] == 0)):
            k2 = f[i-1] + f[i-3] - 2*f[i-2]
            k1 = (f[i-1] - f[i-3])/2
            h = dt / 2
            dis0 = 3/2 - 2*h
            dis1 = 3/2 - h
            dis2 = 3/2
            p0 = f[i-2] + k1*dis0 + k2*dis0*dis0/2
            p1 = f[i-2] + k1*dis1 + k2*dis1*dis1/2
            p2 = f[i-2] + k1*dis2 + k2*dis2*dis2/2
            flux[i] = (p0 + 4*p1 + p2)*h/3
            res[i] = res[i] + flux[i]
            res[i-1] = res[i-1] - flux[i]
        elif ((n[i] == 0) | (b[i] == 1)) & (n[i-1] == 0) & ((n[i-2] == 0) | (b[i-2] == 0)):
            k2 = f[i] + f[i-2] - 2*f[i-1]
            k1 = (f[i] - f[i-2])/2
            h = dt / 2
            dis0 = 1/2 - 2*h
            dis1 = 1/2 - h
            dis2 = 1/2
            p0 = f[i-1] + k1*dis0 + k2*dis0*dis0/2
            p1 = f[i-1] + k1*dis1 + k2*dis1*dis1/2
            p2 = f[i-1] + k1*dis2 + k2*dis2*dis2/2
            flux[i] = (p0 + 4*p1 + p2)*h/3
            res[i] = res[i] + flux[i]
            res[i-1] = res[i-1] - flux[i]
        elif ((n[i+1] == 0) | (b[i+1] == 1)) & (n[i] == 0) & ((n[i-1] == 0) | (b[i-1] == 0)):
            k2 = f[i+1] + f[i-1] - 2*f[i]
            k1 = (f[i+1] - f[i-1])/2
            h = dt / 2
            dis0 = -1/2 - 2*h
            dis1 = -1/2 - h
            dis2 = -1/2
            p0 = f[i] + k1*dis0 + k2*dis0*dis0/2
            p1 = f[i] + k1*dis1 + k2*dis1*dis1/2
            p2 = f[i] + k1*dis2 + k2*dis2*dis2/2
            flux[i] = (p0 + 4*p1 + p2)*h/3
            res[i] = res[i] + flux[i]
            res[i-1] = res[i-1] - flux[i]
        elif ((n[i+2] == 0) | (b[i+2] == 1)) & (n[i+1] == 0) & ((n[i] == 0) | (b[i] == 0)):
            k2 = f[i+2] + f[i] - 2*f[i+1]
            k1 = (f[i+2] - f[i])/2
            h = dt / 2
            dis0 = -3/2 - 2*h
            dis1 = -3/2 - h
            dis2 = -3/2
            p0 = f[i+1] + k1*dis0 + k2*dis0*dis0/2
            p1 = f[i+1] + k1*dis1 + k2*dis1*dis1/2
            p2 = f[i+1] + k1*dis2 + k2*dis2*dis2/2
            flux[i] = (p0 + 4*p1 + p2)*h/3
            res[i] = res[i] + flux[i]
            res[i-1] = res[i-1] - flux[i]
            
    res[0] = res[1]  = res[2] = res[3]
    res[nmax-1] = res[nmax-2] = res[nmax-3] = res[nmax-4]
    for i in range(0,nmax-1):
        if n[i] == 1:
            b[i] = b[i] + u[i]*dt
            if b[i] > 1:
                b[i+1] = b[i] - 1 - u[i]*dt
                b[i] = 0
                n[i] = 0
                n[i+1] = 1
    f = f + res
    ft[:,t+1] = f