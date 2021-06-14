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
        f[i] = i
    else:
        f[i] = 3*i
n[7] = 1
b[7] = 1
tmax = 10
u[:] = 1
ft = np.zeros((nmax,tmax+1))
flux = np.zeros((nmax))
res = np.zeros((nmax))
ft[:,0] = f[:]
dt = 0.1
cx = np.zeros((nmax))
for i in range(nmax):
    cx[i] = i
for t in range(tmax):
    res[:] = 0
    for i in range(1,nmax-1):
        k1 = f[i] - f[i-1]
        k2 = f[i+1] - f[i]
        if k1 * k2 > 0:
            k[i] = min(abs(k1),abs(k2))
        qL[i] = f[i] - k[i]/2
        qR[i] = f[i] + k[i]/2
    qL[0] = qL[1]
    qL[nmax-1] = qL[nmax-2]
    qR[0] = qR[1]
    qR[nmax-1] = qR[nmax-2]
    for i in range(1,nmax-1):
        if n[i] == 1:
            flux[i] = (qR[i-2] + qR[i-2] - dt*k[i-2] + 2*k[i-2])/2*dt
        elif (n[i-1] == 1) & (b[i-1] < 1):
            flux[i] = (qR[i] + qR[i] - dt*k[i] - 2*k[i])/2*dt
        else:
            flux[i] = (qR[i-1] + qR[i-1] - dt*k[i-1])/2*dt
        # flux[i] = (qR[i-1] + qR[i-1] - dt*k[i-1])/4
        res[i-1] -= flux[i]
        res[i] += flux[i]
    res[0] = res[1]  = res[2]
    res[nmax-1] = res[nmax-2] = res[nmax-3]
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
    
from matplotlib.animation import FuncAnimation
fig = plt.figure()
def update(frame):
    return ft[:,frame]

anim = FuncAnimation(fig, update, frames=range(tmax+1), interval=1000)
anim.save("c02.gif", writer='imagemagick')
plt.show()
    