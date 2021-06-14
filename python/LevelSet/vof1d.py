import numpy as np
nmax = 128
tmax = 500
f = np.zeros((nmax))
normal = np.zeros((nmax)) # 0 水平 ， 1 垂直即间断
f[:] = 1
f[4:16] = 2
u = np.ones((nmax))
k = np.zeros((nmax))
flux = np.zeros((nmax))
res = np.zeros((nmax))
qL = np.zeros((nmax))
qR = np.zeros((nmax))
ft = np.zeros((nmax,tmax+1))
ft[:,0] = f[:]
def superbee(a,b,c):
    a0 = np.sign(a)
    b0 = np.sign(b)
    c0 = np.sign(c)
    re = 0
    if (a0 == b0) & (b0 == c0):
        re = a0 * min(abs(a),abs(b),abs(c))
    return re
def minmod(a,b):
    a0 = np.sign(a)
    b0 = np.sign(b)
    re = 0
    if (a0 == b0):
        re = a0 * min(abs(a),abs(b))
    return re  

def RungeKutta():
    res[:] = 0
    for i in range(1,nmax-1):
        # k[i] = superbee(f[i] - f[i-1],(f[i+1] - f[i-1])/2,f[i+1] - f[i])
        k[i] = minmod(f[i] - f[i-1],f[i+1] - f[i])
    k[0] = k[1]
    k[nmax-1] = k[nmax-2]
    for i in range(0,nmax-1):
        qL[i] = f[i] + k[i]/2
        qR[i] = f[i+1] - k[i+1]/2
        flux[i] = (u[i]*f[i] + u[i+1]*f[i+1] + (qL[i] - qR[i]))/2
        res[i] = res[i] + flux[i]
        res[i+1] = res[i+1] - flux[i]
    res[0] = 0
    u[nmax-1] = -1
    
import matplotlib.pyplot as plt
cx = np.zeros((nmax))
for i in range(nmax):
    cx[i] = i
for t in range(0,tmax):
    dt = 0.8
    rk0 = f.copy()
    RungeKutta()
    f = f - dt*res
    rk1 = f.copy()
    RungeKutta()
    f = f - dt*res
    rk2 = f.copy()
    RungeKutta()
    f = f - dt*res
    rk3 = f.copy()
    RungeKutta()
    f = f - dt*res
    rk4 = f.copy()
    f = rk4/24 + rk2/4 + rk1/3 + rk0*9/24
    ft[:,t+1] = f[:]
    plt.plot(cx,f)
    plt.ylim([0,3])
    plt.pause(0.1)
    