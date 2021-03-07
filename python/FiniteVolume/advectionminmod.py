import numpy as np
#minmod
nmax = 20
tmax = 10
U = np.zeros((tmax,nmax))
U2 = np.zeros((tmax,nmax))
U3 = np.zeros((tmax,nmax))
k = np.zeros((tmax,nmax))
phi = np.zeros((tmax,nmax))
a = 1#速度
dx = 1
b = np.zeros((2,nmax+1))

for i in range(0,nmax+1):
    b[0,i] = -i*i/2 + i*10
    b[1,i] = -i*i*i/6 + i*i*5
    j = i + 0.5
    #b[1,i] = -j*j*j/6 + j*j*5
for i in range(0,nmax):
    if i <= 8:
        U[0,i] = -i*i/16 + i
    elif i >= 16:
        U[0,i] = -(i-8)*(i-8)/16 + (i-8)
    else:
        U[0,i] = 4
    
def minmod(duR,duL):
    s = (np.sign(duR) + np.sign(duL))/2
    if abs(s) == 1:
        return s*min(abs(duR),abs(duL))
    else:
        return 0
    
def rk2():
    for i in range(1,nmax-1):
        duR = U[t,i+1] - U[t,i]
        duL = U[t,i] - U[t,i-1]
        k[t,i] = minmod(duR,duL)
    k[t,0] = k[t,1] + (k[t,1] - k[t,2])
    k[t,nmax-1] = k[t,nmax-2] + (k[t,nmax-2] - k[t,nmax-3])

for t in range(0,tmax-1):

    U2[t,:] = U[t,:]
    rk2()
    for i in range(0,nmax):
        #U{i+1/2} = U{i} + k[t,i]/2
        #U{i-1/2} = U{i} - k[t,i]/2
        ul = U[t,i] - k[t,i]/2
        ur = U[t,i] + k[t,i]/2
        U[t,i] = U[t,i] + a*(ul - ur)
    U3[t,:] = U[t,:]
    rk2()
    for i in range(0,nmax):
        ul = U[t,i] - k[t,i]/2
        ur = U[t,i] + k[t,i]/2
        U[t+1,i] = (U2[t,i] + U[t,i] + a*(ul - ur))/2
    U[t,:] = U2[t,:]
    
    
    