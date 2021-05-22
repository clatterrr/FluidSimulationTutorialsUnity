import numpy as np
nmax = 256
tmax =100
dt = 0.4/nmax/nmax
x = np.zeros((nmax))
for i in range(0,nmax):
    x[i] = 2*np.pi/nmax*(i - nmax/2)
A = 25
B = 16
u = np.zeros((tmax,nmax))
for i in range(0,nmax):
    u[0,i] = 3*A*A*1/np.cosh(0.5*(A*(x[i]+2)))**2 
    u[0,i] += 3*B*B*1/np.cosh(0.5*(B*(x[i]+1)))**2
v = np.fft.fft(u[0,:])
k = np.zeros((nmax))
for i in range(0,nmax):
    k[i] = i
    if i == nmax//2:
        k[i] = 0
    elif i > nmax//2:
        k[i] = i - nmax
udata = u
tdata = 0
for t in range(0,tmax-1):
    g = -0.5*1j*dt*k
    ik3 = 1j * k * k * k
    te = np.zeros((nmax),dtype = complex)
    me = np.zeros((nmax),dtype = complex)
    for i in range(0,nmax):
        te[i] = dt/2*1j*k[i]**3
        me[i] = np.exp(te[i])
    E2 = me*me
    rk0 = g * np.fft.fft(np.fft.ifft(v).real**2)
    rk1 = g * np.fft.fft(np.fft.ifft(me*(v + rk0/2)).real**2)
    rk2 = g * np.fft.fft(np.fft.ifft(me*v + rk1/2).real**2)
    rk3 = g * np.fft.fft(np.fft.ifft(E2*v + me*rk2).real**2)
    v = E2 * v + (E2 * rk0 + 2 * me * (rk1 + rk2) + rk3)/6
    test = 1
    u[t+1,:] = np.fft.ifft(v).real
    
