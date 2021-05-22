import numpy as np

def DFT(fx):
    fx = np.asarray(fx, dtype=complex)
    nx = fx.shape[0]
    fu = np.zeros(nx, dtype=complex)
    for x in range(0,nx):
        summ = 0
        for y in range(0,nx):
            summ += fx[y]*np.exp(-2j*np.pi*x*y/nx)
        fu[x] = summ
    return fu

def IDFT(fu):
    fu = np.asarray(fu, dtype=complex)
    nx = fu.shape[0]
    fx = np.zeros(nx, dtype=complex)

    for x in range(0,nx):
        summ =0
        for y in range(0,nx):
            summ += fu[y]*np.exp(2j*np.pi*x*y/nx)
        fx[x] = summ /nx
    return fx


def FFT(fx):
    fx = np.asarray(fx, dtype=complex)
    M = fx.shape[0]
    minDivideSize = 4

    if M % 2 != 0:
        raise ValueError("the input size must be 2^n")

    if M <= minDivideSize:
        return DFT(fx)
    else:
        fx_even = FFT(fx[::2]) 
        fx_odd = FFT(fx[1::2])  
        W_ux_2k = np.exp(-2j * np.pi * np.arange(M) / M)
        f_u = fx_even + fx_odd * W_ux_2k[:M//2]
        f_u_plus_k = fx_even + fx_odd * W_ux_2k[M//2:]
        fu = np.concatenate([f_u, f_u_plus_k])

    return fu


def IFFT(fu):
    fu = np.asarray(fu, dtype=complex)
    fu_conjugate = np.conjugate(fu)

    fx = FFT(fu_conjugate)

    fx = np.conjugate(fx)
    fx = fx / fu.shape[0]

    return fx

#dh/dt + dh/dx = 0
nmax = 64
L = 1
dx = L/nmax
nu = 1/4
dt = nu*dx
M = np.zeros((nmax))
tmax = 32 # 时间总数
h = np.zeros((tmax,nmax))
for i in range(0,nmax):
    h[0,i] = 1/(4 + 3*np.cos(2*np.pi*i*dx)) # 原始波形
    # h[0,i] = np.sin(2*np.pi*i*dx)
    M[i] = i
    if (i >= nmax/2):
        M[i] -= nmax
k = 2*np.pi*M/L#波数
hhat = np.zeros((tmax,nmax),dtype = complex)
rk = np.zeros((5,nmax),dtype=complex)
hhat[0,:] = FFT(h[0,:])
for t in range(0,tmax-1):
    #Runge Kutta 4阶
    rk[0,:] = hhat[t,:]
    rk[1,:] = rk[0,:] - dt*1j*k[:]*rk[0,:]
    rk[2,:] = rk[1,:] - dt*1j*k[:]*rk[1,:]
    rk[3,:] = rk[2,:] - dt*1j*k[:]*rk[2,:]
    rk[4,:] = rk[3,:] - dt*1j*k[:]*rk[3,:]
    hhat[t+1,:] = rk[0,:]*9/24 + rk[1,:]/3 + rk[2,:]/4 + rk[4,:]/24
    h[t+1,:] = IFFT(hhat[t+1,:]).real
    '''
    # Runge Kutta 4阶也可写成下面这样子，结果是一样的，但我觉很难理解不方便记忆
    rk[0,:] = -dt*1j*k[:]*hhat[t,:]
    rk[1,:] = -dt*1j*k[:]*(hhat[t,:] + 0.5*rk[0,:])
    rk[2,:] = -dt*1j*k[:]*(hhat[t,:] + 0.5*rk[1,:])
    rk[3,:] = -dt*1j*k[:]*(hhat[t,:] + rk[2,:])
    hhat[t+1,:] = hhat[t,:] + (rk[0,:] 
            + 2*rk[1,:] + 2*rk[2,:] + rk[3,:])/6 
    # Runge Kuuta 4阶也可写成下面这样，就是泰勒公式
    rk[0,:] = hhat[t,:]
    rk[1,:] = -dt*1j*k[:]*rk[0,:]
    rk[2,:] = -dt*1j*k[:]*rk[1,:]
    rk[3,:] = -dt*1j*k[:]*rk[2,:]
    rk[4,:] = -dt*1j*k[:]*rk[3,:]
    hhat[t+1,:] = rk[0,:] + rk[1,:] + rk[2,:]/2 + rk[3,:]/6 + rk[4,:]/24
    '''
    h[t+1,:] = IFFT(hhat[t+1,:]).real