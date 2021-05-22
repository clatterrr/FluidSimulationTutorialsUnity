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


def derivation(q):
    qhat = FFT(q)
    qh = np.fft.fft(q)
    qx = IFFT(1j*k[:]*qhat[:]).real
    qxx = IFFT(-k[:]**2*qhat[:]).real
    qxxx = IFFT(-1j*k[:]**3*qhat[:]).real
    a = 1 # 系数
    res = a*qx # 对流方程
    # res = - a*qxx # 扩散方程
    res = 6*q*qx + qxxx # kdv 方程
    # res = q*qx # 无粘性burgers 方程
    # res = q*qx - a*qxx # 粘性burgers 方程
    return res

L = 16
nmax = 32
tmax = 100
amp = 2
b = np.sqrt(amp/2)  # 孤子宽度的倒数
h = np.zeros((tmax,nmax))
M = np.zeros((nmax))
rk = np.zeros((4,nmax))
dt = 0.01
for i in range(0,nmax):
    h[0,i] = amp*1/np.cosh(b*i-4)**2
    M[i] = i
    if (i >= nmax/2):
        M[i] -= nmax
k = 2*np.pi*M/L#波数
for t in range(0,tmax-1):
    rk[0,:] = -dt*derivation(h[t,:])
    rk[1,:] = -dt*derivation(h[t,:] + 0.5*rk[0,:])
    rk[2,:] = -dt*derivation(h[t,:] + 0.5*rk[1,:])
    rk[3,:] = -dt*derivation(h[t,:] + rk[2,:])
    h[t+1,:] = h[t,:] + (rk[0,:] 
            + 2*rk[1,:] + 2*rk[2,:] + rk[3,:])/6 