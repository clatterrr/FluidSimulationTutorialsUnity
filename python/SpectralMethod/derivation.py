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
    qx = IFFT(1j*k[:]*qhat[:]).real
    qxx = IFFT(-k[:]**2*qhat[:]).real
    qxxx = IFFT(-1j*k[:]**3*qhat[:]).real
    a = 0.1 # 扩散系数
    res = - a*qxx # 扩散方程
    res = 6*q*qx + qxxx # kdv 方程
    res = q*qx - a*qxx # 粘性burgers 方程
    return res

L = 2
N = 32
x = np.zeros((N))
k = np.zeros((N))
exact = np.zeros((3,N))
fourier = np.zeros((3,N))
error = np.zeros((3,N))
for i in range(0,N):
    x[i] = L/N*(i - N/2)
    if i < N/2:
        k[i] = 2*np.pi/L*i
    else:
        k[i] = 2*np.pi/L*(i - N)
u = np.exp(np.sin(np.pi*x))
ut = FFT(u)
exact[0,:] = np.pi*np.cos(np.pi*x)*u#解析一阶导
exact[1,:] = np.pi*np.pi*(np.cos(np.pi*x)**2 - np.sin(np.pi*x))*u
exact[2,:] = np.pi**3 * np.cos(np.pi*x)*(np.cos(np.pi*x)**2 - 3*np.sin(np.pi*x)-1)*u
fourier[0,:] = IFFT(complex(0,1)*k*ut).real#谱方法一阶导
fourier[1,:] = IFFT((complex(0,1)*k)**2*ut).real
fourier[2,:] = IFFT((complex(0,1)*k)**3*ut).real
error[:,:] = exact[:,:] - fourier[:,:]#误差