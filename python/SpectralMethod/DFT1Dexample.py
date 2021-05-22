import numpy as np
nmax = 4
arr = np.zeros((nmax),dtype=complex)
for i in range(0,nmax):
    arr[i] = i + 1
libres = np.fft.fft(arr)#调用库函数的结果,最终和dftres一样
idftres = np.zeros((nmax),dtype=complex)#正变换后加反变换的结果，应该和原始arr一样
dftres = np.zeros((nmax),dtype=complex)
# DFT
for x0 in range(0,nmax):
    for x1 in range(0,nmax):
        dftres[x0] += arr[x1]*np.exp(-2j*np.pi*(x0 * x1 / nmax))
# IDFT    
for x0 in range(0,nmax):
    for x1 in range(0,nmax):
        idftres[x0] += dftres[x1]*np.exp(-2j*np.pi*(x0 * x1 / nmax))/nmax