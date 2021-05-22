import numpy as np
nmax = 4
mmax = 4
arr = np.zeros((nmax,mmax),dtype=complex)
for i in range(0,nmax):
    for j in range(0,mmax):
        arr[i,j] = i + j + 1
libres = np.fft.fft2(arr)#调用库函数的结果,最终和dftres一样
idftres = np.zeros((nmax,mmax),dtype=complex)#正变换后加反变换的结果，应该和原始arr一样
dftres = np.zeros((nmax,mmax),dtype=complex)
# 二维离散傅里叶变换
for y0 in range(0,mmax):
    for x0 in range(0,nmax):
        for y1 in range(0,mmax):
            for x1 in range(0,nmax):
                dftres[x0,y0] += arr[x1,y1]*np.exp(-2j*np.pi*(x0 * x1 / nmax + y0 * y1 / mmax))
# 二维离散傅里叶逆变换
for y0 in range(0,mmax):
    for x0 in range(0,nmax):
        for y1 in range(0,mmax):
            for x1 in range(0,nmax):
                idftres[x0,y0] += dftres[x1,y1]*np.exp(2j*np.pi*(x0 * x1 / nmax + y0 * y1 / mmax))/nmax/mmax