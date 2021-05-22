import numpy as np
level = 3 # 需要将矩阵扩大3次，1->2,2->4,4->8
nmax = 2**level
arr = np.zeros((nmax),dtype=complex)
arr2 = np.zeros((nmax),dtype=complex)
for i in range(0,nmax):
    arr[i] = i + 1
libres = np.fft.fft(arr)#调用库函数的结果,最终和dftres一样
dftres = arr.copy()
# DFT
for level0 in range(0,level):
    # 现在算出了2**level0的矩阵，要计算2**(level0 + 1)的矩阵
    even = 0 # 这部分算得相同
    odd = 0 # 这部分算得不同
    maxi = 2**level0
    interval = nmax//maxi
    for i in range(0,maxi): 
        even =  dftres[i] # 蝴蝶变换中，偶数项上个level已经计算好，直接复制就行
        odd = 0
        # 蝴蝶变换
        for j in range(0,nmax,interval):
            idx = j + interval//2
            # 奇数项需要重新计算，也就是在上个level中计算好的两列正中间的那一列
            odd += arr[idx] * np.exp(-2j*np.pi*idx*i/nmax)
        # 消去引理
        dftres[i] = even + odd
        dftres[i + maxi] = even - odd

# IDFT
idftres = dftres.copy()/nmax
for level0 in range(0,level):
    maxi = 2**level0
    for i in range(0,maxi): 
        odd = 0
        for j in range(0,nmax,nmax//maxi):
            idx = j + nmax//maxi//2
            odd += dftres[idx] * np.exp(2j*np.pi*idx*i/nmax)/nmax
        idftres[i] += odd
        idftres[i + maxi] = idftres[i] - 2*odd