import numpy as np
import time
level = 6
nmax = 2**level
total = nmax * nmax
arr = np.zeros((nmax,nmax),dtype=complex)
for i in range(0,nmax):
    for j in range(0,nmax):
        arr[i,j] = i + j + 1
start =time.process_time()
libres = np.fft.fft2(arr)#调用库函数的结果,最终和dftres一样
end =time.process_time()
print('库函数运行时间: %s 秒'%(end-start))
idftres = np.zeros((nmax,nmax),dtype=complex)#正变换后加反变换的结果，应该和原始arr一样
dftres = np.zeros((nmax,nmax),dtype=complex)
start =time.process_time()
# 二维离散快速傅里叶变换
dftres = arr.copy()
# 先对x轴做一遍，此时y轴等于零就行了
for level0 in range(0,level):
    maxi = 2**level0
    xdis = nmax//maxi
    for i in range(0,maxi):
        x0 = int(i % nmax)
        y0 = int(i / nmax)
        even =  dftres[x0,y0] # 蝴蝶变换中，偶数项已经计算好，因此之间复制就行
        odd = 0
        # 蝴蝶变换
        for j in range(0,nmax,xdis):
            idx = j + xdis//2
            x1 = int(idx % nmax)
            y1 = int(idx / nmax)
            odd += arr[x1,y1] * np.exp(-2j*np.pi*x0*x1/nmax)#奇数项需要重新计算
        # 消去引理
        dftres[x0,y0] = even + odd
        dftres[x0+maxi,y0] = even - odd

# 然后再对y轴做一遍
for level0 in range(0,level):
    maxi = 2**level0
    xdis = total//maxi # 方块之间一下跳跃的距离
    maxii = maxi * nmax # 目前已有的行数
    for i in range(0,maxii):
        x0 = int(i % nmax)
        y0 = int(i / nmax)
        even =  dftres[x0,y0] # 蝴蝶变换中，偶数项已经计算好，因此之间复制就行
        odd = 0
        # 蝴蝶变换
        for j0 in range(0,total,xdis):
            for j1 in range(0,nmax):
                idx = j0 + xdis//2 + j1
                x1 = int(idx % nmax)
                y1 = int(idx / nmax)
                odd += arr[x1,y1] * np.exp(-2j*np.pi*(x0*x1 + y0*y1)/nmax)#奇数项需要重新计算
        # 消去引理
        dftres[x0,y0] = even + odd
        dftres[x0,y0+maxi] = even - odd
end =time.process_time()
print('自己的函数运行时间: %s 秒'%(end-start))
# 二维离散快速傅里叶逆变换
idftres = dftres.copy()/total
# 先对x轴做一遍
for level0 in range(0,level):
    maxi = 2**level0
    xdis = nmax//maxi
    for i in range(0,maxi):
        x0 = int(i % nmax)
        y0 = int(i / nmax)
        even =  idftres[x0,y0] # 蝴蝶变换中，偶数项已经计算好，因此之间复制就行
        odd = 0
        # 蝴蝶变换
        for j in range(0,nmax,xdis):
            idx = j + xdis//2
            x1 = int(idx % nmax)
            y1 = int(idx / nmax)
            odd += dftres[x1,y1] * np.exp(2j*np.pi*x0*x1/nmax)/total#奇数项需要重新计算
        # 消去引理
        idftres[x0,y0] = even + odd
        idftres[x0+maxi,y0] = even - odd

# 然后再对y轴做一遍
for level0 in range(0,level):
    maxi = 2**level0
    xdis = total//maxi # 方块之间一下跳跃的距离
    maxii = maxi * nmax # 目前已有的行数
    for i in range(0,maxii):
        x0 = int(i % nmax)
        y0 = int(i / nmax)
        even =  idftres[x0,y0] # 蝴蝶变换中，偶数项已经计算好，因此之间复制就行
        odd = 0
        # 蝴蝶变换
        for j0 in range(0,total,xdis):
            for j1 in range(0,nmax):
                idx = j0 + xdis//2 + j1
                x1 = int(idx % nmax)
                y1 = int(idx / nmax)
                odd += dftres[x1,y1] * np.exp(2j*np.pi*(x0*x1 + y0*y1)/nmax)/total#奇数项需要重新计算
        # 消去引理
        idftres[x0,y0] = even + odd
        idftres[x0,y0+maxi] = even - odd    


