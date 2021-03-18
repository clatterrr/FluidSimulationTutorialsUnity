import numpy as np
"""
有限体积法解二维一元函数，教程及程序作者：光影帽子
https://zhuanlan.zhihu.com/p/358033368
解的偏微分方程为 dU/dt + adU/dx = 0
U为未知变量，a为常数
"""

#可以修改的初始值
a = 2#即上式中的a
num = 4#行数和列数

#接下来交给硅基生物完成，你们这些碳基生物快一边凉快去ε=ε=ε=(~￣▽￣)~
nmax = num*num#总顶点数
kmat = np.zeros((nmax,nmax))
kinv = np.zeros((nmax,nmax))
bmat = np.zeros((nmax))
res = np.zeros((num,num))
for i in range(0,num-1):
    for j in range(0,num-1):
        a0 = j*num + i
        a1 = a0 + 1
        a2 = a0 + num
        a3 = a2 + 1
        '''
        a2 --- a3
        |       |
        a0 --- a1
        
        '''
        kmat[a0,a0] = kmat[a0,a0] +  (-1/6) + a*(-1/6)
        kmat[a0,a1] = kmat[a0,a1] +  1/6 + a*(-1/12)
        kmat[a0,a2] = kmat[a0,a2] +  (-1/12) + a*(1/6)
        kmat[a0,a3] = kmat[a0,a3] +  (1/12) + a*(1/12)
        
        kmat[a1,a0] = kmat[a1,a0] +  (-1/6) + a*(-1/12)
        kmat[a1,a1] = kmat[a1,a1] +  1/6 + a*(-1/6)
        kmat[a1,a2] = kmat[a1,a2] +  (-1/12) + a*(1/12)
        kmat[a1,a3] = kmat[a1,a3] +  (1/12) + a*(1/6)
        
        kmat[a2,a0] = kmat[a2,a0] +  (-1/12) + a*(-1/6)
        kmat[a2,a1] = kmat[a2,a1] +  1/12 + a*(-1/12)
        kmat[a2,a2] = kmat[a2,a2] +  (-1/6) + a*(1/6)
        kmat[a2,a3] = kmat[a2,a3] +  (1/6) + a*(1/12)
        
        kmat[a3,a0] = kmat[a3,a0] + (-1/12) + a*(-1/12)
        kmat[a3,a1] = kmat[a3,a1] + 1/12 + a*(-1/6)
        kmat[a3,a2] = kmat[a3,a2] + (-1/6) + a*(1/12)
        kmat[a3,a3] = kmat[a3,a3] + (1/6) + a*(1/6)
        
for i in range(0,nmax):
    kinv[i,i] = 1
        
#应用初始条件修改矩阵
for i in range(0,num):
    kmat[i,:] = 0
    kmat[i,i] = 1
    bmat[i] = i
kmat2 = kmat.copy()
temprow = np.zeros((nmax))
#高斯消元法求逆矩阵
for i in range(0,nmax):
    if abs(kmat[i,i]) < 1e-10:
        for j in range(i+1,nmax):
            if kmat[j,i] != 0:
                temprow[:] = kmat[i,:]
                kmat[i,:] = kmat[j,:]
                kmat[j,:] = temprow[:]
                temprow[:] = kinv[i,:]
                kinv[i,:] = kinv[j,:]
                kinv[j,:] = temprow[:]
                break
    temp = kmat[i,i]
    for j in range(0,nmax):
        kmat[i,j] = kmat[i,j]/temp
        kinv[i,j] = kinv[i,j]/temp
    for k in range(0,nmax):
        if i == k:
            continue
        if kmat[k,i] == 0:
            continue
        mul = kmat[k,i]
        for j in range(0,nmax):
            kmat[k,j] -= mul*kmat[i,j]
            kinv[k,j] -= mul*kinv[i,j]
# 求解最后结果
for i in range(0,nmax):
    temp = 0
    for j in range(0,nmax):
        temp += kinv[i,j]*bmat[j]
    res[(int)(i/num),(int)(i%num)] = temp
