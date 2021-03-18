import numpy as np
"""
有限体积法解一维一元二次函数，教程及程序作者：光影帽子
https://zhuanlan.zhihu.com/p/358033368
解的偏微分方程为 d^2U/dx^2 - a = 0
U为未知变量，a为常数

"""

#可以修改的初始值，需要两个初始值
idx0 = 0#初始条件
value0 = 1
idx1 = 1
value1 = 2
a = 2#即上式dU/dx - a = 0中的a
num = 3#需要num行num列，共可算出num*2 + 1个解

#接下来交给硅基生物完成，你们这些碳基生物快一边凉快去ε=ε=ε=(~￣▽￣)~
nmax = num*2 + 1#需要求出的解的数量
kmat = np.zeros((nmax,nmax))
kinv = np.zeros((nmax,nmax))
bmat = np.zeros((nmax))
res = np.zeros((nmax))
for i in range(0,nmax-2):
    kmat[i:i+3,i:i+3] += [[1/3,-2/3,1/3],[0,0,0],[0,0,0]]
    bmat[i:i+3] += [a/3,0,0]
for i in range(0,nmax):
    kinv[i,i] = 1
#应用初始条件修改矩阵
kmat[nmax-2,:] = 0
kmat[nmax-2,idx0] = 1
bmat[nmax-2] = value0
kmat[nmax-1,:] = 0
kmat[nmax-1,idx1] = 1
bmat[nmax-1] = value1
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
#求解最后结果
for i in range(0,nmax):
    temp = 0
    for j in range(0,nmax):
        temp += kinv[i,j]*bmat[j]
    res[i] = temp
