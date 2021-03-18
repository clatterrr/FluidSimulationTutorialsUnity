import numpy as np
"""
有限体积法解一维一元函数，教程及程序作者：光影帽子
https://zhuanlan.zhihu.com/p/358033368
解的偏微分方程为 dU/dx - a = 0
U为未知变量，a为常数

"""

#可以修改的初始值
idx = 1#初始条件，下面两个的意思是U[idx] = value，注意是从零开始
value = 0
a = 2#即上式dU/dx - a = 0中的a
nmax = 8#需要求出的解的数量

#接下来交给硅基生物完成，你们这些碳基生物快一边凉快去ε=ε=ε=(~￣▽￣)~
kmat = np.zeros((nmax,nmax))
kinv = np.zeros((nmax,nmax))
bmat = np.zeros((nmax))
res = np.zeros((nmax))
for i in range(0,nmax-1):
    kmat[i:i+2,i:i+2] += [[-0.5,0.5],[-0.5,0.5]]
    bmat[i:i+2] += [a*0.5,a*0.5]
    kinv[i,i] = 1
#应用初始条件修改矩阵
kmat[0,:] = 0
kmat[0,idx] = 1
bmat[0] = value
#原矩阵主对角线上大部分都是零，为了方便求解逆矩阵
#因此将原矩阵整体上移循环一位
kmat2 = np.zeros((nmax,nmax))
kmat = kmat.copy()
kmat2[0:nmax-1,:] = kmat[1:nmax,:]
kmat2[nmax-1,:] = kmat[0,:]
kmat = kmat2.copy()
temp = bmat[0]
for i in range(0,nmax-1):
    bmat[i] = bmat[i+1]
bmat[nmax-1] = temp
#消元法求解逆矩阵
for i in range(0,nmax):
    temp = kmat[i,i]
    for j in range(0,nmax):
        kmat[i,j] = kmat[i,j]/temp
        kinv[i,j] = kinv[i,j]/temp
    for k in range(0,nmax):
        if i == k:
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
