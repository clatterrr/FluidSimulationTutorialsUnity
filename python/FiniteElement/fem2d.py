import numpy as np
"""
有限体积法解二维一元函数，教程及程序作者：光影帽子
https://zhuanlan.zhihu.com/p/358033368
解的偏微分方程为 dU/dx + adU/dy + b= 0
U为未知变量，ab为常数
例如对于 U = x + y + 1 那么a = 1，b = -2，或者a = 2且b = -3
为了正确求解逆矩阵，需要三个初始条件
"""

#可以修改的初始值，至少需要三个不同初始条件才能正确求解逆矩阵
idx0 = 0
value0 = 1
idx1 = 1
value1 = 2
idx2 = 2
value2 = 3
a = 2#即上式dU/dx + adU/dy + b= 0中的a
b = -3#即上式dU/dx + adU/dy + b= 0中的b
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
        
        bmat[a0] = bmat[a0] - b*1/4
        bmat[a1] = bmat[a1] - b*1/4
        bmat[a2] = bmat[a2] - b*1/4
        bmat[a3] = bmat[a3] - b*1/4
for i in range(0,nmax):
    kinv[i,i] = 1
        
#应用初始条件修改矩阵
kmat[idx0,:] = 0
kmat[idx0,idx0] = 1
bmat[idx0] = value0
kmat[idx1,:] = 0
kmat[idx1,idx1] = 1
bmat[idx1] = value1
kmat[idx2,:] = 0
kmat[idx2,idx2] = 1
bmat[idx2] = value2
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
