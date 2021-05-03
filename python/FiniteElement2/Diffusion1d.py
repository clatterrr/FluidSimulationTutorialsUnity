import numpy as np
'''
Galerkin法解压力泊松方程

'''

def pf(n,t,x):
    if n == 0:
        return (t*t/2 - t/2)*(x*x/2 - x/2)
    elif n == 1:
        return (t*t/2 - t/2)*(1 - x*x)
    elif n == 2:
        return (t*t/2 - t/2)*(x*x/2 + x/2)
    elif n == 3:
        return (1 - t*t)*(x*x/2 - x/2)
    elif n == 4:
        return (1 - t*t)*(1 - x*x)
    elif n == 5:
        return (1 - t*t)*(x*x/2 + x/2)
    elif n == 6:
        return (t*t/2 + t/2)*(x*x/2 - x/2)
    elif n == 7:
        return (t*t/2 + t/2)*(1 - x*x)
    elif n == 8:
        return (t*t/2 + t/2)*(x*x/2 + x/2)
        
def bool5(point0,point1,point2,point3,point4,h):#五阶精度
    return 2*h/45*(7*point0 + 32*point1 + 
                   12*point2 + 32*point3 + 7*point4)


nmax = 100 # 时间
mmax = 8 # x轴
total = nmax*mmax
num = 9#每个正方形有多少个顶点，x轴两个，t轴两个，共四个
Dx2mat = np.zeros((num,num))#用于组装的矩阵之扩散矩阵
Tmat = np.zeros((num,num))
kmat = np.zeros((total,total))
bmat = np.zeros((total))
btmat = np.zeros((total))
Inte0 = np.zeros((5))
Inte1 = np.zeros((5))
res = np.zeros((total))
U = np.zeros((nmax,mmax))
        
    
for i in range(0,num):
    for j in range(0,num):
        for tr in range(0,5):
            t0 = tr/2 - 1
            for xr in range(0,5):
                x0 = xr/2 - 1
                #注意i和j
                Inte0[xr] = -(pf(i,t0,x0+1) - pf(i,t0,x0-1))/2*(pf(j,t0,x0+1) - pf(j,t0,x0-1))/2
            Inte1[tr] = bool5(Inte0[0],Inte0[1],Inte0[2],Inte0[3],Inte0[4],0.5)
            Inte1[tr] = Inte1[tr] + pf(i,t0,1)*(pf(j,t0,2) - pf(j,t0,0))/2 - pf(i,t0,-1)*(pf(j,t0,0) - pf(j,t0,-2))/2 
        Dx2mat[i,j] = bool5(Inte1[0],Inte1[1],Inte1[2],Inte1[3],Inte1[4],0.5)

for i in range(0,num):
    for j in range(0,num):
        for xr in range(0,5):
            x0 = xr/2 - 1
            for tr in range(0,5):
                t0 = tr/2 - 1
                Inte0[tr] = pf(i,t0,x0)*(pf(j,t0+1,x0) - pf(j,t0-1,x0))/2
            Inte1[xr] = bool5(Inte0[0],Inte0[1],Inte0[2],Inte0[3],Inte0[4],0.5)
        Tmat[i,j] = bool5(Inte1[0],Inte1[1],Inte1[2],Inte1[3],Inte1[4],0.5)

idxk = np.zeros((9))
for k in range(0,total):
    # 组装矩阵
    colu = int(k/mmax) # t轴
    row = int(k%mmax) # x轴
    for t0 in range(0,3):
        for x0 in range(0,3):
            idx = (colu + t0)*mmax + row + x0
            idxx = int(t0*3 + x0)
            idxk[idxx] = idx
    if(colu >= (nmax - 2))|(row >= (mmax - 2)):
        continue
    for i in range(0,num):
        for j in range(0,num):
            i0 = int(idxk[i])
            j0 = int(idxk[j])#u1的索引
            kmat[i0,j0] = kmat[i0,j0] + Tmat[i,j] - Dx2mat[i,j] 

test = np.zeros((total))
for i in range(0,mmax):
    for j in range(0,nmax):
        idx = j *mmax + i
        test[idx] = i*i + 2*j
    
for i in range(0,nmax):
    summ = 0
    for j in range(0,nmax):
        summ += kmat[i,j]*test[j]
    btmat[i] = summ  
# 设置初始条件
for i in range(0,mmax):
    idx = i
    kmat[idx,:] = 0
    kmat[idx,idx] = 1
    bmat[idx] = i*4
for i in range(0,nmax):
    idx = i*mmax
    kmat[idx,:] = 0
    kmat[idx,idx] = 1
    bmat[idx] = 0
    idx = i*mmax + mmax - 1
    kmat[idx,:] = 0
    kmat[idx,idx] = 1
    bmat[idx] = mmax

kinv = np.linalg.inv(kmat)
for i in range(0,total):
    summ = 0
    for j in range(0,total):
        summ += kinv[i,j]*bmat[j]
    res[i] = summ  
for t0 in range(0,nmax):
    for x0 in range(0,mmax):
            U[t0,x0] = res[t0*mmax + x0]
