import numpy as np
'''
Galerkin法解压力泊松方程

'''

def pf(n,t,x,y):
    if n == 0:
        return (t*t/2 - t/2)*(x*x/2 - x/2)*(y*y/2 - y/2)
    elif n == 1:
        return (t*t/2 - t/2)*(x*x/2 - x/2)*(1 - y*y)
    elif n == 2:
        return (t*t/2 - t/2)*(x*x/2 - x/2)*(y*y/2 + y/2)
    elif n == 3:
        return (t*t/2 - t/2)*(1 - x*x)*(y*y/2 - y/2)
    elif n == 4:
        return (t*t/2 - t/2)*(1 - x*x)*(1 - y*y)
    elif n == 5:
        return (t*t/2 - t/2)*(1 - x*x)*(y*y/2 + y/2)
    elif n == 6:
        return (t*t/2 - t/2)*(x*x/2 + x/2)*(y*y/2 - y/2)
    elif n == 7:
        return (t*t/2 - t/2)*(x*x/2 + x/2)*(1 - y*y)
    elif n == 8:
        return (t*t/2 - t/2)*(x*x/2 + x/2)*(y*y/2 + y/2)
    elif n == 9:
        return (1 - t*t)*(x*x/2 - x/2)*(y*y/2 - y/2)
    elif n == 10:
        return (1 - t*t)*(x*x/2 - x/2)*(1 - y*y)
    elif n == 11:
        return (1 - t*t)*(x*x/2 - x/2)*(y*y/2 + y/2)
    elif n == 12:
        return (1 - t*t)*(1 - x*x)*(y*y/2 - y/2)
    elif n == 13:
        return (1 - t*t)*(1 - x*x)*(1 - y*y)
    elif n == 14:
        return (1 - t*t)*(1 - x*x)*(y*y/2 + y/2)
    elif n == 15:
        return (1 - t*t)*(x*x/2 + x/2)*(y*y/2 - y/2)
    elif n == 16:
        return (1 - t*t)*(x*x/2 + x/2)*(1 - y*y)
    elif n == 17:
        return (1 - t*t)*(x*x/2 + x/2)*(y*y/2 + y/2)
    elif n == 18:
        return (t*t/2 + t/2)*(x*x/2 - x/2)*(y*y/2 - y/2)
    elif n == 19:
        return (t*t/2 + t/2)*(x*x/2 - x/2)*(1 - y*y)
    elif n == 20:
        return (t*t/2 + t/2)*(x*x/2 - x/2)*(y*y/2 + y/2)
    elif n == 21:
        return (t*t/2 + t/2)*(1 - x*x)*(y*y/2 - y/2)
    elif n == 22:
        return (t*t/2 + t/2)*(1 - x*x)*(1 - y*y)
    elif n == 23:
        return (t*t/2 + t/2)*(1 - x*x)*(y*y/2 + y/2)
    elif n == 24:
        return (t*t/2 + t/2)*(x*x/2 + x/2)*(y*y/2 - y/2)
    elif n == 25:
        return (t*t/2 + t/2)*(x*x/2 + x/2)*(1 - y*y)
    elif n == 26:
        return (t*t/2 + t/2)*(x*x/2 + x/2)*(y*y/2 + y/2)
    
def bool5(point0,point1,point2,point3,point4,h):#五阶精度
    return 2*h/45*(7*point0 + 32*point1 + 
                   12*point2 + 32*point3 + 7*point4)

tmax = 100 # t轴
nmax = 4 # x轴
mmax = 4 # y轴
total = tmax*nmax*mmax
num = 27#每个正方形有多少个顶点，x轴两个，t轴两个，共四个
Dx2mat = np.zeros((num,num))#用于组装的矩阵之扩散矩阵
Dy2mat = np.zeros((num,num))
Tmat = np.zeros((num,num))
rhs = np.zeros((num))#等式右边项，用于组装
kmat = np.zeros((total,total))
bmat = np.zeros((total))
btmat = np.zeros((total))
Inte0 = np.zeros((5))
Inte1 = np.zeros((5))
Inte2 = np.zeros((5))
res = np.zeros((total))
U = np.zeros((tmax,nmax,mmax))

for i in range(0,num):
    for j in range(0,num):
        for tr in range(0,5):
            t0 = tr/2 - 1
            for xr in range(0,5):
                x0 = xr/2 - 1
                for yr in range(0,5):
                    y0 = yr/2 - 1
                    Inte0[yr] = -(pf(i,t0,x0,y0+1) - pf(i,t0,x0,y0-1))/2*(pf(j,t0,x0,y0+1) - pf(j,t0,x0,y0-1))/2
                Inte1[xr] = bool5(Inte0[0],Inte0[1],Inte0[2],Inte0[3],Inte0[4],0.5)
                Inte1[xr] = Inte1[xr] + pf(i,t0,x0,1)*(pf(j,t0,x0,2) - pf(j,t0,x0,0))/2 - pf(i,t0,x0,-1)*(pf(j,t0,x0,0) - pf(j,t0,x0,-2))/2 
            Inte2[tr] = bool5(Inte1[0],Inte1[1],Inte1[2],Inte1[3],Inte1[4],0.5)
        Dy2mat[i,j] = bool5(Inte2[0],Inte2[1],Inte2[2],Inte2[3],Inte2[4],0.5)
        
    
for i in range(0,num):
    for j in range(0,num):
        for tr in range(0,5):
            t0 = tr/2 - 1
            for yr in range(0,5):
                y0 = yr/2 - 1
                for xr in range(0,5):
                    x0 = xr/2 - 1
                    Inte0[xr] = -(pf(i,t0,x0+1,y0) - pf(i,t0,x0-1,y0))/2*(pf(j,t0,x0+1,y0) - pf(j,t0,x0-1,y0))/2
                Inte1[yr] = bool5(Inte0[0],Inte0[1],Inte0[2],Inte0[3],Inte0[4],0.5)
                Inte1[yr] = Inte1[yr] + pf(i,t0,1,y0)*(pf(j,t0,2,y0) - pf(j,t0,0,y0))/2 - pf(i,t0,-1,y0)*(pf(j,t0,0,y0) - pf(j,t0,-2,y0))/2 
            Inte2[tr] = bool5(Inte1[0],Inte1[1],Inte1[2],Inte1[3],Inte1[4],0.5)
        Dx2mat[i,j] = bool5(Inte2[0],Inte2[1],Inte2[2],Inte2[3],Inte2[4],0.5)

for i in range(0,num):
    for j in range(0,num):
        for tr in range(0,5):
            t0 = tr/2 - 1
            for yr in range(0,5):
                y0 = yr/2 - 1
                for xr in range(0,5):
                    x0 = xr/2 - 1
                    Inte0[xr] = pf(i,t0,x0,y0)*(pf(j,t0+1,x0,y0) - pf(j,t0-1,x0,y0))/2
                Inte1[yr] = bool5(Inte0[0],Inte0[1],Inte0[2],Inte0[3],Inte0[4],0.5)
            Inte2[tr] = bool5(Inte1[0],Inte1[1],Inte1[2],Inte1[3],Inte1[4],0.5)
        Tmat[i,j] = bool5(Inte2[0],Inte2[1],Inte2[2],Inte2[3],Inte2[4],0.5)

#装配刚度矩阵K
idxk = np.zeros((27))
for k in range(0,total):
    layer = int(k/nmax/mmax) # 时间层
    colu = int(k%(nmax*mmax)/mmax) # x轴
    row = int(k%mmax) # y轴
    for t0 in range(0,3):
        for x0 in range(0,3):
            for y0 in range(0,3):
                idx = (layer + t0)*nmax*mmax + (colu + x0)*mmax + row + y0
                idxx = int(t0*9 + x0*3 + y0)
                idxk[idxx] = idx
        
    if(colu >= nmax - 2)|(row >= mmax - 2)|(layer >= tmax - 2):
        continue
    for i in range(0,num):
        for j in range(0,num):
            i0 = int(idxk[i])
            j0 = int(idxk[j])#u1的索引
            kmat[i0,j0] = kmat[i0,j0] + Tmat[i,j] - (Dx2mat[i,j] + Dy2mat[i,j])

#给初始条件，每个时间层的第一个和最后一个网格温度不变
for t0 in range(0,tmax):
    idx = (t0 + 1)*nmax*mmax - 1
    kmat[idx,:] = 0
    kmat[idx,idx] = 1
    bmat[idx] = 8
    
for t0 in range(0,tmax):
    idx = t0*nmax*mmax
    kmat[idx,:] = 0
    kmat[idx,idx] = 1
    bmat[idx] = 0
    
#没有非线性和牛顿迭代法的代码写起来太舒服了
kinv = np.linalg.inv(kmat)

#得出结果,res是一列结果，U是重新按照txy轴排列的结果
for i in range(0,total):
    summ = 0
    for j in range(0,total):
        summ += kinv[i,j]*bmat[j]
    res[i] = summ  
for t0 in range(0,tmax):
    for x0 in range(0,nmax):
        for y0 in range(0,mmax):
            U[t0,x0,y0] = res[t0*nmax*mmax + x0*mmax + y0]