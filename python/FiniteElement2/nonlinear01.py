import numpy as np
'''
Galerkin法求解非线性方程udu/dx = d1x + d2
设解为u = c1x + c2
那么d1x + d2 = c1c1x + c1c2
'''

c1 = 2
c2 = 1
d1 = c1 * c1
d2 = c1 * c2
# 接下来用有限元Galerkin法解
def phi(n,x):
    if n == 0:
        return 1 - x
    else:
        return x
def simpon(p0,p1,p2,h):
    return h/3*(p0 + 4*p1 + p2)
num = 2 # 权函数的数量
rhs = np.zeros((2,2))
inte = np.zeros((3))
Fcoeff = np.zeros((num,num,num)) # 非线性方程组系数
for i in range(0,num):
    rhs[0,i] = (phi(i,0) + phi(i,1))/2 #梯形公式算数值积分 , 权函数本身的积分
    for xr in range(0,3):
        x0 = xr / 2
        inte[xr] = x0*phi(i,x0) # 权函数乘上x的积分
    rhs[1,i] = simpon(inte[0], inte[1], inte[2], 0.5)
for k in range(0,num):
    for i in range(0,num):
        for j in range(0,num):
            for xr in range(0,3):
                x0 = xr / 2
                inte[xr] = phi(k,x0)*phi(i,x0)*(phi(j,x0+1) - phi(j,x0-1))/2
            Fcoeff[k,i,j] = simpon(inte[0], inte[1], inte[2], 0.5)

tmax = 100 # 牛顿迭代的时间
nmax = 6 # 待求解变量的个数
res = np.zeros((tmax,nmax)) # 待求解的变量
# 初始猜的值不能全是零，否则之后的雅可比矩阵是奇异的，无法求逆
for i in range(0,nmax):
    res[0,i] = i + 1
Fmat = np.zeros((nmax))
Jmat = np.zeros((nmax,nmax))
# 开始牛顿迭代法
idxk = np.zeros((num))# 2 个权函数一次可求解两个变量
test = np.zeros((nmax))
for t in range(0,tmax-1):
    Fmat[:] = 0
    Jmat[:,:] = 0
    for k in range(0,nmax - 1):
        # 现在要求解k和k + 1处的变量
        idxk[0] = k 
        idxk[1] = k + 1
        for i in range(0,num):
            idx = int(idxk[i])
            '''
            d1(x + k) + d2 = d1x + (d2 + d1 * k)
            '''
            tempd1 = d1
            tempd2 = d2 + d1 * idxk[0]
            Fmat[idx] = Fmat[idx]  - tempd2*rhs[0,i] -  tempd1*rhs[1,i]
        for p in range(0,num):
            for i in range(0,num):
                for j in range(0,num):
                    k0 = int(idxk[p])# 求解的位置
                    i0 = int(idxk[i])# u1的索引
                    j0 = int(idxk[j])# u2的索引
                    Fmat[k0] = Fmat[k0] + Fcoeff[p,i,j]*res[t,i0]*res[t,j0]
        for i in range(0,num):
            for j in range(0,num):
                i0 = int(idxk[i])
                j0 = int(idxk[j])
                for p in range(0,num):
                    for q in range(0,num):
                        if p == j:
                            Jmat[i0,j0] += Fcoeff[i,p,q]*res[t,int(idxk[q])]
                        if q == j:
                            Jmat[i0,j0] += Fcoeff[i,p,q]*res[t,int(idxk[p])]
    for i in range(0,nmax):
        summ = 0
        for j in range(0,nmax):
            summ += Jmat[i,j]*res[t,j]
        test[i] = summ  
    Jinv = np.linalg.inv(Jmat)
    for i in range(0,nmax):
        summ = 0
        for j in range(0,nmax):
            summ += Jinv[i,j]*(- Fmat[j])
        res[t+1,i] = res[t,i] + summ   