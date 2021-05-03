import numpy as np
'''
Galerkin法求解非线性方程du/dt + udu/dx
设u = c0x + c1t + c2
'''

c = np.zeros((3))
c[:] = 1
d = np.zeros((3))
d[0] = c[0]*c[0]
d[1] = c[0]*c[1]
d[2] = c[0]*c[2] + c[1]
# 接下来用有限元Galerkin法解
def phi(n,t,x):
    if n == 0:
        return (1 - t)*(1 - x)
    elif n == 1:
        return (1 - t)*x
    elif n == 2:
        return t*(1 - x)
    elif n == 3:
        return t*x
    
def simpon(p0,p1,p2,h):
    return h/3*(p0 + 4*p1 + p2)
def bool5(point0,point1,point2,point3,point4,h):#五阶精度数值积分
    return 2*h/45*(7*point0 + 32*point1 + 
                   12*point2 + 32*point3 + 7*point4)
num = 4 # 权函数的数量
rhs = np.zeros((3,num))
inte = np.zeros((5))
inte1 = np.zeros((5))
Fcoeff = np.zeros((num,num,num)) # 非线性方程组系数
Fdt = np.zeros((num,num))
for i in range(0,num):
    for tr in range(0,3):
        t0 = tr / 2
        for xr in range(0,3):
            x0 = xr / 2
            inte[xr] = x0*phi(i,t0,x0)
        inte1[tr] = simpon(inte[0], inte[1], inte[2], 0.5)
    rhs[0,i] = simpon(inte1[0], inte1[1], inte1[2], 0.5)
    for tr in range(0,3):
        t0 = tr / 2
        for xr in range(0,3):
            x0 = xr / 2
            inte[xr] = t0*phi(i,t0,x0)
        inte1[tr] = simpon(inte[0], inte[1], inte[2], 0.5)
    rhs[1,i] = simpon(inte1[0], inte1[1], inte1[2], 0.5)
    for tr in range(0,3):
        t0 = tr / 2
        for xr in range(0,3):
            x0 = xr / 2
            inte[xr] = phi(i,t0,x0)
        inte1[tr] = simpon(inte[0], inte[1], inte[2], 0.5)
    rhs[2,i] = simpon(inte1[0], inte1[1], inte1[2], 0.5)

for i in range(0,num):
    for j in range(0,num):
        for tr in range(0,3):
            t0 = tr / 2
            for xr in range(0,3):
                x0 = xr / 2
                inte[xr] = phi(i,t0,x0)*(phi(j,t0+1,x0) - phi(j,t0-1,x0))/2
            inte1[tr] = simpon(inte[0], inte[1], inte[2], 0.5)
        Fdt[i,j] = simpon(inte1[0], inte1[1], inte1[2], 0.5)

for k in range(0,num):
    for i in range(0,num):
        for j in range(0,num):
            for tr in range(0,5):
                t0 = tr / 4
                for xr in range(0,5):
                    x0 = xr / 4
                    inte[xr] = phi(k,t0,x0)*phi(i,t0,x0)*(phi(j,t0,x0+1) - phi(j,t0,x0-1))/2
                inte1[tr] = bool5(inte[0],inte[1],inte[2],inte[3],inte[4],0.25)
            Fcoeff[k,i,j] = bool5(inte1[0],inte1[1],inte1[2],inte1[3],inte1[4],0.25)

itemax = 100
tmax = 3 # 牛顿迭代的时间
nmax = 3 # 待求解变量的个数
total = tmax * nmax
res = np.zeros((itemax,total)) # 待求解的变量
cons = np.zeros((total))
cons[0] = cons[1] = 1
# 初始猜的值不能全是零，否则之后的雅可比矩阵是奇异的，无法求逆
for i in range(0,tmax):
    for j in range(0,nmax):
        x0 = j 
        t0 = i 
        idx = i * nmax + j
        if cons[idx] == 1:
            res[:,idx] = c[0]*x0 + c[1]*t0 + c[2]
        else:
            res[0,idx] = j**3 + 5

Fmat = np.zeros((total))
Fmat2 = np.zeros((total))
Fmat3 = np.zeros((total))
Jmat = np.zeros((total,total))
# 开始牛顿迭代法
idxk = np.zeros((num))# 2 个权函数一次可求解两个变量
test = np.zeros((total))
for t in range(0,itemax-1):
    Fmat[:] = Fmat2[:] = Fmat3[:] = 0
    Jmat[:,:] = 0
    for k in range(0,tmax*nmax):
        # 现在要求解k和k + 1处的变量
        colu = int(k / nmax)
        row = k % nmax
        idxk[0] = colu * nmax + row
        idxk[1] = idxk[0] + 1
        idxk[2] = idxk[0] + nmax
        idxk[3] = idxk[2] + 1
        if(colu >= tmax - 1)|(row >= nmax - 1):
            continue
        for i in range(0,num):
            idx = int(idxk[i])
            td0 = d[0]
            td1 = d[1]
            td2 = d[0]*row + d[1]*colu + d[2]
            Fmat[idx] = Fmat[idx]  - td0*rhs[0,i] - td1*rhs[1,i] - td2*rhs[2,i]
        for i in range(0,num):
            for j in range(0,num):
                i0 = int(idxk[i])
                j0 = int(idxk[j])
                Fmat3[i0] = Fmat3[i0] + Fdt[i,j]*res[t,j0]
                Jmat[i0,j0] += Fdt[i,j]
        for p in range(0,num):
            for i in range(0,num):
                for j in range(0,num):
                    k0 = int(idxk[p])# 求解的位置
                    i0 = int(idxk[i])# u1的索引
                    j0 = int(idxk[j])# u2的索引
                    Fmat2[k0] = Fmat2[k0] + Fcoeff[p,i,j]*res[t,i0]*res[t,j0]
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
    for i in range(0,total):
        summ = 0
        for j in range(0,total):
            summ += Jmat[i,j]*res[t,j]
        test[i] = summ
        
    # J du = F
    # 这里du_0 应是零，也就是u_0一直作为设定好的初始条件一直不变
    # 当算法最终收敛的时候，Fmat + Fmat2为零，这也是调试时重要的条件
    for i in range(0,total):
        if cons[i] == 1:
            Jmat[i,:] = 0
            Jmat[i,i] = 1
            Fmat[i] = Fmat2[i] = Fmat3[i] =  0
    
    
    
    Jinv = np.linalg.inv(Jmat)
    for i in range(0,total):
        summ = 0
        for j in range(0,total):
            summ += Jinv[i,j]*(- Fmat[j] - Fmat2[j] - Fmat3[j])
        res[t+1,i] = res[t,i] + summ   