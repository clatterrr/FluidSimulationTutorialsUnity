import numpy as np
'''
Galerkin法求解非线性方程udu/dx = d1x^3 + d2x^2 + d3x + d4
设解为u = c1x^2 + c2x + c3
'''

c1 = 3
c2 = 5
c3 = 2
d1 = 2 * c1 * c1
d2 = 3 * c1 * c2
d3 = 2 * c1 * c3 + c2 * c2
d4 = c2 * c3
# 接下来用有限元Galerkin法解
def phi(n,x):
    if n == 0:
        return x*x/2 - x/2
    elif n == 1:
        return 1 - x*x
    elif n == 2:
        return x*x/2 + x/2
def simpon(p0,p1,p2,h):
    return h/3*(p0 + 4*p1 + p2)
def Integrate7(po0,po1,po2,po3,po4,po5,po6,h):
    return h*(po0*41/140 + po1*54/35 + po2*27/140
             + po3*68/35 + po4*27/140 + po5*54/35 + po6*41/140)
num = 3 # 权函数的数量
rhs = np.zeros((4,3))
inte = np.zeros((7))
Fcoeff = np.zeros((num,num,num)) # 非线性方程组系数
for i in range(0,num):
    for xr in range(0,7):
        x0 = xr / 3 - 1
        inte[xr] = x0*x0*x0*phi(i,x0)
    rhs[0,i] = Integrate7(inte[0],inte[1],inte[2],inte[3],inte[4],inte[5],inte[6],1/3)
    for xr in range(0,7):
        x0 = xr / 3 - 1
        inte[xr] = x0*x0*phi(i,x0)
    rhs[1,i] = Integrate7(inte[0],inte[1],inte[2],inte[3],inte[4],inte[5],inte[6],1/3)
    for xr in range(0,7):
        x0 = xr / 3 - 1
        inte[xr] = x0*phi(i,x0)
    rhs[2,i] = Integrate7(inte[0],inte[1],inte[2],inte[3],inte[4],inte[5],inte[6],1/3) 
    for xr in range(0,3):
        x0 = xr - 1
        inte[xr] = phi(i,x0)
    rhs[3,i] = simpon(inte[0],inte[1],inte[2],1)

for k in range(0,num):
    for i in range(0,num):
        for j in range(0,num):
            for xr in range(0,7):
                x0 = xr / 3 - 1
                inte[xr] = phi(k,x0)*phi(i,x0)*(phi(j,x0+1) - phi(j,x0-1))/2
            Fcoeff[k,i,j] = Integrate7(inte[0],inte[1],inte[2],inte[3],inte[4],inte[5],inte[6],1/3)

tmax = 100 # 牛顿迭代的时间
nmax = 8 # 待求解变量的个数
res = np.zeros((tmax,nmax)) # 待求解的变量
# 初始猜的值不能全是零，否则之后的雅可比矩阵是奇异的，无法求逆
for i in range(0,nmax):
    res[0,i] = c1*j*j + c2*j + c3
    res[0,i] = i*i*i + 5
res[:,0] = c1*(-1)*(-1) + c2*(-1) + c3 # 初始条件
Fmat = np.zeros((nmax))
Fmat2 = np.zeros((nmax))
Jmat = np.zeros((nmax,nmax))
# 开始牛顿迭代法
idxk = np.zeros((num))# 2 个权函数一次可求解两个变量
test = np.zeros((nmax))
for t in range(0,tmax-1):
    Fmat[:] = Fmat2[:] =0
    Jmat[:,:] = 0
    for k in range(0,nmax - 2):
        # 现在要求解k和k + 1处的变量
        idxk[0] = k 
        idxk[1] = k + 1
        idxk[2] = k + 2
        for i in range(0,num):
            idx = int(idxk[i])
            '''
            d1(x + k)^3 + d2(x+k)^2  + d3(x + k) + d4
            
            '''
            tempd1 = d1
            tempd2 = d1*3*k + d2
            tempd3 = d1*3*k*k + d2*2*k + d3
            tempd4 = d1*k*k*k + d2*k*k + d3*k + d4
            Fmat[idx] = Fmat[idx]  - tempd1*rhs[0,i] - tempd2*rhs[1,i] - tempd3*rhs[2,i] - tempd4*rhs[3,i]
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
    for i in range(0,nmax):
        summ = 0
        for j in range(0,nmax):
            summ += Jmat[i,j]*res[t,j]
        test[i] = summ
        
    # J du = F
    # 这里du_0 应是零，也就是u_0一直作为设定好的初始条件一直不变
    # 当算法最终收敛的时候，Fmat + Fmat2为零，这也是调试时重要的条件
    Jmat[0,:] = 0
    Jmat[0,0] = 1
    Fmat[0] = Fmat2[0] = 0
    
    
    Jinv = np.linalg.inv(Jmat)
    for i in range(0,nmax):
        summ = 0
        for j in range(0,nmax):
            summ += Jinv[i,j]*(- Fmat[j] - Fmat2[j])
        res[t+1,i] = res[t,i] + summ   