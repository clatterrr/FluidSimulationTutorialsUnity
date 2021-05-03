import numpy as np
'''
Galerkin法求解非线性方程du/dt + udu/dx 
设u = c0x^2 + c1x + c2t^2 + c3t + c4
'''

c = np.zeros((5))
c[:] = 1
d = np.zeros((8))
d[0] = 2*c[0]*c[0] # x^3
d[1] = 3*c[0]*c[1] # x^2
d[2] = 2*c[0]*c[4] + c[1]*c[1] # x
d[3] = c[1]*c[2] # t^2
d[4] = c[1]*c[3] + 2*c[2] # t
d[5] = 2*c[0]*c[2] # xt^2
d[6] = 2*c[0]*c[3] # xt
d[7] = c[1]*c[4] + c[3] # 1
# 接下来用有限元Galerkin法解
def phi(n,t,x):
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
def simpon(p0,p1,p2,h):
    return h/3*(p0 + 4*p1 + p2)
def bool5(point0,point1,point2,point3,point4,h):#五阶精度数值积分
    return 2*h/45*(7*point0 + 32*point1 + 
                   12*point2 + 32*point3 + 7*point4)
def Integrate7(po0,po1,po2,po3,po4,po5,po6,h):
    return h*(po0*41/140 + po1*54/35 + po2*27/140
             + po3*68/35 + po4*27/140 + po5*54/35 + po6*41/140)
num = 9 # 权函数的数量
rhs = np.zeros((8,num))
inte = np.zeros((7))
inte1 = np.zeros((7))
Fcoeff = np.zeros((num,num,num)) # 非线性方程组系数
Tmat = np.zeros((num,num))
for i in range(0,num):
    for tr in range(0,7):
        t0 = tr / 3 - 1
        for xr in range(0,7):
            x0 = xr / 3 - 1
            inte[xr] = x0*x0*x0*phi(i,t0,x0)
        inte1[tr] = Integrate7(inte[0],inte[1],inte[2],inte[3],inte[4],inte[5],inte[6],1/3)
    rhs[0,i] = Integrate7(inte1[0],inte1[1],inte1[2],inte1[3],inte1[4],inte1[5],inte1[6],1/3)
    for tr in range(0,7):
        t0 = tr / 3 - 1
        for xr in range(0,7):
            x0 = xr / 3 - 1
            inte[xr] = x0*x0*phi(i,t0,x0)
        inte1[tr] = Integrate7(inte[0],inte[1],inte[2],inte[3],inte[4],inte[5],inte[6],1/3)
    rhs[1,i] = Integrate7(inte1[0],inte1[1],inte1[2],inte1[3],inte1[4],inte1[5],inte1[6],1/3)
    for tr in range(0,7):
        t0 = tr / 3 - 1
        for xr in range(0,7):
            x0 = xr / 3 - 1
            inte[xr] = x0*phi(i,t0,x0)
        inte1[tr] = Integrate7(inte[0],inte[1],inte[2],inte[3],inte[4],inte[5],inte[6],1/3)
    rhs[2,i] = Integrate7(inte1[0],inte1[1],inte1[2],inte1[3],inte1[4],inte1[5],inte1[6],1/3)
    for tr in range(0,7):
        t0 = tr / 3 - 1
        for xr in range(0,7):
            x0 = xr / 3 - 1
            inte[xr] = t0*t0*phi(i,t0,x0)
        inte1[tr] = Integrate7(inte[0],inte[1],inte[2],inte[3],inte[4],inte[5],inte[6],1/3)
    rhs[3,i] = Integrate7(inte1[0],inte1[1],inte1[2],inte1[3],inte1[4],inte1[5],inte1[6],1/3)
    for tr in range(0,7):
        t0 = tr / 3 - 1
        for xr in range(0,7):
            x0 = xr / 3 - 1
            inte[xr] = t0*phi(i,t0,x0)
        inte1[tr] = Integrate7(inte[0],inte[1],inte[2],inte[3],inte[4],inte[5],inte[6],1/3)
    rhs[4,i] = Integrate7(inte1[0],inte1[1],inte1[2],inte1[3],inte1[4],inte1[5],inte1[6],1/3)
    for tr in range(0,7):
        t0 = tr / 3 - 1
        for xr in range(0,7):
            x0 = xr / 3 - 1
            inte[xr] = x0*t0*t0*phi(i,t0,x0)
        inte1[tr] = Integrate7(inte[0],inte[1],inte[2],inte[3],inte[4],inte[5],inte[6],1/3)
    rhs[5,i] = Integrate7(inte1[0],inte1[1],inte1[2],inte1[3],inte1[4],inte1[5],inte1[6],1/3)
    for tr in range(0,7):
        t0 = tr / 3 - 1
        for xr in range(0,7):
            x0 = xr / 3 - 1
            inte[xr] = x0*t0*phi(i,t0,x0)
        inte1[tr] = Integrate7(inte[0],inte[1],inte[2],inte[3],inte[4],inte[5],inte[6],1/3)
    rhs[6,i] = Integrate7(inte1[0],inte1[1],inte1[2],inte1[3],inte1[4],inte1[5],inte1[6],1/3)
    for tr in range(0,5):
        t0 = tr / 2 - 1
        for xr in range(0,5):
            x0 = xr / 2 - 1
            inte[xr] = phi(i,t0,x0)
        inte1[tr] = bool5(inte[0],inte[1],inte[2],inte[3],inte[4],0.5)
    rhs[7,i] = bool5(inte1[0],inte1[1],inte1[2],inte1[3],inte1[4],0.5)

for i in range(0,num):
    for j in range(0,num):
        for tr in range(0,5):
            t0 = tr / 2 - 1
            for xr in range(0,5):
                x0 = xr / 2 - 1
                inte[xr] = phi(i,t0,x0)*(phi(j,t0+1,x0) - phi(j,t0-1,x0))/2
            inte1[tr] = bool5(inte[0], inte[1], inte[2],inte[3],inte[4],0.5)
        Tmat[i,j] = bool5(inte1[0],inte1[1],inte1[2],inte1[3],inte1[4],0.5)

for k in range(0,num):
    for i in range(0,num):
        for j in range(0,num):
            for tr in range(0,7):
                t0 = tr / 3 - 1
                for xr in range(0,7):
                    x0 = xr / 3 - 1
                    inte[xr] = phi(k,t0,x0)*phi(i,t0,x0)*(phi(j,t0,x0+1) - phi(j,t0,x0-1))/2
                inte1[tr] = Integrate7(inte[0],inte[1],inte[2],inte[3],inte[4],inte[5],inte[6],1/3)
            Fcoeff[k,i,j] = Integrate7(inte1[0],inte1[1],inte1[2],inte1[3],inte1[4],inte1[5],inte1[6],1/3)

itemax = 2
tmax = 4 # 牛顿迭代的时间
nmax = 4 # 
total = tmax*nmax
res = np.zeros((itemax,total)) # 待求解的变量
# 初始猜的值不能全是零，否则之后的雅可比矩阵是奇异的，无法求逆
for i in range(0,tmax):
    for j in range(0,nmax):
        x0 = j - 1
        t0 = i - 1
        idx = i * nmax + j
        res[0,idx] = c[0]*x0*x0 + c[1]*x0 + c[2]*t0*t0 + c[3]*t0 + c[4] + 1
Fmat = np.zeros((total))
Fmat2 = np.zeros((total))
Fmat3 = np.zeros((total))
Jmat = np.zeros((total,total))
# 开始牛顿迭代法
idxk = np.zeros((num))# 2 个权函数一次可求解两个变量
test = np.zeros((total))
for ite in range(0,itemax-1):
    Fmat[:] = Fmat2[:] = Fmat3[:] = 0
    Jmat[:,:] = 0
    for k in range(0,total):
        # 现在要求解k和k + 1处的变量
        colu = int(k / nmax)
        row = k % nmax
        idxk[0] = colu * nmax + row
        idxk[1] = idxk[0] + 1
        idxk[2] = idxk[1] + 1
        idxk[3] = idxk[0] + nmax
        idxk[4] = idxk[3] + 1
        idxk[5] = idxk[4] + 1
        idxk[6] = idxk[3] + nmax
        idxk[7] = idxk[6] + 1
        idxk[8] = idxk[7] + 1
        if(colu >= tmax - 2)|(row >= nmax - 2):
            continue
        for i in range(0,num):
            idx = int(idxk[i])
            td0 = d[0] # x^3
            td1 = 3*row*d[0] + d[1] # x^2
            td2 = 3*row*row*d[0] + 2*row*d[1] + d[2] + colu*colu*d[5] + colu*d[6] # x
            td3 = d[3] + row*d[5] # t^2
            td4 = d[3]*2*colu + d[4] + d[5]*2*row*colu + d[6]*row # t
            td5 = d[5] # xt^2
            td6 = d[5]*2*colu + d[6] # xt
            td7 = d[0]*row*row*row + d[1]*row*row + d[2]*row + (d[3]*colu*colu
                + d[4]*colu + d[5]*row*colu*colu + d[6]*row*colu + d[7])
            Fmat[idx] = Fmat[idx]  - (td0*rhs[0,i] + td1*rhs[1,i] + td2*rhs[2,i] + td3*rhs[3,i] + 
                                        td4*rhs[4,i] + td5*rhs[5,i] + td6*rhs[6,i] + td7*rhs[7,i])
        for p in range(0,num):
            for i in range(0,num):
                for j in range(0,num):
                    k0 = int(idxk[p])# 求解的位置
                    i0 = int(idxk[i])# u1的索引
                    j0 = int(idxk[j])# u2的索引
                    Fmat2[k0] = Fmat2[k0] + Fcoeff[p,i,j]*res[ite,i0]*res[ite,j0]
        for i in range(0,num):
            for j in range(0,num):
                i0 = int(idxk[i])
                j0 = int(idxk[j])
                Fmat3[i0] = Fmat3[i0] + Tmat[i,j]*res[ite,j0]
                Jmat[i0,j0] += Tmat[i,j]
        for i in range(0,num):
            for j in range(0,num):
                i0 = int(idxk[i])
                j0 = int(idxk[j])
                for p in range(0,num):
                    for q in range(0,num):
                        if p == j:
                            Jmat[i0,j0] += Fcoeff[i,p,q]*res[ite,int(idxk[q])]
                        if q == j:
                            Jmat[i0,j0] += Fcoeff[i,p,q]*res[ite,int(idxk[p])]
    for i in range(0,total):
        summ = 0
        for j in range(0,total):
            summ += Jmat[i,j]*res[ite,j]
        test[i] = summ
        
    # J du = F
    # 这里du_0 应是零，也就是u_0一直作为设定好的初始条件一直不变
    # 当算法最终收敛的时候，Fmat + Fmat2为零，这也是调试时重要的条件
    # Jmat[0,:] = 0
    # Jmat[0,0] = 1
    # Fmat[0] = Fmat2[0] = 0
    
    
    Jinv = np.linalg.inv(Jmat)
    for i in range(0,total):
        summ = 0
        for j in range(0,total):
            summ += Jinv[i,j]*(- Fmat[j] - Fmat2[j] - Fmat3[j])
        res[ite+1,i] = res[ite,i] + summ   