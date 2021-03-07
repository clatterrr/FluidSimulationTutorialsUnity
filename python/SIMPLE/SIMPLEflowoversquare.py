import numpy as np

# 使用SIMPLE算法解决方块绕流问题，交错网格
# 教程 https://zhuanlan.zhihu.com/p/347410166
# 并没有卡门涡街现象出现，原因未知


nmax = 16  # 有多长
mmax = 6  # 有多高
tmax = 20
u = np.zeros((nmax + 1, mmax + 2))  # 每回合初始猜的，和最终值
u_star = np.zeros((nmax + 1, mmax + 2))  # 用于迭代

v = np.zeros((nmax + 2, mmax + 1))  # 每回合初始猜的，和最终值
v_star = np.zeros((nmax + 2, mmax + 1))  # 用于迭代

p = np.zeros((nmax + 2, mmax + 2))  # 每回合初始猜的，和最终值
p_star = np.zeros((nmax + 2, mmax + 2))  # 用于迭代
p_source = np.zeros((nmax + 2, mmax + 2))  # 其实这东西就是散度，用于解泊松压力方程

Ae = np.zeros((nmax + 2, mmax + 2))
Aw = np.zeros((nmax + 2, mmax + 2))
An = np.zeros((nmax + 2, mmax + 2))
As = np.zeros((nmax + 2, mmax + 2))
Apu = np.zeros((nmax + 2, mmax + 2))
Apv = np.zeros((nmax + 2, mmax + 2))

# 用于解压力的动量方程
Ape = np.zeros((nmax + 2, mmax + 2))
Apw = np.zeros((nmax + 2, mmax + 2))
Apn = np.zeros((nmax + 2, mmax + 2))
Aps = np.zeros((nmax + 2, mmax + 2))
App = np.zeros((nmax + 2, mmax + 2))

# 用于验算每个时间步后速度是否是无散的
div = np.zeros((nmax + 2, mmax + 2))
# 寻找涡量
omega = np.zeros((tmax,nmax + 2, mmax + 2))

Re = 200
dt = 0.01
dx = 1 / nmax
dy = 1 / mmax
rho = 1

stx = 4#左下角
sty = 3#左下角
edx = 6#右上角
edy = 5#右上角
u[0, 1:mmax+1] = 1
u[nmax,1:mmax+1] = 1
block = np.zeros((nmax+2,mmax+2))
for i in range(stx,edx):
    for j in range(sty,edy):
        block[i,j] = 1
block[:,0] = block[:,-1] = 1

alpha_uv = 0.3  # 速度动量方程中的松弛系数
alpha_p = 0.2  # 校正压力时的松弛系数

vt = np.zeros((tmax,nmax + 2, mmax + 1))


for t in range(0, tmax):
    u_star = u.copy()
    v_star = v.copy()
    p_star = p.copy()
    # 第二步，计算u动量方程
    Ae[:, :] = Aw[:, :] = An[:, :] = As[:, :] = 0
    for i in range(1, nmax):
        for j in range(1, mmax + 1):
            Fe = 0.5 * (u[i, j] + u[i + 1, j]) * dy * rho
            Fw = 0.5 * (u[i, j] + u[i - 1, j]) * dy * rho
            Fn = 0.5 * (v[i, j] + v[i + 1, j]) * dx * rho
            Fs = 0.5 * (v[i, j - 1] + v[i + 1, j - 1]) * dx * rho
            Ae[i, j] = max(0, -Fe) + dy / dx / Re
            Aw[i, j] = max(0, Fw) + dy / dx / Re
            An[i, j] = max(0, -Fn) + dx / dy / Re
            As[i, j] = max(0, Fs) + dx / dy / Re
            Apu[i, j] = Ae[i, j] + Aw[i, j] + An[i, j] + As[i, j] + Fe - Fw + Fn - Fs

    # 由于交错网格的设置，最后一格速度到墙壁只有半格的距离，也就是dx/2或dy/2
    for i in range(1, nmax):
        As[i, 1] += dy / dx / Re
        An[i, mmax] += dy / dx / Re
        Apu[i, 1] += dy / dx / Re
        Apu[i, mmax] += dy / dx / Re

    for k in range(0,100):
        for i in range(1, nmax):
            for j in range(1, mmax + 1):
                if(block[i,j] == 1):
                    continue
                if(block[i+1,j] == 1):
                    continue
                term = Ae[i, j] * u_star[i + 1, j] + Aw[i, j] * u_star[i - 1, j]
                term += As[i, j] * u_star[i, j - 1] + An[i, j] * u_star[i, j + 1]
                pressure_term = (p[i, j] - p[i + 1, j]) * dy
                u_star[i, j] = u_star[i, j] + alpha_uv * ((term + pressure_term) / Apu[i, j] - u_star[i, j])

    test = 1
    # 第二步，计算v动量方程
    Ae[:, :] = Aw[:, :] = An[:, :] = As[:, :] = 0
    for i in range(1, nmax + 1):
        for j in range(1, mmax):
            Fe = 0.5 * (u[i, j] + u[i, j + 1]) * dy * rho
            Fw = 0.5 * (u[i - 1, j] + u[i - 1, j + 1]) * dy * rho
            Fn = 0.5 * (v[i, j] + v[i, j + 1]) * dx * rho
            Fs = 0.5 * (v[i, j] + v[i, j - 1]) * dx * rho
            Ae[i, j] = max(0, -Fe) + dy / dx / Re
            Aw[i, j] = max(0, Fw) + dy / dx / Re
            An[i, j] = max(0, -Fn) + dx / dy / Re
            As[i, j] = max(0, Fs) + dx / dy / Re
            Apv[i, j] = Ae[i, j] + Aw[i, j] + An[i, j] + As[i, j] + Fe - Fw + Fn - Fs

    for j in range(1, mmax):
        Ae[nmax, j] += dx / dy / Re
        Aw[1, j] += dx / dy / Re
        Apv[nmax, j] += dx / dy / Re
        Apv[1, j] += dx / dy / Re

    for k in range(0, 100):
        for i in range(1, nmax + 1):
            for j in range(1, mmax):
                if(block[i,j] == 1):
                    continue
                if(block[i,j+1] == 1):
                    continue
                term = Ae[i, j] * v_star[i + 1, j] + Aw[i, j] * v_star[i - 1, j]
                term += As[i, j] * v_star[i, j - 1] + An[i, j] * v_star[i, j + 1]
                pressure_term = (p[i, j] - p[i, j + 1]) * dx
                v_star[i, j] = v_star[i, j] + alpha_uv * ((term + pressure_term) / Apv[i, j] - v_star[i, j])

    total = 0
    for i in range(1, nmax + 1):
        for j in range(1, mmax + 1):
            if i < nmax:
                Ape[i, j] = rho * dy * dy / Apu[i, j]
            if i > 1:
                Apw[i, j] = rho * dy * dy / Apu[i - 1, j]
            if j < mmax:
                Apn[i, j] = rho * dx * dx / Apv[i, j]
            if j > 1:
                Aps[i, j] = rho * dx * dx / Apv[i, j - 1]
            # Ape[i,j] = Apw[i,j] = Apn[i,j] = Aps[i,j] = 1
            App[i, j] = Ape[i, j] + Apw[i, j] + Apn[i, j] + Aps[i, j]
            p_source[i, j] = (u_star[i, j] - u_star[i - 1, j]) * dy + (v_star[i, j] - v_star[i, j - 1]) * dx
            total += p_source[i, j]
    # 迭代次数越多，最后算出的散度越低，越接近无压缩流体的性质
    # 其实这里更应该用变化值来决定是否继续下去，当变化幅度很小，就停止迭代
    for k in range(0, 12000):
        for i in range(1, nmax + 1):
            for j in range(1, mmax + 1):
                if(block[i,j] == 1):
                    continue
                
                if(block[i+1,j] == 1):
                    p_star[i+1,j] = p_star[i,j]
                if(block[i-1,j] == 1):
                    p_star[i-1,j] = p_star[i,j]
                term = Ape[i, j] * p_star[i + 1, j] + Apw[i, j] * p_star[i - 1, j]
                
                if(block[i,j+1] == 1):
                    p_star[i,j+1] = p_star[i,j]
                if(block[i,j-1] == 1):
                    p_star[i,j-1] = p_star[i,j]
                term += Apn[i, j] * p_star[i, j + 1] + Aps[i, j] * p_star[i, j - 1]
                term = (term - p_source[i, j]) / App[i, j] # 就是公式那个d
                p_star[i, j] = p_star[i, j] + 0.8 * (term - p_star[i, j])
    # 校正压力和速度
    totalp = 0
    for i in range(1, nmax + 1):
        for j in range(1, mmax + 1):
            if(block[i,j] == 1):
                p_star[i,j] = 0
                continue
            p[i, j] = p[i, j] + alpha_p * (p_star[i, j] - p[i,j])
            totalp += p[i,j]

    # 压力的边界条件，
    p[:, 0] = p[:, 1]
    p[:, mmax + 1] = p[:, mmax]
    p[0,:] = 0  # 自由出入口的压力值是固定值，要根据实际的物理情况来设置
    p[nmax + 1, :] = 0  # 不过这里我就随便设置了两个

    for i in range(1, nmax):
        for j in range(1, mmax + 1):
            # 必须保证Apw[i + 1, j] == Ape[i,j]
            if(block[i,j] == 1):
                continue
            if(block[i+1,j] == 1):
                continue
            u[i, j] = u_star[i, j] + Ape[i, j] * (p_star[i, j] - p_star[i + 1, j]) / dy

    for i in range(1, nmax + 1):
        for j in range(1, mmax):
            # 必须保证Aps[i, j+1] == Apn[i,j]
            if(block[i,j] == 1):
                continue
            if(block[i,j+1] == 1):
                continue
            v[i, j] = v_star[i, j] + Apn[i, j] * (p_star[i, j] - p_star[i, j + 1]) / dx

    totald = 0
    for i in range(1, nmax + 1):
        for j in range(1, mmax + 1):
            # 由于最后算出的速度要是无散的，所以每次都要保证散度为零
            div[i, j] = (v[i, j] - v[i, j - 1]) * dx + (u[i, j] - u[i - 1, j]) * dy
            omega[t,i,j] = (v[i, j] - v[i, j - 1]) * dy + (u[i, j] - u[i - 1, j]) * dx
            totald += div[i,j]
    
    print(t)
    vt[t,:,:] = v[:,:]