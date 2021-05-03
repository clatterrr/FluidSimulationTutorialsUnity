import numpy as np
# 有限元法解决一维结构力学问题，粒子排布在一垂直线上，仅仅受到垂直方向上重力与弹力的作用
nmax = 10 # 多少个结点
tmax = 1000
force = np.zeros((nmax)) # 受到的力
Kmat = np.zeros((nmax,nmax)) # 刚度矩阵
g = 0.1
A = 1 # 截面面积
E = 3 # 杨氏模量
L = 1.5 # 长度
dt = 0.15
damp = 0
mass = 1
fixed = np.zeros((nmax))# 0为可以任意移动 1 为位置不变，2为仅受到重力影响
fixed[0] = 1
constrain = np.zeros((nmax,2)) # 自己与哪些节点连接，2意味着左右连接
pos = np.zeros((tmax + 1,nmax)) # 位移
vel = np.zeros((tmax + 1,nmax)) # 位移
pos[0,0] = pos[1,0] = 0
for i in range(1,nmax):
    pos[0,i] = pos[1,i] = pos[0,i-1] + L - (nmax - i)*g/E + 2
for i in range(0,nmax):
    if i == 0:
        constrain[i,0] = -1
    else:
        constrain[i,0] = i - 1
    if i == nmax-1:
        constrain[i,1] = -1
    else:
        constrain[i,1] = i + 1
for t in range(1,tmax):
    Kmat[:,:] = 0
    force[:] = 0
    for i in range(0,nmax):
        if fixed[i] == 1:
            Kmat[i,i] = 1
            force[i] = 0
        elif fixed[i] == 2:
            Kmat[i,i] = 1
            force[i] = -g
        else:
            # x 方向上
            idx = int(i)
            idxdown = int(constrain[i,0])
            idxup = int(constrain[i,1])
            if idxdown >= 0:
                Kmat[idx,idx] += A*E/L
                Kmat[idx,idxdown] += -A*E/L
                if idxdown == 0:
                    force[idx] += (L*0.02 - (pos[t,idx] - pos[t,idxdown]))*E/100
                else:
                    force[idx] += (L - (pos[t,idx] - pos[t,idxdown]))*E
            if idxup >= 0:
                Kmat[idx,idx] += A*E/L
                Kmat[idx,idxup] += -A*E/L
                force[idx] += ((pos[t,idxup] - pos[t,idx]) - L)*E
            force[idx] -= g
    Kinv = np.linalg.inv(Kmat)
    for i in range(0,nmax):
        summ = 0
        for j in range(0,nmax):
            summ += Kinv[i,j]*force[j]
        # 这个0.8就是阻尼矩阵，作用于速度上。如果为一，那么杆模型将一直上下振动
        # vel[t+1,i] = (1-damp)*(vel[t,i] + dt*summ) 
        # pos[t+1,i] = pos[t,i] + dt*vel[t+1,i]
        pos[t+1,i] = pos[t,i] + (1 - damp)*(pos[t,i] - pos[t-1,i])*dt + dt*dt*summ
        if pos[t+1,i] < 0:
            pos[t+1,i] = 0
        
