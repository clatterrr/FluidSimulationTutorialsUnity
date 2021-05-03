import numpy as np
# 有限元法解决一维结构力学问题, 粒子排布在一水平条线上，仅仅受到垂直方向上重力与弹力的作用
nmax = 5 # 多少个结点
tmax = 100
force = np.zeros((nmax)) # 受到的力
pos = np.zeros((tmax + 1,nmax)) # 位移
vel = np.zeros((tmax + 1,nmax)) # 位移
for i in range(0,nmax):
    pos[0,i] = i
Kmat = np.zeros((nmax,nmax)) # 刚度矩阵
gravity = 1
A = 1 # 面积
E = 1 # 杨氏模量
L = 1 # 长度
Damp = 0.2
fixed = np.zeros((nmax))
fixed[0] = 1 # 1 号节点固定
constrain = np.zeros((nmax,2)) # 自己与哪些节点连接，2意味着左右连接
for i in range(0,nmax):
    if i == 0:
        constrain[i,0] = -1
    else:
        constrain[i,0] = i - 1
    if i == nmax-1:
        constrain[i,1] = -1
    else:
        constrain[i,1] = i + 1
for t in range(0,tmax):
    Kmat[:,:] = 0
    force[:] = 0
    for i in range(0,nmax):
        if fixed[i] == 1:
            Kmat[i,i] = 1
            force[i] = 0
        else:
            idx = int(i)
            idxleft = int(constrain[i,0]) 
            idxright = int(constrain[i,1])
            if idxleft >= 0:
                Kmat[idx,idx] += A*E/L
                Kmat[idx,idxleft] += -A*E/L
                force[idx] += (pos[t,idxleft] - pos[t,idx])*E
            if idxright >= 0:
                Kmat[idx,idx] += A*E/L
                Kmat[idx,idxright] += -A*E/L
                force[idx] += (pos[t,idxright] - pos[t,idx])*E
            force[idx] -= 1
    Kinv = np.linalg.inv(Kmat)
    for i in range(0,nmax):
        summ = 0
        for j in range(0,nmax):
            summ += Kinv[i,j]*force[j]
        # 这个0.8就是阻尼矩阵，作用于速度上。如果为一，那么杆模型将一直上下振动
        vel[t+1,i] = (1-Damp)*vel[t,i] + summ 
        pos[t+1,i] = pos[t,i] + vel[t+1,i]
        
