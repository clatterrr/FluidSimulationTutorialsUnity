import numpy as np
"""
参考：有限元法与matlab程序设计 郭吉坦 薛齐文 编著

本地代码 D:\FluidSim\FluidSim\FEMNEW\有限元与MTALAB程序设计程序文件\有限元与MTALAB程序设计程序文件\matlabfile

第七章 杆系结构
"""
E = 2e8 # 弹性模量 kN / m^2
A = 0.006 # 截面面积 m^2
nd = 4 # 节点总数
pos = np.array([[0,0],[2,0],[0,2],[2,2]],dtype = int)
ne = 6 # 单元总数
ele = np.array([[0,1,0],[0,2,0],[0,3,0],[2,1,0],[2,3,0],[1,3,0]])
ng = 3 # 约束总数
# 节点0的xy的方向都被约束，节点1的y方向被约束
boundary = np.array([[0,1],[0,0],[1,0]])
nj = 1 # 集中力个数
QJ = np.array([4,1,2])
Kmat = np.zeros((nd*2,nd*2))
for i in range(0,ne):
    dx = pos[ele[i,1],0] - pos[ele[i,0],0]
    dy = pos[ele[i,1],1] - pos[ele[i,0],1]
    dis = np.sqrt(dx*dx + dy*dy)
    c = dx / dis
    s = dy / dis
    Ke = np.array([[c*c,c*s,-c*c,-c*s],[c*s,s*s,-c*s,-s*s],
                   [-c*c,-c*s,c*c,c*s],[-c*s,s*s,c*s,s*s]])*E*A/dis
    idx = np.zeros((4),dtype = int)
    idx[0] = ele[i,0]*2 # 节点0的x轴
    idx[1] = ele[i,0]*2 + 1 # 节点0的y轴
    idx[2] = ele[i,1]*2 # 节点1的x轴
    idx[3] = ele[i,1]*2 + 1 # 节点1的y轴
    for j in range(0,4):
        for k in range(0,4):
            Kmat[idx[j],idx[k]] += Ke[j,k]
            