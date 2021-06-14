import numpy as np
# D:\FluidSim\FluidSim\Immersed\IB-MATLAB-master\IB_matlab_2D
L = 1 # 求解区域长度
N = 64 # 网格数量
h = L / N # 单个网格长度
Nb = int(np.ceil(np.pi*(L/2)/(h/2))) # 固体边界上粒子的数量
dtheta = 2 * np.pi / Nb
X = np.zeros((Nb,2)) # 粒子坐标
U = np.zeros((Nb,2)) # 粒子速度
F = np.zeros((Nb,2)) # 粒子力
for i in range(0,Nb):
    X[i,0] = L/2 + L/4*np.cos(dtheta*i)
    X[i,1] = L/2 + L/4*np.sin(dtheta*i)
    
u = np.zeros((N,N)) # 流体速度x轴
v = np.zeros((N,N)) # 流体速度y轴
for i in range(N):
    for j in range(N):
        v[i,j] = np.sin(2*np.pi*i*h/L)
    
pos = np.zeros((Nb,Nb))
w = np.zeros((4,Nb,4))
w1 = np.zeros((4,Nb,4))
w2 = np.zeros((4,Nb,4))
w3 = np.zeros((Nb,4,4))
dt = 0.01
elastic = 1
for t in range(1):
    
    # 第一大步，更新粒子的位置
    pos = X / h
    posi = np.asarray(pos,dtype = int)
    r = pos - posi
    
    q = np.sqrt(1 + 4*r*(1 - r))
    w1[3,:,0] = w1[3,:,1] = w1[3,:,2] = w1[3,:,3] = (1 + 2*r[:,1] - q[:,1])/8
    w1[2,:,0] = w1[2,:,1] = w1[2,:,2] = w1[2,:,3] = (1 + 2*r[:,1] + q[:,1])/8
    w1[1,:,0] = w1[1,:,1] = w1[1,:,2] = w1[1,:,3] = (3 - 2*r[:,1] + q[:,1])/8
    w1[0,:,0] = w1[0,:,1] = w1[0,:,2] = w1[0,:,3] = (3 - 2*r[:,1] - q[:,1])/8
    
    w2[0,:,3] = w2[1,:,3] = w2[2,:,3] = w2[3,:,3] = (1 + 2*r[:,0] - q[:,0])/8
    w2[0,:,2] = w2[1,:,2] = w2[2,:,2] = w2[3,:,2] = (1 + 2*r[:,0] + q[:,0])/8
    w2[0,:,1] = w2[1,:,1] = w2[2,:,1] = w2[3,:,1] = (3 - 2*r[:,0] + q[:,0])/8
    w2[0,:,0] = w2[1,:,0] = w2[2,:,0] = w2[3,:,0] = (3 - 2*r[:,0] - q[:,0])/8
    w = w1*w2
    for k in range(0,Nb):
        ww = w[:,k,:]
        U[k,0] = U[k,1] = 0
        for i in range(0,4):
            for j in range(0,4):
                x0 = int(posi[k,0] + i - 1)
                y0 = int(posi[k,1] + j - 1)
                U[k,0] += ww[i,j]*u[x0,y0]
                U[k,1] += ww[i,j]*v[x0,y0]
    
    newX = X + dt/2*U
    
    # 第二大步，根据粒子的位置，计算粒子弹力
    for i in range(0,Nb):
        x0 = np.mod(i - 1 + Nb,Nb)
        x1 = np.mod(i + 1,Nb)
        F[i,:] = elastic*(newX[x0,:] + newX[x1,:] 
                         - 2*newX[i,:])/(dtheta*dtheta)
    
    # 第三大步，根据粒子的弹力，计算网格上的弹力
    c = dtheta / h / h
    fx = np.zeros((N,N))
    fy = np.zeros((N,N))
    
    pos = newX / h
    posi = np.asarray(pos,dtype = int)
    r = pos - posi
    
    q = np.sqrt(1 + 4*r*(1 - r))
    w1[3,:,0] = w1[3,:,1] = w1[3,:,2] = w1[3,:,3] = (1 + 2*r[:,1] - q[:,1])/8
    w1[2,:,0] = w1[2,:,1] = w1[2,:,2] = w1[2,:,3] = (1 + 2*r[:,1] + q[:,1])/8
    w1[1,:,0] = w1[1,:,1] = w1[1,:,2] = w1[1,:,3] = (3 - 2*r[:,1] + q[:,1])/8
    w1[0,:,0] = w1[0,:,1] = w1[0,:,2] = w1[0,:,3] = (3 - 2*r[:,1] - q[:,1])/8
    
    w2[0,:,3] = w2[1,:,3] = w2[2,:,3] = w2[3,:,3] = (1 + 2*r[:,0] - q[:,0])/8
    w2[0,:,2] = w2[1,:,2] = w2[2,:,2] = w2[3,:,2] = (1 + 2*r[:,0] + q[:,0])/8
    w2[0,:,1] = w2[1,:,1] = w2[2,:,1] = w2[3,:,1] = (3 - 2*r[:,0] + q[:,0])/8
    w2[0,:,0] = w2[1,:,0] = w2[2,:,0] = w2[3,:,0] = (3 - 2*r[:,0] - q[:,0])/8
    
    for k in range(0,Nb):
        ww = w[:,k,:]
        for i in range(0,4):
            for j in range(0,4):
                x0 = int(posi[k,0] + i - 1)
                y0 = int(posi[k,1] + j - 1)
                fx[x0,y0] += ww[i,j]*F[k,0]*c
                fy[x0,y0] += ww[i,j]*F[k,1]*c
    
    # 第四大步，根据网格上的力，计算流体速度
    skewu = np.zeros((N,N))
    skewv = np.zeros((N,N))
    wu = np.zeros((N,N))
    wv = np.zeros((N,N))
    rho = 1
    for i in range(0,N):
        for j in range(0,N):
            x0 = np.mod(i - 1 + N,N)
            x1 = np.mod(i + 1,N)
            y0 = np.mod(j - 1 + N,N)
            y1 = np.mod(j + 1,N)
            skewu[i,j] = ((u[x1,j] + u[i,j])*u[x1,j] 
              - (u[x0,j] + u[i,j])*u[x0,j]
              + (v[i,y1] + v[i,j])*u[i,y1] 
              - (v[i,y0] + v[i,j])*u[i,y0])/(4*h)
            skewv[i,j] = ((u[x1,j] + u[i,j])*v[x1,j] 
              - (u[x0,j] + u[i,j])*v[x0,j]
              + (v[i,y1] + v[i,j])*v[i,y1] 
              - (v[i,y0] + v[i,j])*v[i,y0])/(4*h)
    wu = u - dt/2*skewu + dt/2/rho*fx
    wv = v - dt/2*skewv + dt/2/rho*fy
    
    # 第五大步，根据网格上的速度，计算粒子的速度
    
    