import numpy as np
'''
参考资料：https://www.pplusplus.lima-city.de/femfluid.html Pressure Solve with Finite Elements  很好的matlab库

代码地址：https://www.pplusplus.lima-city.de/lib/data/femfluid/FEM%20Fluid%20Source.zip

本地地址：D:\FluidSim\FluidSim\FEMNEW\FEM Fluid Source\FEM Fluid

我所做的就是把matlab代码翻译成python代码，过程和结果一模一样，翻译过程中学到了很多知识
'''
xcells = 99
ycells = 99
Lx = 1
Ly = 1
dx = Lx / xcells
dy = Ly / ycells
nmax = xcells + 1
mmax = ycells + 1
dt = 0.0083
disspation = 1
pressure = np.zeros((nmax,mmax)) # 压力

obstacle = np.zeros((nmax,mmax)) # 障碍物
obstacle[63:84,39:60] = 1
obstacle[nmax-2:nmax,:] = 1
obstacle[:,0:2] = 1
obstacle[:,nmax-2:nmax] = 1

import scipy.io as scio
fluidnodes = scio.loadmat('fluidnodes.mat')['fluidNodes']
distanceField = np.zeros((nmax,mmax))
distancesToAir = scio.loadmat('distancesToAir.mat')['distancesToAir']
distancesToWater = scio.loadmat('distancesToWater.mat')['distancesToWater']
for i in range(nmax):
    for j in range(mmax):
        if fluidnodes[i,j] == 0:
            distanceField[i,j] = distancesToWater[i,j]
        else:
            distanceField[i,j] = distancesToAir[i,j]

velocityFieldU = np.zeros((nmax,mmax))
velocityFieldV = np.ones((nmax,mmax))*0.58175
ForceFieldX = np.ones((nmax,mmax))*9.81
ForceFieldY = np.ones((nmax,mmax))*9.81
for i in range(nmax):
    for j in range(mmax):
        if obstacle[i,j] == 1:
            velocityFieldU[i,j] = 0
            velocityFieldV[i,j] = 0

# 不参与实际运算
frame = 1200
divergence = np.zeros((frame))
fluidPixels = np.zeros((frame))

for t in range(0,frame):
    # 第一大步，施加速度的影响
    newVelocityFieldU = velocityFieldU + ForceFieldX * dt
    newVelocityFieldV = velocityFieldV + ForceFieldY * dt
    
    # 障碍物之中不应该有速度
    for i in range(nmax):
        for j in range(mmax):
            if obstacle[i,j] == 1:
                velocityFieldU[i,j] = 0
                velocityFieldV[i,j] = 0
    
    # 第二大步，解算压力，获得无散速度场
    s = np.array([[-np.sqrt(1/3),-np.sqrt(1/3)],
                  [-np.sqrt(1/3),+np.sqrt(1/3)],
                  [+np.sqrt(1/3),-np.sqrt(1/3)],
                  [+np.sqrt(1/3),+np.sqrt(1/3)]])
    zeta = s[:,0]
    eta = s[:,1]
    
    # 第2.1步，混合流体部分
    fluidCells = np.zeros((xcells,ycells))
    for i in range(xcells):
        for j in range(ycells):
            flag0 = distanceField[i,j] >= 0
            flag1 = distanceField[i+1,j] >= 0
            flag2 = distanceField[i,j+1] >= 0
            flag3 = distanceField[i+1,j+1] >= 0
            fluidCells[i,j] = flag0 | flag1 | flag2 | flag3
            
    fluid = np.zeros((nmax,mmax),dtype = bool)
    for i in range(nmax):
        for j in range(mmax):
            fluid[i,j] = distanceField[i,j] >= 0
            
    # 第2.2 步，计算边界点
    # a node is interface node if:
    # it is not in an obstacle
    # at least one neighbor's fluid state is different than the middle
    # that neighbor node is not in an obstacle
    paddedfluid = np.zeros((nmax+2,mmax+2),dtype = bool)
    paddedobstacle = np.zeros((nmax+2,mmax+2),dtype = bool)
    paddedfluid[1:nmax+1,1:mmax+1] = fluid[:,:]
    paddedobstacle[1:nmax+1,1:mmax+1] = obstacle[:,:]
    interfaceNodes = np.zeros((nmax,mmax),dtype = bool)
    
    # 其实不用新建这么多变量，但为了方便调试和理解代码，就这样吧
    fluidLeft = np.zeros((nmax,mmax),dtype = bool)
    fluidRight = np.zeros((nmax,mmax),dtype = bool)
    fluidUp = np.zeros((nmax,mmax),dtype = bool)
    fluidDown = np.zeros((nmax,mmax),dtype = bool)
    fluidDownLeft = np.zeros((nmax,mmax),dtype = bool)
    fluidDownRight = np.zeros((nmax,mmax),dtype = bool)
    fluidUpLeft = np.zeros((nmax,mmax),dtype = bool)
    fluidUpRight = np.zeros((nmax,mmax),dtype = bool)
            
    obstacleLeft = np.zeros((nmax,mmax),dtype = bool)
    obstacleRight = np.zeros((nmax,mmax),dtype = bool)
    obstacleUp = np.zeros((nmax,mmax),dtype = bool)
    obstacleDown = np.zeros((nmax,mmax),dtype = bool)
    obstacleDownLeft = np.zeros((nmax,mmax),dtype = bool)
    obstacleDownRight = np.zeros((nmax,mmax),dtype = bool)
    obstacleUpLeft = np.zeros((nmax,mmax),dtype = bool)
    obstacleUpRight = np.zeros((nmax,mmax),dtype = bool)
    
    for i in range(1,nmax+1):
        for j in range(1,mmax+1):
            fluidLeft[i-1,j-1] = paddedfluid[i,j + 1]
            fluidRight[i-1,j-1] = paddedfluid[i,j - 1]
            fluidUp[i-1,j-1] = paddedfluid[i + 1,j]
            fluidDown[i-1,j-1] = paddedfluid[i - 1,j]
            fluidDownLeft[i-1,j-1] = paddedfluid[i - 1,j + 1]
            fluidDownRight[i-1,j-1] = paddedfluid[i - 1,j - 1]
            fluidUpLeft[i-1,j-1] = paddedfluid[i + 1,j + 1]
            fluidUpRight[i-1,j-1] = paddedfluid[i + 1,j - 1]
            
            obstacleLeft[i-1,j-1] = paddedobstacle[i,j + 1]
            obstacleRight[i-1,j-1] = paddedobstacle[i,j - 1]
            obstacleUp[i-1,j-1] = paddedobstacle[i + 1,j]
            obstacleDown[i-1,j-1] = paddedobstacle[i - 1,j]
            obstacleDownLeft[i-1,j-1] = paddedobstacle[i - 1,j + 1]
            obstacleDownRight[i-1,j-1] = paddedobstacle[i - 1,j - 1]
            obstacleUpLeft[i-1,j-1] = paddedobstacle[i + 1,j + 1]
            obstacleUpRight[i-1,j-1] = paddedobstacle[i + 1,j - 1]
            
    for i in range(0,nmax):
        for j in range(0,mmax):
            if obstacle[i,j] == 1:
                continue
            
            flag0 = fluid[i,j] != (fluidUp[i,j] & ~obstacleUp[i,j])
            flag1 = (fluid[i,j] != (fluidDown[i,j] & ~obstacleDown[i,j]))
            flag2 = (fluid[i,j] != (fluidRight[i,j] & ~obstacleRight[i,j]))
            flag3 = (fluid[i,j] != (fluidLeft[i,j] & ~obstacleLeft[i,j]))
            interfaceNodes[i,j] = interfaceNodes[i,j] | flag0 | flag1 | flag2 | flag3
            
            flag0 = fluid[i,j] != (fluidDownLeft[i,j] & ~obstacleDownLeft[i,j])
            flag1 = (fluid[i,j] != (fluidDownRight[i,j] & ~obstacleDownRight[i,j]))
            flag2 = (fluid[i,j] != (fluidUpLeft[i,j] & ~obstacleUpLeft[i,j]))
            flag3 = (fluid[i,j] != (fluidUpRight[i,j] & ~obstacleUpRight[i,j]))
            interfaceNodes[i,j] = interfaceNodes[i,j] | flag0 | flag1 | flag2 | flag3
            
    # 第2.3步，组装矩阵
    nodeIndex = 0
    nodeCoordToIndex = np.ones((nmax,mmax)) *(-1)
    nodeIndexToCoord = np.zeros((nmax*mmax,2))
    nodeOffsets = np.array([[0,0],[0,1],[1,1],[1,0]])
    total = 0
    for j in range(ycells):
        for i in range(xcells):
            if fluidCells[j,i] == 0:
                continue
            total += 1
    Kx = np.zeros((total * 16))
    Ky = np.zeros((total * 16))
    Kval = np.zeros((total * 16))
    fx = np.zeros((total * 4))
    fval = np.zeros((total * 4))
    total = -1
    for j in range(ycells):
        for i in range(xcells):
            if fluidCells[j,i] == 0:
                continue
            total += 1
            wTildeEX = np.zeros((4))
            wTildeEY = np.zeros((4))
            for corner in range(4):
                corner1 = i + nodeOffsets[corner,0]
                corner2 = j + nodeOffsets[corner,1]
                if nodeCoordToIndex[corner2,corner1] == -1:
                    nodeCoordToIndex[corner2,corner1] = nodeIndex
                    nodeIndexToCoord[nodeIndex,0] = corner1 
                    nodeIndexToCoord[nodeIndex,1] = corner2
                    nodeIndex += 1
                    
            Ke = np.zeros((4,4))
            Fe = np.zeros((4))
            for corner in range(4):
                zeta0 = zeta[corner]
                eta0 = eta[corner]
                # 求导后的矩阵
                currG = np.array([[eta0-1,zeta0-1],[-eta0-1,-zeta0+1],
                                  [eta0+1,zeta0+1],[-eta0+1,-zeta0-1]])/4
                Xtilde = np.array([[i,i,i+1,i+1],[j,j+1,j+1,j]])
                J = np.dot(Xtilde,currG)
                Jinv = np.linalg.inv(J)
                Bsi = np.dot(currG,Jinv) # 求导后的矩阵
                JsiDet = np.linalg.det(J)
                ShapeFunction = np.zeros((4))
                ShapeFunction[0] = (1 - zeta0)*(1 - eta0)/4
                ShapeFunction[1] = (1 - zeta0)*(1 + eta0)/4
                ShapeFunction[2] = (1 + zeta0)*(1 + eta0)/4
                ShapeFunction[3] = (1 + zeta0)*(1 - eta0)/4
                
                wi = 1
                Ke = Ke + wi * JsiDet * (np.dot(Bsi,np.transpose(Bsi)))
                Fe = Fe + wi * JsiDet * np.dot(
                    np.dot(Bsi[:,0],wTildeEX) + np.dot(Bsi[:,1],wTildeEY),ShapeFunction)          
                test =1
            newKx = np.zeros((16))
            newKy = np.zeros((16))
            newKval = np.zeros((16))
            newFx = np.zeros((4))
            vectorIndex = 0
            for corner in range(4):
                fromCornerCoord1 = i + nodeOffsets[corner,0]
                fromCornerCoord2 = j + nodeOffsets[corner,1]
                fromCornerIndex = int(nodeCoordToIndex[fromCornerCoord2,fromCornerCoord1])
                for ToCorner in range(4):
                    toCornerCoord1 = i + nodeOffsets[ToCorner,0]
                    toCornerCoord2 = j + nodeOffsets[ToCorner,1]
                    toCornerIndex = int(nodeCoordToIndex[toCornerCoord2,toCornerCoord1])
                    newKx[vectorIndex] = fromCornerIndex
                    newKy[vectorIndex] = toCornerIndex
                    newKval[vectorIndex] = Ke[corner,ToCorner]
                    vectorIndex += 1
                newFx[corner] = fromCornerIndex
            
            stidx = total*16
            edidx = total*16 + 16
            Kx[stidx:edidx] = newKx[:]
            Ky[stidx:edidx] = newKy[:]
            Kval[stidx:edidx] = newKval[:]
            stidx = total*4
            edidx = total*4 + 4
            fx[stidx:edidx] = newFx[:]
            fval[stidx:edidx] = Fe[:]

    # 第2.4步，移除边界节点
    maxidx = len(Kx)
    sparseInBoundary = np.zeros((maxidx),dtype = bool)
    fromNodeInBoundary = np.zeros((maxidx),dtype = bool)
    toNodeInBoundary = np.zeros((maxidx),dtype = bool)
    for i in range(0,maxidx):
        fromidx = int(Kx[i])
        x0 = int(nodeIndexToCoord[fromidx,0])
        y0 = int(nodeIndexToCoord[fromidx,1])
        if interfaceNodes[y0,x0] == 1:
            fromNodeInBoundary[i] = 1
            
        toidx = int(Ky[i])
        x0 = int(nodeIndexToCoord[toidx,0])
        y0 = int(nodeIndexToCoord[toidx,1])
        if interfaceNodes[y0,x0] == 1:
            toNodeInBoundary[i] = 1
        
        sparseInBoundary[i] = fromNodeInBoundary[i] | toNodeInBoundary[i]
        if sparseInBoundary[i] == 1:
            Kx[i] = -1
            Ky[i] = -1
            Kval[i] = 0
    maxidx = len(fx)
    sparseInBoundary = np.zeros((maxidx))
    for i in range(0,maxidx):
        idx = int(fx[i])
        x0 = int(nodeIndexToCoord[idx,0])
        y0 = int(nodeIndexToCoord[idx,1])
        if interfaceNodes[y0,x0] == 1:
            sparseInBoundary[i] = 1
            
        if sparseInBoundary[i] == 1:
            fx[i] = 0
            fval[i] = 0
    test = 1