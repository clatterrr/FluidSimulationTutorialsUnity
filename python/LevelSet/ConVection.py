import numpy as np
from matplotlib import pyplot as plt
tmax = 1
t0 = 0
nmax = 201
data = np.zeros((nmax,nmax))
gxs = np.zeros((nmax,nmax))
gys = np.zeros((nmax,nmax))
ydot = np.zeros((nmax,nmax))
radius = 0.35
dx = 0.01
for i in range(0,nmax):
    for j in range(0,nmax):
        gxs[i,j] = dx * i - 1
        gys[i,j] = dx * j - 1
        data[i,j] = np.sqrt((gxs[i,j] - ( -0.4))**2 + (gys[i,j] - 0)**2) - radius
        
# 三阶精度
def RungeKutta3():
    
    gdata = np.zeros((nmax + 6,nmax))
    gdata[3:nmax+3,:] = data[:,:]
    gradient = gdata[3,:] - gdata[4,:]
    gdata[2,:] = gdata[3,:] + gradient
    gdata[1,:] = gdata[2,:] + gradient
    gdata[0,:] = gdata[1,:] + gradient
    gradient = gdata[nmax+2,:] - gdata[nmax+1,:]
    gdata[nmax+3,:] = gdata[nmax+2,:] + gradient
    gdata[nmax+4,:] = gdata[nmax+3,:] + gradient
    gdata[nmax+5,:] = gdata[nmax+4,:] + gradient
    
    D1 = 1 / dx * (gdata[1:nmax+6,:] - gdata[0:nmax+5,:])
    D2 = 0.5 / dx *(D1[1:nmax+5,:] - D1[0:nmax+4,:])
    D3 = 1 / 3 / dx *(D2[1:nmax+4,:] - D2[0:nmax+3,:])
    
    dL1 = D1[2:nmax+2,:].copy()
    dL2 = D1[2:nmax+2,:].copy()
    dL3 = D1[2:nmax+2,:].copy()
    dR1 = D1[3:nmax+3,:].copy()
    dR2 = D1[3:nmax+3,:].copy()
    dR3 = D1[3:nmax+3,:].copy()
    
    coeffL = 1 *  dx
    coeffR = -1 * dx
    dL1 = dL1 + coeffL * D2[1:nmax+1,:]
    dL2 = dL2 + coeffL * D2[1:nmax+1,:]
    dL3 = dL3 + coeffL * D2[2:nmax+2,:]
    dR1 = dR1 + coeffR * D2[2:nmax+2,:]
    dR2 = dR2 + coeffR * D2[2:nmax+2,:]
    dR3 = dR3 + coeffR * D2[3:nmax+3,:]
    
    coeffLL = 2 * dx * dx
    coeffLR = - dx * dx
    coeffRL = - dx * dx
    coeffRR = 2 * dx * dx
    # 原代码这里似乎写错了
    dL1 = dL1 + coeffLL * D3[0:nmax,:]
    dL2 = dL2 + coeffLL * D3[1:nmax+1,:]
    dL3 = dL3 + coeffLR * D3[2:nmax+2,:]
    dR1 = dR1 + coeffRL * D3[1:nmax+1,:]
    dR2 = dR2 + coeffRR * D3[2:nmax+2,:]
    dR3 = dR3 + coeffRR * D3[3:nmax+3,:]
    
    derivL = np.zeros((nmax,nmax))
    derivR = np.zeros((nmax,nmax))
    smallerL = np.zeros((nmax + 1,nmax), dtype=bool)
    smallerR = np.zeros((nmax + 1,nmax), dtype=bool)
    
    for i in range(0,nmax+1):
        for j in range(0,nmax):
            smallerL[i,j] = (abs(D2[i+1,j]) < abs(D2[i+2,j]))
            smallerR[i,j] = 1 - smallerL[i,j]
    smallerTemp = np.zeros((nmax + 2,nmax), dtype=bool)
    for i in range(0,nmax + 2):
        for j in range(0,nmax):
            smallerTemp[i,j] = (abs(D3[i,j]) < abs(D3[i+1,j]))
    smallerLL = np.zeros((nmax+1,nmax), dtype=bool)
    smallerRL = np.zeros((nmax+1,nmax), dtype=bool)
    smallerLR = np.zeros((nmax+1,nmax), dtype=bool)
    smallerRR = np.zeros((nmax+1,nmax), dtype=bool)
    smallerM = np.zeros((nmax+1,nmax), dtype=bool)
    
    for i in range(0,nmax+1):
        for j in range(0,nmax):
            temp1 = smallerTemp[i,j]
            temp2 = smallerTemp[i+1,j]
            smallerLL[i,j] = temp1 & smallerL[i,j]
            smallerRL[i,j] = temp2 & smallerR[i,j]
            
            temp1 = ~temp1
            temp2 = ~temp2
            smallerLR[i,j] = temp1 & smallerL[i,j]
            smallerRR[i,j] = temp2 & smallerR[i,j]
            smallerM[i,j] = smallerLR[i,j] | smallerRL[i,j]
            
    for i in range(0,nmax):
        for j in range(0,nmax):
            derivL[i,j] = (dL1[i,j]*smallerLL[i,j] + 
                dL2[i,j]*smallerM[i,j] + dL3[i,j]*smallerRR[i,j])
            derivR[i,j] = (dR1[i,j]*smallerLL[i+1,j] + 
                dR2[i,j]*smallerM[i+1,j] + dR3[i,j]*smallerRR[i+1,j])
    
    delta = np.zeros((nmax,nmax))
    stepBoundInv = 0
    vx = 2
    flowL = vx < 0
    flowR = vx > 0
    deriv = derivL * flowL + derivR * flowR
    delta = delta + deriv * vx
    
    gdata = np.zeros((nmax,nmax+6))
    gdata[:,3:nmax+3] = data[:,:]
    gradient = gdata[:,3] - gdata[:,4]
    gdata[:,2] = gdata[:,3] + gradient
    gdata[:,1] = gdata[:,2] + gradient
    gdata[:,0] = gdata[:,1] + gradient
    gradient = gdata[:,nmax+2] - gdata[:,nmax+1]
    gdata[:,nmax+3] = gdata[:,nmax+2] + gradient
    gdata[:,nmax+4] = gdata[:,nmax+3] + gradient
    gdata[:,nmax+5] = gdata[:,nmax+4] + gradient
    
    D1 = 1 / dx * (gdata[:,1:nmax+6] - gdata[:,0:nmax+5])
    D2 = 0.5 / dx *(D1[:,1:nmax+5] - D1[:,0:nmax+4])
    D3 = 1 / 3 / dx *(D2[:,1:nmax+4] - D2[:,0:nmax+3])
    
    dL1 = D1[:,2:nmax+2].copy()
    dL2 = D1[:,2:nmax+2].copy()
    dL3 = D1[:,2:nmax+2].copy()
    dR1 = D1[:,3:nmax+3].copy()
    dR2 = D1[:,3:nmax+3].copy()
    dR3 = D1[:,3:nmax+3].copy()
    
    coeffL = 1 *  dx
    coeffR = -1 * dx
    dL1 = dL1 + coeffL * D2[:,1:nmax+1]
    dL2 = dL2 + coeffL * D2[:,1:nmax+1]
    dL3 = dL3 + coeffL * D2[:,2:nmax+2]
    dR1 = dR1 + coeffR * D2[:,2:nmax+2]
    dR2 = dR2 + coeffR * D2[:,2:nmax+2]
    dR3 = dR3 + coeffR * D2[:,3:nmax+3]
    
    coeffLL = 2 * dx * dx
    coeffLR = - dx * dx
    coeffRL = - dx * dx
    coeffRR = 2 * dx * dx
    # 原代码这里似乎写错了
    dL1 = dL1 + coeffLL * D3[:,0:nmax]
    dL2 = dL2 + coeffLL * D3[:,1:nmax+1]
    dL3 = dL3 + coeffLR * D3[:,2:nmax+2]
    dR1 = dR1 + coeffRL * D3[:,1:nmax+1]
    dR2 = dR2 + coeffRR * D3[:,2:nmax+2]
    dR3 = dR3 + coeffRR * D3[:,3:nmax+3]
    
    derivL = np.zeros((nmax,nmax))
    derivR = np.zeros((nmax,nmax))
    smallerL = np.zeros((nmax,nmax+1), dtype=bool)
    smallerR = np.zeros((nmax,nmax+1), dtype=bool)
    
    for i in range(0,nmax):
        for j in range(0,nmax+1):
            smallerL[i,j] = (abs(D2[i,j+1]) < abs(D2[i,j+2]))
            smallerR[i,j] = 1 - smallerL[i,j]
    smallerTemp = np.zeros((nmax,nmax+2), dtype=bool)
    for i in range(0,nmax):
        for j in range(0,nmax+2):
            smallerTemp[i,j] = (abs(D3[i,j]) < abs(D3[i,j+1]))
    smallerLL = np.zeros((nmax,nmax+1), dtype=bool)
    smallerRL = np.zeros((nmax,nmax+1), dtype=bool)
    smallerLR = np.zeros((nmax,nmax+1), dtype=bool)
    smallerRR = np.zeros((nmax,nmax+1), dtype=bool)
    smallerM = np.zeros((nmax,nmax+1), dtype=bool)
    
    for i in range(0,nmax):
        for j in range(0,nmax+1):
            temp1 = smallerTemp[i,j]
            temp2 = smallerTemp[i,j+1]
            smallerLL[i,j] = temp1 & smallerL[i,j]
            smallerRL[i,j] = temp2 & smallerR[i,j]
            
            temp1 = ~temp1
            temp2 = ~temp2
            smallerLR[i,j] = temp1 & smallerL[i,j]
            smallerRR[i,j] = temp2 & smallerR[i,j]
            smallerM[i,j] = smallerLR[i,j] | smallerRL[i,j]
            
    for i in range(0,nmax):
        for j in range(0,nmax):
            derivL[i,j] = (dL1[i,j]*smallerLL[i,j] + 
                dL2[i,j]*smallerM[i,j] + dL3[i,j]*smallerRR[i,j])
            derivR[i,j] = (dR1[i,j]*smallerLL[i,j+1] + 
                dR2[i,j]*smallerM[i,j+1] + dR3[i,j]*smallerRR[i,j+1])
    
    vy = 0
    flowL = vy < 0
    flowR = vy > 0
    deriv = derivL * flowL + derivR * flowR
    delta = delta + deriv * vx
    stepBoundInv = stepBoundInv + vx / dx
    
    return 1 / stepBoundInv,- delta

for t in range(0,100):
    cfl = 0.5
    Data0 = data.copy()
    
    TimeDelta,ydot = RungeKutta3()
    data = data + cfl * TimeDelta * ydot
    Data1 = data.copy()
    
    TimeDelta,ydot = RungeKutta3()
    data = data + cfl * TimeDelta * ydot
    Data2 = data.copy()
    
    TimeDelta,ydot = RungeKutta3()
    data = data + cfl * TimeDelta * ydot
    Data3 = data.copy()    
    
    data = Data0 / 3 + Data1 / 2 + Data3 / 6
    plt.imshow(data, interpolation='none')
    plt.show()


