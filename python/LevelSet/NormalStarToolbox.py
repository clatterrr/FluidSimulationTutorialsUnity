import numpy as np
import math
'''
参考自MatlabLevelSetToolbox
normalStarDemo
'''
speed = 0.25
useTimeDependent = 0
tMax = 1 # 最终时间
plotSteps = 9 # 绘图时间
t0 = 0 # 开始时间
singleStep = 0 

tPlot = (tMax - t0) / (plotSteps - 1)
small = 1e-10

level = 0

periodic = 0 # 是否使用周期边界条件
nmax = 101
dx = 0.02
gxs0 = np.zeros((nmax,nmax))
gxs1 = np.zeros((nmax,nmax))
theta = np.zeros((nmax,nmax))
rad = np.zeros((nmax,nmax))
data0 = np.zeros((nmax,nmax))
scale = 0.2
points = 7
shift = 2.5
for i in range(0,nmax):
    for j in range(0,nmax):
        gxs0[i,j] = dx * i - 1
        gxs1[i,j] = dx * j - 1
        theta[i,j] = math.atan2(gxs1[i,j],gxs0[i,j])
        rad[i,j] = np.sqrt(gxs1[i,j]**2 + gxs0[i,j]**2)
        data0[i,j] = rad[i,j] - scale * (np.cos(points * theta[i,j]) + shift)

'''  upwind First ENO2 Start'''
data1 = np.zeros((nmax+4,nmax))
data1[2:nmax+2,:] = data0[:,:]
data1[1,:] = data1[2,:] + data1[2,:] - data1[3,:]
data1[0,:] = data1[1,:] + data1[1,:] - data1[2,:]

data1[nmax+2,:] = 2 * data1[nmax+1,:] - data1[nmax,:]
data1[nmax+3,:] = 2 * data1[nmax+2,:] - data1[nmax+1,:]

D1 = np.zeros((nmax + 3,nmax))
D2 = np.zeros((nmax + 2,nmax))
D1 = (data1[1:nmax+4,:] - data1[0:nmax+3,:]) / dx
D2 = 0.5 / dx * (D1[1:nmax+3,:] - D1[0:nmax+2,:])

dL1 = D1[1:nmax+1,:]
dL2 = D1[1:nmax+1,:]
dR1 = D1[2:nmax+2,:]
dR2 = D1[2:nmax+2,:]
dL1 = dL1 + dx * D2[0:nmax,:]
dL2 = dL2 + dx * D2[1:nmax+1,:]
dR1 = dR1 - dx * D2[1:nmax+1,:]
dR2 = dR2 - dx * D2[2:nmax+2,:]

derivL = np.zeros((nmax,nmax))
derivR = np.zeros((nmax,nmax))
for i in range(0,nmax):
    for j in range(0,nmax):
        smallerL = (abs(D2[i,j]) < abs(D2[i+1,j]))
        smallerR = 1 - smallerL
        derivL[i,j] = dL1[i,j] * smallerL + dL2[i,j] * smallerR
        
        smallerL = (abs(D2[i+1,j]) < abs(D2[i+2,j]))
        smallerR = 1 - smallerL
        derivR[i,j] = dR1[i,j] * smallerL + dR2[i,j] * smallerR
'''  upwind First ENO2 END
     termNormal Start
'''
flowL = np.zeros((nmax,nmax))
flowR = np.zeros((nmax,nmax))
magnitude = np.zeros((nmax,nmax))
stepBoundinv = np.zeros((nmax,nmax))
for i in range(0,nmax):
    for j in range(0,nmax):
        prodL = speed * derivL[i,j]
        prodR = speed * derivR[i,j]
        magL = abs(prodL)
        magR = abs(prodR)
        flowL[i,j] = ((prodL >= 0 )&(prodR >=0)) |  ((prodL >= 0) & (prodR <= 0) & (magL >= magR))
        flowR[i,j] = ((prodL <= 0 )&(prodR <=0)) |  ((prodL >= 0) & (prodR <= 0) & (magL < magR))
        magnitude[i,j] = magnitude[i,j] + derivL[i,j]**2*flowL[i,j] + derivR[i,j]**2*flowR[i,j]
        stepBoundinv[i,j] = (magL * flowL[i,j] + magR * flowR[i,j])/dx
'''
    termNormal End
    UpwindFirst Another Direction i.e. y-axis
'''
data1 = np.zeros((nmax,nmax+4))
data1[:,2:nmax+2] = data0[:,:]
data1[:,1] = data1[:,2] + data1[:,2] - data1[:,3]
data1[:,0] = data1[:,1] + data1[:,1] - data1[:,2]

data1[:,nmax+2] = 2 * data1[:,nmax+1] - data1[:,nmax]
data1[:,nmax+3] = 2 * data1[:,nmax+2] - data1[:,nmax+1]

D1 = np.zeros((nmax,nmax+3))
D2 = np.zeros((nmax,nmax+2))
D1 = (data1[:,1:nmax+4] - data1[:,0:nmax+3]) / dx
D2 = 0.5 / dx * (D1[:,1:nmax+3] - D1[:,0:nmax+2])

dL1 = D1[:,1:nmax+1]
dL2 = D1[:,1:nmax+1]
dR1 = D1[:,2:nmax+2]
dR2 = D1[:,2:nmax+2]
dL1 = dL1 + dx * D2[:,0:nmax]
dL2 = dL2 + dx * D2[:,1:nmax+1]
dR1 = dR1 - dx * D2[:,1:nmax+1]
dR2 = dR2 - dx * D2[:,2:nmax+2]

derivL = np.zeros((nmax,nmax))
derivR = np.zeros((nmax,nmax))
for i in range(0,nmax):
    for j in range(0,nmax):
        smallerL = (abs(D2[i,j]) < abs(D2[i,j+1]))
        smallerR = 1 - smallerL
        derivL[i,j] = dL1[i,j] * smallerL + dL2[i,j] * smallerR
        
        smallerL = (abs(D2[i,j+1]) < abs(D2[i,j+2]))
        smallerR = 1 - smallerL
        derivR[i,j] = dR1[i,j] * smallerL + dR2[i,j] * smallerR
'''  upwind First ENO2 END
     termNormal Start
'''
flowL = np.zeros((nmax,nmax))
flowR = np.zeros((nmax,nmax))
ydot = np.zeros((nmax,nmax))
stepBound = 0
for i in range(0,nmax):
    for j in range(0,nmax):
        prodL = speed * derivL[i,j]
        prodR = speed * derivR[i,j]
        magL = abs(prodL)
        magR = abs(prodR)
        flowL[i,j] = ((prodL >= 0 )&(prodR >=0)) |  ((prodL >= 0) & (prodR <= 0) & (magL >= magR))
        flowR[i,j] = ((prodL <= 0 )&(prodR <=0)) |  ((prodL >= 0) & (prodR <= 0) & (magL < magR))
        magnitude[i,j] = np.sqrt(magnitude[i,j] + derivL[i,j]**2*flowL[i,j] + derivR[i,j]**2*flowR[i,j])
        ydot[i,j] = - speed * magnitude[i,j]
        stepBoundinv[i,j] = stepBoundinv[i,j] + (magL * flowL[i,j] + magR * flowR[i,j])/dx
        if magnitude[i,j] != 0:
            stepBound = max(stepBound,stepBoundinv[i,j] / magnitude[i,j])
stepBound = 1 / stepBound
cfl = 0.5
TimeDelta = min(cfl * stepBound,1)