import numpy as np
#Lax-Friedriches For 1d Shock Tube
tmax = 2000
nmax = 20
#F:\Projects\FluidSim\FVM\ApproximateRiemannSolvers-master\MUSCL_TVD
#成功了，但一点也不理想


rho = np.zeros((tmax,nmax))#密度
u = np.zeros((tmax,nmax))#速度
E = np.zeros((tmax,nmax))#能量
p = np.zeros((tmax,nmax))#压力
e = np.zeros((tmax,nmax))
c = np.zeros((tmax,nmax))#声速
H = np.zeros((tmax,nmax))#焓变

U1 = np.zeros((tmax,nmax))#rho
U2 = np.zeros((tmax,nmax))#rho*u
U3 = np.zeros((tmax,nmax))#rho*E
rhos = np.zeros((tmax,nmax))
Utemp = np.zeros((3,nmax))
Res = np.zeros((3,nmax))

gamma = 1.4
dt = 0.042
dx = 1/(nmax-1)
cfl = 0.5

for i in range(0,nmax):
    if(i < nmax/2):
        U1[0,i] = 1
        U2[0,i] = 0.9
        U3[0,i] = 2.5
    else:
        U1[0,i] = 0.125
        U3[0,i] = 0.25


def minmod(a,b,c):
    s = (np.sign(a) + np.sign(b) + np.sign(c))/3
    mm = 0
    if(abs(s) == 1):
        mm = s*min(abs(a),min(abs(b),abs(c)))
    return mm

def LFflux(qL,qR):
    #就是dU/dt + dF/dx = 0中的F
    rL = qL[0]
    uL = qL[1]/rL
    pL = (gamma-1)*(qL[2] - rL*uL*uL/2)
    HL = (qL[2] + pL)/rL
    
    rR = qR[0]
    uR = qR[1]/rR
    pR = (gamma-1)*(qR[2] - rR*uR*uR/2)
    HR = (qR[2] + pR)/rR
    
    smax = 1
    fx = np.zeros((3))
    fx[0] = 0.5*(rL*uL + rR*uR + smax*(qL[0] - qR[0]))
    fx[1] = 0.5*(rL*uL*uL+pL + rR*uR*uR+pR + smax*(qL[1] - qR[1]))
    fx[2] = 0.5*(uL*(rL*HL)  + uR*(rR*HR) + smax*(qL[2] - qR[2]))
    
    return fx 

def HLLflux(qL,qR):
    rL = qL[0]
    uL = qL[1]/rL
    pL = (gamma-1)*(qL[2] - rL*uL*uL/2)
    aL = np.sqrt(gamma * pL/rL)
    
    rR = qR[0]
    uR = qR[1]/rR
    pR = (gamma-1)*(qR[2] - rR*uR*uR/2)
    aR = np.sqrt(gamma * pR/rR)
    
    #Compute Wave Speed
    s1 = min(uL - aL,uR - aR)
    s3 = max(uL + aL,uR + aR)
    s2 = (pR - pL + rL*uL*(s1 - uL) - rR*uR*(s3 - uR))/(rL*(s1-uL)-rR*(s3-uR))
    
    uL_vec = np.vstack(uL,uL,uL)
    nL = np.vstack(0,1,uL)
    uR_vec = np.vstack(uR,uR,uR)
    nR = np.vstack(0,1,uR)
    n_s = np.vstack(0,1,s2)
    p_s = rL*(uL - s1)*(uL - s2) + pL
    qL_s = ((uL_vec-s1)*qL + (pL*nL - p_s*n_s))/(s2-s1)
    qR_s = ((uR_vec-s3)*qR + (pR*nR - p_s*n_s))/(s2-s3)
    
    #Compute jumps
    W1 = qL_s - qL
    W2 = qR_s - qL_s
    W3 = qR - qR_s
    
    
def MUSCL():
    
    ex = np.vstack((U1[t,:],U2[t,:],U3[t,:]))#就是w
    
    umod = 2
    if umod == 1:
        ex[:,0] = ex[:,1]
        ex[:,nmax-1] = ex[:,nmax-2]
    elif umod == 2:
        ex[:,0] = ex[:,1]
        ex[:,nmax-1] = ex[:,nmax-2]
        if ex[1,0] < 0:
            ex[1,0] = -ex[1,1]
        if ex[1,nmax-1] > 0:
            ex[1,nmax-1] = -ex[1,nmax-2]
    
    dw = np.zeros((3,nmax))
    wR = np.zeros((3,nmax))
    wL = np.zeros((3,nmax))
    flux = np.zeros((3,nmax))
    res = np.zeros((3,nmax))
    
    for i in range(0,3):
        for j in range(1,nmax-1):
            dwR = 2*(ex[i,j+1] - ex[i,j])
            dwL = 2*(ex[i,j] - ex[i,j-1])
            dwC = (ex[i,j+1] - ex[i,j-1])/2
            dw[i,j] = minmod(dwR,dwL,dwC)
    
    for i in range(1,nmax-1):
        wL[:,i] = ex[:,i]+ dw[:,i]/2
        wR[:,i] = ex[:,i+1] - dw[:,i+1]/2
    
    #wR[:,0] = wL[:,0] = wL[:,1] - dw[:,1]/2
    #wR[:,nmax-1]  = wL[:,nmax-1] = wL[:,nmax-2] + dw[:,nmax-2]/2
    mod = 1
    wR[:,0] = ex[:,1] - dw[:,1]/2
    wL[:,nmax-1] = ex[:,nmax-2] + dw[:,nmax-2]/2
    if mod == 1:
        wL[:,0] = wL[:,1]
        wR[:,nmax-1]  = wR[:,nmax-2]
    elif mod == 2:
        wL[0,0] = wR[0,1]
        wL[1,0] = -wR[1,1]
        wL[2,0] = wR[2,1]
        wR[0,nmax-1] = wR[0,nmax-2]
        wR[1,nmax-1] = -wR[1,nmax-1]
        wR[2,nmax-1] = wR[1,nmax-1]
    for i in range(0,nmax):
        flux[:,i] = LFflux(wL[:,i],wR[:,i])
        if(i > 0):
            res[:,i] = res[:,i] + flux[:,i]/dx
        if(i < nmax-1):
            res[:,i+1] = res[:,i+1] - flux[:,i]/dx
    
    return res

tstart = 0
tend = 10
import matplotlib.pyplot as plt
coordx = np.zeros((nmax))
for i in range(nmax):
    coordx[i] = i
for t in range(0,tmax-2):
    
    for i in range(0,nmax):
        rho[t,i] = U1[t,i]
        u[t,i] = U2[t,i]/rho[t,i]
        E[t,i] = U3[t,i]/rho[t,i]
        p[t,i] = (gamma-1)*(E[t,i] - u[t,i]*u[t,i]/2)*rho[t,i]
    
    dt = cfl*dx/max(abs(u[t,:]) + np.sqrt(gamma*p[t,:]/rho[t,:]))
    tstart += dt
    if tstart > tend:
        break
    
    Utemp[0,:] = U1[t,:]
    Utemp[1,:] = U2[t,:]
    Utemp[2,:] = U3[t,:]
    
    Res = MUSCL()
    U1[t,:] -= dt*Res[0,:]
    U2[t,:] -= dt*Res[1,:]
    U3[t,:] -= dt*Res[2,:]
    
    U1[t,0] = U1[t,1]
    U1[t,nmax-1] = U1[t,nmax-2]
    U2[t,0] = U1[t,0]#U2[t,1]
    U2[t,nmax-1] = U2[t,nmax-2]
    U3[t,0] = U3[t,1]
    U3[t,nmax-1] = U3[t,nmax-2]
    
    Res = MUSCL()
    U1[t+1,:] = (Utemp[0,:] + U1[t,:] - dt*Res[0,:])/2
    U2[t+1,:] = (Utemp[1,:] + U2[t,:] - dt*Res[1,:])/2
    U3[t+1,:] = (Utemp[2,:] + U3[t,:] - dt*Res[2,:])/2
    
    U1[t+1,0] = U1[t+1,1]
    U1[t+1,nmax-1] = U1[t+1,nmax-2]
    U2[t+1,0] = U1[t+1,0]#U2[t+1,1]
    U2[t+1,nmax-1] = U2[t+1,nmax-2]
    if U2[t+1,0] < 0:
        U2[t+1,0] *= (-1)
    if U2[t+1,nmax-1] > 0:
        U2[t+1,nmax-1] *= (-1)
    U3[t+1,0] = U3[t+1,1]
    U3[t+1,nmax-1] = U3[t+1,nmax-2]
    
    U1[t,:] = Utemp[0,:]
    U2[t,:] = Utemp[1,:]
    U3[t,:] = Utemp[2,:]
    plt.plot(coordx, U1[t,:], 's',label='original values')
    plt.pause(0.1)
    
    
    