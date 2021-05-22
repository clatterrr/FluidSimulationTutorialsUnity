import numpy as np
#5-th WENO for ShockTube成功
tmax = 1000
nmax = 20



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

F1 = np.zeros((tmax,nmax))
F2 = np.zeros((tmax,nmax))
F3 = np.zeros((tmax,nmax))


Bze = np.zeros((3,nmax))
Bon = np.zeros((3,nmax))
Btw = np.zeros((3,nmax))
wn = np.zeros((3,nmax))
wp = np.zeros((3,nmax))

Utemp = np.zeros((3,nmax))
Res = np.zeros((3,nmax))

gamma = 1.4
dt = 0.042
dx = 1/(nmax-1)
cfl = 0.5

for i in range(0,nmax):
    if(i < 10):
        U1[0,i] = 1
        U3[0,i] = 2.5
    else:
        U1[0,i] = 0.125
        U3[0,i] = 0.25

def LFflux(qL,qR):
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

def WENO5():
    rho2 = np.zeros((nmax))
    u2 = np.zeros((nmax))
    E2 = np.zeros((nmax))
    p2 = np.zeros((nmax))
    for i in range(0,nmax):
        rho2[i] = U1[t,i]
        u2[i] = U2[t,i]/rho2[i]
        E2[i] = U3[t,i]/rho2[i]
        p2[i] = (gamma-1)*(E2[i] - u2[i]*u2[i]/2)*rho2[i]
    ex = np.vstack((rho2[:],u2[:],p2[:]))
    qL = np.zeros((3,nmax))
    qR = np.zeros((3,nmax))
    flux = np.zeros((3,nmax))
    res = np.zeros((3,nmax))
    #WENO5虽然使用5个点，但只有三个sub-stencils
    #第一个为x[i-2],x[i-1],x[i]
    #第二个为x[i-1],x[i],x[i+1]
    #第三个为x[i],x[i+1],x[i+2]
    
    #The suggested Smooth Indicatores of Jiang and Shu，也就是WENO-JS
    for i in range(2,nmax-2):
        Bze[:,i] = 13/12*(ex[:,i-2]-2*ex[:,i-1]+ex[:,i])**2 + 0.25*(ex[:,i-2]-4*ex[:,i-1]+3*ex[:,i])**2
        Bon[:,i] = 13/12*(ex[:,i-1]-2*ex[:,i]+ex[:,i+1])**2 + 0.25*(ex[:,i-1]-ex[:,i+1])**2
        Btw[:,i] = 13/12*(ex[:,i]-2*ex[:,i+1]+ex[:,i+2])**2 + 0.25*(3*ex[:,i]-4*ex[:,i+1]+ex[:,i+2])**2
    
    eps = 1e-6
    #eps = 1e-6是经验值，(eps + Beta)^p，其中p=2也是经验值
    
    #在spyder的变量查看器中,A0只有0~nmax-4个格子是有效的
    #A0只有1~nmax-3个格子是有效的，A0只有2~nmax-2个格子是有效的
    A0 = 1/10/(eps+Bze)/(eps+Bze)
    A1 = 6/10/(eps+Bon)/(eps+Bon)
    A2 = 3/10/(eps+Btw)/(eps+Btw)
    ASum = A0 + A1 + A2
    
    #ENO stencils weights
    for i in range(2,nmax-2):
        wn[:,i-1] = A0[:,i]/ASum[:,i]*(2*ex[:,i-2]-7*ex[:,i-1]+11*ex[:,i])/6 
        wn[:,i-1] += A1[:,i]/ASum[:,i]*(-ex[:,i-1]+5*ex[:,i]+2*ex[:,i+1])/6 
        wn[:,i-1] += A2[:,i]/ASum[:,i]*(2*ex[:,i]+5*ex[:,i+1]-ex[:,i+2])/6 
    
    A0 = 3/10/(eps+Bze)/(eps+Bze)
    A1 = 6/10/(eps+Bon)/(eps+Bon)
    A2 = 1/10/(eps+Btw)/(eps+Btw)
    ASum = A0 + A1 + A2
    
    for i in range(2,nmax-2):
        wp[:,i] = A0[:,i]/ASum[:,i]*(-ex[:,i-2]+5*ex[:,i-1]+2*ex[:,i])/6 
        wp[:,i] += A1[:,i]/ASum[:,i]*(2*ex[:,i-1]+5*ex[:,i]-ex[:,i+1])/6 
        wp[:,i] += A2[:,i]/ASum[:,i]*(11*ex[:,i]-7*ex[:,i+1]+2*ex[:,i+2])/6 
        
    wn[:,0] = wn[:,1] = wn[:,2] = ex[:,2]
    wn[:,nmax-1] = wn[:,nmax-2] = wn[:,nmax-3] = ex[:,nmax-3]
    wp[:,0] = wp[:,1] = wp[:,2] = ex[:,2]
    wp[:,nmax-1] = wp[:,nmax-2] = wp[:,nmax-3] = ex[:,nmax-3]
     
    qR[0,:] = wn[0,:]
    qR[1,:] = wn[1,:]*wn[0,:]
    qR[2,:] = wn[2,:]/(gamma-1) + 0.5*wn[0,:]*wn[1,:]*wn[1,:]
    
    qL[0,:] = wp[0,:]
    qL[1,:] = wp[1,:]*wp[0,:]
    qL[2,:] = wp[2,:]/(gamma-1) + 0.5*wp[0,:]*wp[1,:]*wp[1,:]   
    
    for i in range(0,nmax):
        flux[:,i] = LFflux(qL[:,i],qR[:,i])
        if(i > 0):
            res[:,i] = res[:,i] + flux[:,i]/dx
        if(i < nmax-1):
            res[:,i+1] = res[:,i+1] - flux[:,i]/dx
    
    return res

tstart = 0
tend = 0.15
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
    
    #q = qo-dt*dF; 
    #q = 0.75*qo+0.25*(q-dt*dF);
    #q = (qo+2*(q-dt*dF))/3;
    
    Res = WENO5()
    U1[t,:] = Utemp[0,:] - dt*Res[0,:]
    U2[t,:] = Utemp[1,:] - dt*Res[1,:]
    U3[t,:] = Utemp[2,:] - dt*Res[2,:]
    
    Res = WENO5()
    U1[t,:] = 0.75*Utemp[0,:] + 0.25*(U1[t,:] -  dt*Res[0,:])
    U2[t,:] = 0.75*Utemp[1,:] + 0.25*(U2[t,:] -  dt*Res[1,:])
    U3[t,:] = 0.75*Utemp[2,:] + 0.25*(U3[t,:] -  dt*Res[2,:])
    
    Res = WENO5()
    U1[t+1,:] = (Utemp[0,:] + 2*(U1[t,:] - dt*Res[0,:]))/3
    U2[t+1,:] = (Utemp[1,:] + 2*(U2[t,:] - dt*Res[1,:]))/3
    U3[t+1,:] = (Utemp[2,:] + 2*(U3[t,:] - dt*Res[2,:]))/3
    
    U1[t,:] = Utemp[0,:]
    U2[t,:] = Utemp[1,:]
    U3[t,:] = Utemp[2,:]
    
        