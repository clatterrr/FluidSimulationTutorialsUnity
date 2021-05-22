import numpy as np
#Lax-Friedriches For 1d Shock Tube
tmax = 1000
nmax = 12



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

F1 = np.zeros((tmax,nmax))#rho
F2 = np.zeros((tmax,nmax))#rho*u
F3 = np.zeros((tmax,nmax))#rho*E

f1 = np.zeros((tmax,nmax+1))
f2 = np.zeros((tmax,nmax+1))
f3 = np.zeros((tmax,nmax+1))

rho_bar = np.zeros((tmax,nmax))#Roe平均
u_bar = np.zeros((tmax,nmax))#Roe平均
H_bar = np.zeros((tmax,nmax))#Roe平均
c_bar = np.zeros((tmax,nmax))#Roe平均

lambda_1 = np.zeros((tmax,nmax))#矩阵A的特征值
lambda_2 = np.zeros((tmax,nmax))#矩阵A的特征值
lambda_3 = np.zeros((tmax,nmax))#矩阵A的特征值

e1_1 = np.zeros((tmax,nmax))#矩阵A的特征向量e1
e1_2 = np.zeros((tmax,nmax))#矩阵A的特征向量e1
e1_3 = np.zeros((tmax,nmax))#矩阵A的特征向量e1
e2_1 = np.zeros((tmax,nmax))#矩阵A的特征向量e2
e2_2 = np.zeros((tmax,nmax))#矩阵A的特征向量e2
e2_3 = np.zeros((tmax,nmax))#矩阵A的特征向量e2
e3_1 = np.zeros((tmax,nmax))#矩阵A的特征向量e3
e3_2 = np.zeros((tmax,nmax))#矩阵A的特征向量e3
e3_3 = np.zeros((tmax,nmax))#矩阵A的特征向量e3

gamma = 1.4
dt = 1
dx = 0.2
cfl = 0.9

for i in range(0,nmax):
    if(i < nmax/2):
        rho[0,i] = 1
        p[0,i] = 1
    else:
        rho[0,i] = 0.125
        p[0,i] = 0.1

for i in range(0,nmax):
    U1[0,i] = rho[0,i]
    U2[0,i] = rho[0,i]*u[0,i]
    #E = p ./ ((gamma-1) * rho) + u.^2/2;
    E[0,i] = p[0,i]/((gamma-1)*rho[0,i]) + u[0,i]*u[0,i]/2
    U3[0,i] = rho[0,i]*E[0,i]

for t in range(0,tmax-2):
    
    
    for i in range(0,nmax):
        rho[t,i] = U1[t,i]
        u[t,i] = U2[t,i]/rho[t,i]
        E[t,i] = U3[t,i]/rho[t,i]
        p[t,i] = (gamma-1)*(E[t,i] - u[t,i]*u[t,i]/2)*rho[t,i]
        e[t,i] = p[t,i]/(gamma-1)/rho[t,i]
        c[t,i] = np.sqrt(gamma * p[t,i]/rho[t,i])
        H[t,i] = (p[t,i] * gamma)/((gamma-1) * rho[t,i]) + u[t,i]*u[t,i]/2
        
        F1[t,i] = rho[t,i]*u[t,i]
        F2[t,i] = rho[t,i]*u[t,i]*u[t,i]+p[t,i]
        F3[t,i] = rho[t,i]*u[t,i]*H[t,i]
        
        #a = sqrt(gamma * p ./ rho);
        #S_max = max(max(abs(u) + a));
    
    smax = max(abs(u[t,:] + np.sqrt(gamma * p[t,:]/rho[t,:])))
    dt = dx * cfl/smax
    
    for i in range(0,nmax-1):
        #Roe平均
        rho_bar[t,i] = (np.sqrt(rho[t,i]) + np.sqrt(rho[t,i+1]))/2
        u_bar[t,i] = (np.sqrt(rho[t,i])*u[t,i] + np.sqrt(rho[t,i+1])*u[t,i+1])/(np.sqrt(rho[t,i]) + np.sqrt(rho[t,i+1]))
        H_bar[t,i] = (np.sqrt(rho[t,i])*H[t,i] + np.sqrt(rho[t,i+1])*H[t,i+1])/(np.sqrt(rho[t,i]) + np.sqrt(rho[t,i+1]))
        c_bar[t,i] = np.sqrt((gamma-1)*(H_bar[t,i] - u_bar[t,i]*u_bar[t,i]/2))
        
        #矩阵A的特征值
        lambda_1[t,i] = u_bar[t,i] - c_bar[t,i]
        lambda_2[t,i] = u_bar[t,i] 
        lambda_3[t,i] = u_bar[t,i] + c_bar[t,i]
        
        #矩阵A的特征向量
        e1_1[t,i] = 1
        e1_2[t,i] = u_bar[t,i] - c_bar[t,i]
        e1_3[t,i] = H_bar[t,i] - u_bar[t,i]*c_bar[t,i]
        e2_1[t,i] = 1
        e2_2[t,i] = u_bar[t,i]
        e2_3[t,i] = u_bar[t,i]*u_bar[t,i]/2   
        e3_1[t,i] = 1
        e3_2[t,i] = u_bar[t,i] + c_bar[t,i]
        e3_3[t,i] = H_bar[t,i] + u_bar[t,i]*c_bar[t,i]
        
        
    # RF = 0.5 * (F(:,1:N+2)+F(:,2:N+3))...
    #     - (abs(lambda_bar_1).*alpha_bar_1.*K_bar_1...
    #     + abs(lambda_bar_2).*alpha_bar_2.*K_bar_2...
    #     + abs(lambda_bar_3).*alpha_bar_3.*K_bar_3) * 0.5; %flux acorss cells
    for i in range(0,nmax):
        U1[t+1,i] = U1[t,i] - dt/dx * (f1[t,i+1] - f1[t,i])
        U2[t+1,i] = U2[t,i] - dt/dx * (f2[t,i+1] - f2[t,i])
        U3[t+1,i] = U3[t,i] - dt/dx * (f3[t,i+1] - f3[t,i])