import numpy as np
'''
Based on "Finite Element Methods for flow problems" of
Jean Donea and Antonio Huerta
参考自特别好的非定常代码
https://www.mathworks.com/matlabcentral/fileexchange/?q=profileid:4187051
这份代码是有限元Galerkin法解非定常一维对流扩散方程
暂时不打算继续写下去
如此看来之前不应该用等距数值积分，而应该用高斯数值积分
而且限制条件值得参考
'''
beta = 0
polynomial_degree=1 # Shape functions polynomial degree
dof = polynomial_degree+1 # Number of DOFs per element
nmax = 150
Lei = 0.1
h = Lei / polynomial_degree
J = h / 2 # % Jacobian of the transformation
a = 1 # 对流速度
v = 0.01 / 3 # 扩散速度

GaussN = 3
gxi = np.zeros((GaussN,GaussN))
geta = np.zeros((GaussN,GaussN))
gw = np.zeros((GaussN,GaussN))
gxi[:,0] = -0.7746
gxi[:,2] = 0.7746
geta[0,:] = -0.7746
geta[2,:] = 0.7746
gw[0,:] = gw[2,:] = 5/9
gw[1,:] = 8/9
gw[:,0] *= 5/9
gw[:,2] *= 5/9
gw[:,1] *= 8/9

xgauss = np.zeros((GaussN))
ygauss = np.zeros((GaussN))
for i in range(0,GaussN):
    xgauss[i] = Lei / 2 * (1 + gxi[i])
    ygauss[i] = Lei / 2 * (1 + gxi[i])


Jmat = np.array([Lei / 2,0],[0,Lei / 2])
J = 0.25
# Shape Function
N = np.zeros((2,GaussN))
dN = np.zeros((2,GaussN))
W = np.zeros((2,GaussN))
dW = np.zeros((2,GaussN))
for i in range(0,GaussN):
    N[0,i] = (1 - csi[i])/2
    N[1,i] = (1 + csi[i])/2
    dN[0,i] = -0.5
    dN[1,i] = 0.5
    
    W[0,i] = (1 - csi[i])/2 - 3/4*beta*(1 - csi[i]*csi[i])
    W[1,i] = (1 + csi[i])/2 + 3/4*beta*(1 - csi[i]*csi[i])
    dW[0,i] = -0.5 + 3/2*beta*csi[i]
    dW[1,i] = 0.5 - 3/2*beta*csi[i]

A = np.zeros((nmax,dof))
for i in range(0,nmax):
    for j in range(0,dof):
        A[i,j] = (i - 1) + (dof-1) + j
        
        
Mass = np.zeros((dof,dof))
for i in range(0,dof):
    for j in range(0,dof):
        for k in range(0,GaussN):
            Mass[i,j] = Mass[i,j] + (W[i,k]*N[j,k])*w[k]
        Mass[i,j] = Mass[i,j]*J
        
Convection = np.zeros((dof,dof))
for i in range(0,dof):
    for j in range(0,dof):
        for k in range(0,GaussN):
            Convection[i,j] = Convection[i,j] + (W[i,k]*N[j,k])*w[k]
        Convection[i,j] = Convection[i,j]*a
        
Diffusion = np.zeros((dof,dof))
for i in range(0,dof):
    for j in range(0,dof):
        for k in range(0,GaussN):
            Diffusion[i,j] = Diffusion[i,j] + (W[i,k]*N[j,k])*w[k]
        Diffusion[i,j] = Diffusion[i,j]*v/J
        
        
        
        