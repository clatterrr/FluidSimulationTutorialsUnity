import numpy as np
'''
Nodal Discountines Galerkin
https://github.com/tcew/nodal-dg/tree/master/Codes1.1
D:\FluidSim\FluidSim\CompressibeNewgood\nodal-dg-master\nodal-dg-master\Codes1.1
'''

N = 8
Np = int((N+1)*(N+2)/2)
idx = 0
L1 = np.zeros((Np))
L2 = np.zeros((Np))
L3 = np.zeros((Np))
x = np.zeros((Np))
y = np.zeros((Np))
blend1 = np.zeros((Np))
blend2 = np.zeros((Np))
blend3 = np.zeros((Np))
warp1 = np.zeros((Np))
warp2 = np.zeros((Np))
warp3 = np.zeros((Np))

def gamma(x):
    inta = np.exp(-1)*((-1)**(x - 1))
    intb = np.exp(-11)*((-11)**(x - 1))
    intc = np.exp(-21)*((-21)**(x - 1))
    return 10 / 3 * (inta + 4 * intb + intc)

def JacobiGL(alpha,beta,Nx):
    localx = np.zeros((Nx+1,1))
    if Nx == 1:
        localx[0] = -1
        localx[1] = 1
    return localx
for i in range(0,N+1):
    for j in range(0,N+2-i-1):
        L1[idx] = i / N
        L3[idx] = j / N
        L2[idx] = 1 - L1[idx] - L3[idx]
        x[idx] = -L2[idx] + L3[idx]
        y[idx] = (-L2[idx] - L3[idx] + 2*L1[idx])/np.sqrt(3)
        blend1 = 4 * L2[idx]*L3[idx]
        blend2 = 4 * L1[idx]*L3[idx]
        blend3 = 4 * L1[idx]*L2[idx]
        
        # Vandermode
        req = np.array([-1,1])
        Nx = 1
        VID  = np.zeros((2,Nx+1))
        for k in range(0,Nx+1):
            PL = np.zeros((1,2))
            alpha = 0
            beta = 0
            gamma0 = 2**(alpha + beta + 1)/(alpha + beta + 1)*gamma(alpha + 1)
            gamma0 = gamma0 / gamma(beta + 1)/gamma(alpha + beta + 1)
            PL[0,:] = 1 / np.sqrt(gamma0)
            
            test = 1
        Veq = np.array([[0.7071,-1.2247],[0.7071,1.2247]])
        idx += 1
nmax = 100
Q = np.zeros((nmax,4))
gam = 1.4
rho = gam * np.ones((nmax))
p = np.ones((nmax))
rhou = rho * 3
rhov = np.zeros((nmax))
Ener = p / (gam - 1) + (rhou**2 + rhov**2)/(2*rho)
Q[:,0] = rho.copy()
Q[:,1] = rhou.copy()
Q[:,2] = rhov.copy()
Q[:,3] = Ener.copy()

  # % Build weights used in limiting
  # g1 = (dVdxC1.^2 + dVdyC1.^2); g2 = (dVdxC2.^2 + dVdyC2.^2); g3 = (dVdxC3.^2 + dVdyC3.^2);
  
  # epse = 1e-10; fac = g1.^2 + g2.^2 + g3.^2;
  # w1 = (g2.*g3+epse)./(fac+3*epse);
  # w2 = (g1.*g3+epse)./(fac+3*epse);
  # w3 = (g1.*g2+epse)./(fac+3*epse);
  
  # % Limit gradients
  # LdVdxC0 = w1.*dVdxC1 + w2.*dVdxC2 + w3.*dVdxC3;
  # LdVdyC0 = w1.*dVdyC1 + w2.*dVdyC2 + w3.*dVdyC3;