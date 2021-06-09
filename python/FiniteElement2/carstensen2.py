import numpy as np
'''
LinearElascity

参考 https://www.math.hu-berlin.de/~cc/cc_homepage/software/software.shtml

Matlab implementation of the finite element method in elasticity

代码地址：https://www.math.hu-berlin.de/~cc/cc_homepage/software/code/2002-AJ_CC_FS_KR-Matlab_Implementation_FEM_Elasticity.tar.gz

本地代码：D:\FluidSim\FluidSim\FEMNEW\2002-AJ_CC_FS_KR-Matlab_Implementation_FEM_Elasticity\Software2\fem_lame2d\cooks

代码未完工
'''

E = 2900
nu = 0.4
mu = E / (2 * (1 + nu))
lamb = E * nu /((1 + nu)*(1 - 2*nu))
element = np.loadtxt('datafile/elements2.dat')
neumann = np.loadtxt('datafile/neumann2.dat')
coordinates = np.loadtxt('datafile/coordinates2.dat')
dirichlet = np.loadtxt('datafile/dirichlet2.dat')
nmax = 2 * coordinates.shape[0]
A = np.zeros((nmax,nmax))
b = np.zeros((nmax))
idx = np.zeros((6),dtype = int)
for i in range(element.shape[0]):
    ele0 = int(element[i,0] - 1)
    ele1 = int(element[i,1] - 1)
    ele2 = int(element[i,2] - 1)
    phi = np.ones((3,3))
    phi[1:3,0] = coordinates[ele0,:]
    phi[1:3,1] = coordinates[ele1,:]
    phi[1:3,2] = coordinates[ele2,:]
    pinv = np.linalg.inv(phi)
    R = np.zeros((3,6))
    R[0,0] = R[2,1] = pinv[0,1]
    R[0,2] = R[2,3] = pinv[1,1]
    R[0,4] = R[2,5] = pinv[2,1]
    R[2,0] = R[1,1] = pinv[0,2]
    R[2,2] = R[1,3] = pinv[1,2]
    R[2,4] = R[1,5] = pinv[2,2]
    C1 = mu*np.array([[2,0,0],[0,2,0],[0,0,1]])
    C2 = lamb*np.array([[1,1,0],[1,1,0],[0,0,0]])
    C = C1 + C2
    term1 = np.zeros((6,3))
    for j in range(0,6):
        for k in range(0,3):
            term1[j,:] += R[k,j]*C[k,:]
    term2 = np.zeros((6,6))
    for j in range(0,6):
        for k in range(0,3):
            term2[j,:] += term1[j,k]*R[k,:]    
    stima = np.linalg.det(phi)*term2/2
    idx[0] = 2 * element[i,0] - 2
    idx[1] = idx[0] + 1
    idx[2] = 2 * element[i,1] - 2
    idx[3] = idx[2] + 1
    idx[4] = 2 * element[i,2] - 2
    idx[5] = idx[4] + 1
    for j in range(0,6):
        for k in range(0,6):
            A[idx[j],idx[k]] += stima[j,k]
        
tempn = np.zeros((neumann.shape[0],2))
for i in range(neumann.shape[0]):
    x0 = int(neumann[i,0] - 1)
    x1 = int(neumann[i,1] - 1)
    tempn[i,1] = -(coordinates[x1,0] - coordinates[x0,0])
    tempn[i,0] = coordinates[x1,1] - coordinates[x0,1]
    idx[0] = 2 * neumann[i,0] - 2
    idx[1] = idx[0] + 1
    idx[2] = 2 * neumann[i,1] - 2
    idx[3] = idx[2] + 1
    term = (coordinates[x1,:] + coordinates[x0,:])/2
    norm = np.sqrt(tempn[i,0]**2 + tempn[i,1]**2)
    term1 = tempn[i,:] / norm
    
W = np.zeros((34))
M = np.zeros((34,2))
B = np.zeros((34,578))
for i in range(34):
    if (i % 2) == 0:
        M[i,0] = 1
    else:
        M[i,1] = 1
Dirich = np.array([0,2,6,12,20,54,60,72,90,170,178,188,196,230,238,298,304]);
for i in range(0,2):
    for j in range(0,2):
        for k in range(0,34):
            x0 = int(j + 2*k)
            if x0 >= 34:
                continue
            y0 = int(Dirich[k] + i) 
            B[x0,y0] = M[x0,i]
        tes = 1
A2 = np.zeros((612,612))
A2[0:578,0:578] = A[:,:]
A2[0:578,578:612] = np.transpose((B))
A2[578:612,0:578] = B
b2 = np.zeros((612))
b2[0:578] = b.copy()
b2[578:612] = W[:]