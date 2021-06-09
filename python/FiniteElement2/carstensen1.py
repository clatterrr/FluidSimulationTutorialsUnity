import numpy as np
'''
Elastoviscoplasticity

参考 https://www.math.hu-berlin.de/~cc/cc_homepage/software/software.shtml

Elastoviscoplastic finite element analysis in 100 lines of Matlab

代码地址：https://www.math.hu-berlin.de/~cc/cc_homepage/software/code/2002-CC_KR-Elastoviscoplastic_FE_Analysis_in_Matlab.tar.gz

本地代码：D:\FluidSim\FluidSim\FEMNEW\Software3

'''
element = np.loadtxt('datafile/elements1.dat')
neumann = np.loadtxt('datafile/neumann1.dat')
coordinates = np.loadtxt('datafile/coordinates1.dat')
dirichlet = np.loadtxt('datafile/dirichlet1.dat')
N = coordinates.shape[0]
initvector = np.zeros((N * 3))
NJ = element.shape[0]
sigma = np.zeros((NJ,9))
u0 = np.zeros((3 * N,1))

# 材质参数 
lamb = 1.107438169066076e+05
mu = 8.019379844961240e+04
C1 = lamb + 2 * mu/3
sigma_y=450
nu = 0
th = 1

dt = 0.5  # 时间长度
for t in range(0, 20):
    C2 = nu/(nu/(2*mu)+th*dt)
    C3 = th*dt*sigma_y/(nu/(2*mu)+th*dt)
    tr3 = sigma[:,0] + sigma[:,4] + sigma[:,8]
    res = np.zeros((32,9))
    res[:,0] = res[:,4] = res[:,8] = tr3
    dev3 = sigma - res
    e0 = 1/(9*lamb+6*mu)*res + 1/(2*mu)*dev3;
    # 有限元
    list1 = np.zeros((32,12))
    list2 = np.zeros((32,12))
    list1[:,0] = list1[:,1] = list1[:,2] = 3 * element[:,0]
    list1[:,3] = list1[:,4] = list1[:,5] = 3 * element[:,1]
    list1[:,6] = list1[:,7] = list1[:,8] = 3 * element[:,2]
    list1[:,9] = list1[:,10] = list1[:,11]  = 3 * element[:,3]
    
    list2[:,0] = list2[:,3] = list2[:,6] = list2[:,9] = 2
    list2[:,1] = list2[:,4] = list2[:,7] = list2[:,10] = 1
    list2[:,2] = list2[:,5] = list2[:,8] = list2[:,11] = 0
    list1 = list1 - list2

    Dfun = np.zeros((N*3,N*3))
    Q = np.zeros((N*3))
    for i in range(0,32):
        ele0 = int(element[i,0] - 1) # matlab矩阵是从1开始编号的
        ele1 = int(element[i,1] - 1)
        ele2 = int(element[i,2] - 1)
        ele3 = int(element[i,3] - 1)
        Dmat = np.zeros((4,4))
        Dmat[0,:] = 1
        Dmat[1:4,0] = coordinates[ele0,:]
        Dmat[1:4,1] = coordinates[ele1,:] 
        Dmat[1:4,2] = coordinates[ele2,:] 
        Dmat[1:4,3] = coordinates[ele3,:]
        Dinv = np.linalg.inv(Dmat)
        Data1 = np.zeros((12,9))
        Data2 = np.zeros((12,9))
        Data1[0,0:3] = Data1[1,3:6] = Data1[2,6:9] = Dinv[0,1:4]
        Data1[3,0:3] = Data1[4,3:6] = Data1[5,6:9] = Dinv[1,1:4]
        Data1[6,0:3] = Data1[7,3:6] = Data1[8,6:9] = Dinv[2,1:4]
        Data1[9,0:3] = Data1[10,3:6] = Data1[11,6:9] = Dinv[3,1:4]
        
        Data2[0,0] = Data2[1,1] = Data2[2,2] = Dinv[0,1]
        Data2[0,3] = Data2[1,4] = Data2[2,5] = Dinv[0,2]
        Data2[0,6] = Data2[1,7] = Data2[2,8] = Dinv[0,3]
        Data2[3,0] = Data2[4,1] = Data2[5,2] = Dinv[1,1]
        Data2[3,3] = Data2[4,4] = Data2[5,5] = Dinv[1,2]
        Data2[3,6] = Data2[4,7] = Data2[5,8] = Dinv[1,3]
        Data2[6,0] = Data2[7,1] = Data2[8,2] = Dinv[2,1]
        Data2[6,3] = Data2[7,4] = Data2[8,5] = Dinv[2,2]
        Data2[6,6] = Data2[7,7] = Data2[8,8] = Dinv[2,3]
        Data2[9,0] = Data2[10,1] = Data2[11,2] = Dinv[3,1]
        Data2[9,3] = Data2[10,4] = Data2[11,5] = Dinv[3,2]
        Data2[9,6] = Data2[10,7] = Data2[11,8] = Dinv[3,3]
        
        u = np.zeros((12))
        ut = np.zeros((12))
        for j in range((12)):
            ut[j] = initvector[int(list1[i,j] - 1)]
            u[j] = u0[int(list1[i,j] - 1)]
        eps = (Data1 + Data2) / 2
        v = np.dot(np.transpose(ut - u),eps) + e0[i,:]
        T = np.linalg.det(Dmat) / 6
        
        tr3v = v[0] + v[4] + v[8]
        devv = v.copy()
        devv[0] -= tr3v
        devv[4] -= tr3v
        devv[8] -= tr3v
        norm = np.linalg.norm(devv)
        
        tr3eps = eps[:,0] + eps[:,4] + eps[:,8]
        deveps = eps.copy()
        deveps[:,0] -= tr3eps
        deveps[:,4] -= tr3eps
        deveps[:,8] -= tr3eps
        if norm - sigma_y / (2 * mu) > 0:
            C5 = C2 + C3 / norm
            C6 = C3 / norm**3 * np.dot(deveps,np.transpose(devv))
        else:
            C5 = 2 * mu
            C6 = np.zeros((12))
       
        dm1 = np.zeros((12,12))
        dm2 = np.zeros((12,12))
        dm3 = np.zeros((12,12))
        df1 = np.zeros((12))
        df2 = np.zeros((12))
        for j in range(12):
            dm1[j,:] = C1 * tr3eps[j] * tr3eps[:]
            for k in range(9):
                dm2[j,:] += C5 * deveps[j,k] * eps[:,k]
                df2[j] += C5 * devv[k] * eps[j,k]
            df1[j] = C1 * tr3v * tr3eps[j]
             
        for j in range(9):
            dm3[j,:] = devv[j] * eps[:,j]
        
        
        for j in range(12):
            for k in range(12):
                i0 = int(list1[i,j] - 1)
                j0 = int(list1[i,k] - 1)
                Dfun[i0,j0] += T * (dm1[j,k] + dm2[j,k] - dm3[j,k])
            Q[i0] += T * (df1[j] + df2[j])
    
        test = 1
    test = 1
