import numpy as np
# Bi Conju
A = np.array([[1,2,3],[4,5,6],[7,8,9]],dtype = float)
M = np.array([[1,0,0],[0,1,0],[0,0,1]],dtype = float)
b = np.array([14,32,50],dtype = float)
nmax = 3
x = np.zeros((nmax))
p = np.zeros((nmax))
ptld = np.zeros((nmax))
q = np.zeros((nmax))
qtld = np.zeros((nmax))

res = b.copy() # 残差
for i in range(0,nmax):
    for j in range(0,nmax):
        res[i] -= A[i,j]*x[j]
restld = res.copy()
tmax = 100


rho1 = 1
for t in range(0,tmax):
    # z = M / r # 这玩意咋除？
    z = res.copy()
    # z = M' / restld # 这玩意咋除？
    ztld = restld.copy()
    rho = 0
    for i in range(0,nmax):
        rho += z[i]*restld[i]
    if rho < 1e-8:
        break
    
    beta = 0
    if t > 0:
        beta = rho / rho1
    p = z + beta * p
    ptld = ztld + beta * ptld
    
    for i in range(0,nmax):
        for j in range(0,nmax):
            q[i] += A[i,j]*p[j]
            qtld[i] += A[j,i]*ptld[j]
            
    alpha = 0
    for i in range(0,nmax):
        alpha += ptld[i]*q[i]
    alpha = rho / alpha
            
    x = x + alpha * p
    res = res - alpha * q
    restld = restld - alpha*qtld
    
    rho1 = rho
    