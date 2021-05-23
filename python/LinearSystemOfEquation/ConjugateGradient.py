import numpy as np
# http://www.netlib.org/templates/matlab/
A = np.array([[1,2,3],[4,5,6],[7,8,9]],dtype = float)
M = np.array([[1,0,0],[0,1,0],[0,0,1]],dtype = float)
Minv = np.linalg.inv(M)
b = np.array([14,32,50],dtype = float)
nmax = 3
x = np.zeros((nmax))
p = np.zeros((nmax))
res = np.dot(Minv,b - np.dot(A,x)) # 残差
bnrm2 = np.linalg.norm(b)
tmax = 100
restld = res.copy()
rho1 = 1
for t in range(0,tmax):
    z = np.dot(Minv,res) # 预处理方程组的残差
    rho = np.dot(np.transpose(res),z)
    if rho == 0:
        break
    beta = 0
    if t > 0:
        beta = rho / rho1
    p = z + beta*p
    Ap = np.dot(A,p)
    alpha = rho / np.dot(np.transpose(p),Ap)
    x = x + alpha * p
    res = res - alpha * Ap
    error = np.linalg.norm(res) / bnrm2
    if error < 1e-8:
        break
    rho1 = rho
    
    
    