import numpy as np
# http://www.netlib.org/templates/matlab/
A = np.array([[1,2,3],[4,5,6],[7,8,9]],dtype = float)
M = np.array([[1,0,0],[0,1,0],[0,0,1]],dtype = float)
Minv = np.linalg.inv(M)
b = np.array([14,32,50],dtype = float)
nmax = 3
x = np.zeros((nmax))
u = np.zeros((nmax))
p = np.zeros((nmax))
q = np.zeros((nmax))

res = np.dot(Minv,b - np.dot(A,x)) # 残差
bnrm2 = np.linalg.norm(b)
tmax = 100
restld = res.copy()
rho1 = 1
for t in range(0,tmax):
    rho = np.dot(restld,res)
    if rho == 0:
        break
    beta = 0
    if t > 0:
        beta = rho / rho1
    u = res + beta * q
    p = u + beta*(q + beta * p)
    
    phat = np.dot(Minv,p)
    vhat = np.dot(A,phat)
    alpha = rho / (np.dot(restld,vhat))
    q = u - alpha*vhat
    uhat = np.dot(Minv,u + q)
    x = x + alpha * uhat
    res = res - alpha * np.dot(A,uhat)
    error = np.linalg.norm(res) / bnrm2
    if error < 1e-8:
        break
    rho1 = rho
    
    
    