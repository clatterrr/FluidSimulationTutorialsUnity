import numpy as np
# CGLS
A = np.array([[1,2,3],[4,5,6],[7,8,9]],dtype = float)
b = np.array([14,32,50],dtype = float)
nmax = 3
x = np.zeros((nmax)) # 待求解的值
res = b - np.dot(A,x) # 残差
p = np.dot(np.transpose(A),res)
s = p.copy()
gamma = np.linalg.norm(s)**2
tmax = 100
xt = np.zeros((tmax,nmax))
for t in range(0,tmax):
    Ap = np.dot(A,p)
    alpha  = gamma / np.linalg.norm(Ap)**2
    x = x + alpha * p
    xt[t,:] = x
    res = res - alpha * Ap
    s = np.dot(np.transpose(A),res)
    gamma1 = np.linalg.norm(s)**2
    beta = gamma1 / gamma
    gamma = gamma1
    p = s + beta * p
    if gamma < 1e-8:
        break
    