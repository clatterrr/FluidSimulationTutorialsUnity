import numpy as np

A = np.array([[1,2,3],[4,5,6],[7,8,9]],dtype = float)
b = np.array([14,32,50],dtype = float)
x0 = np.array([1,2,3],dtype = float)
nmax = 3

grad = - b.copy()
for i in range(0,nmax):
    for j in range(0,nmax):
        grad[i] += A[i,j]*x0[j]
direction = - grad.copy()
tmax = 100


for t in range(0,tmax):
    # 计算步长alpha
    numerator = denominator = 0
    for i in range(0,nmax):
        numerator += grad[i]*grad[i]
    arr = np.zeros((nmax))
    for i in range(0,nmax):
        arr[:] += direction[i]*A[i,:]
    for i in range(0,nmax):
        denominator += arr[i]*direction[i]
    if denominator < 1e-8:
        break
    alpha = numerator / denominator
    
    x0 = x0 + alpha*direction
    
    grad2 = grad.copy()
    for i in range(0,nmax):
        for j in range(0,nmax):
            grad2[i] += alpha*A[i,j]*direction[j]
    denominator = numerator
    numerator = 0
    for i in range(0,nmax):
        numerator += grad2[i]*grad2[i]
    beta = numerator / denominator
    
    direction = beta*direction - grad2
    
    grad = grad2
    