import numpy as np
'''
Nodal Discountines Galerkin
https://github.com/tcew/nodal-dg/tree/master/Codes1.1
'''

N = 8
Np = int((N+1)*(N+2)/2)
idx = 0
L1 = np.zeros((Np))
L2 = np.zeros((Np))
L3 = np.zeros((Np))
x = np.zeros((Np))
y = np.zeros((Np))
for i in range(0,N+1):
    for j in range(0,N+2-i-1):
        L1[idx] = i / N
        L3[idx] = j / N
        L2[idx] = 1 - L1[idx] - L3[idx]
        x[idx] = -L2[idx] + L3[idx]
        y[idx] = (-L2[idx] - L3[idx] + 2*L1[idx])/np.sqrt(3)
        idx += 1
