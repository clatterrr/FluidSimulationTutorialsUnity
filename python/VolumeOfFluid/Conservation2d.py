import numpy as np
nmax = 8
mmax = 8
f = np.zeros((nmax,mmax))

n = np.zeros((nmax,mmax))
fx = np.zeros((nmax,mmax))
fy = np.zeros((nmax,mmax))

tmax = 1
alpha = 2
for t in range(0,tmax):
    
    for i in range(1,nmax-1):
        for j in range(1,mmax-1):
            fw = (f[i-1,j-1] + alpha*f[i-1,j] + f[i-1,j+1])/(2 + alpha)
            fe = (f[i+1,j-1] + alpha*f[i+1,j] + f[i+1,j+1])/(2 + alpha)
            fn = (f[i-1,j+1] + alpha*f[i,j+1] + f[i+1,j+1])/(2 + alpha)
            fs = (f[i-1,j-1] + alpha*f[i,j-1] + f[i+1,j-1])/(2 + alpha)
            fx[i,j] = (fw - fe)/2
            fy[i,j] = (fn - fs)/2
            