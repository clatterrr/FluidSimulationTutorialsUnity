import numpy as np
v0 = np.array([1.508870,-1.531271,25.46091])
sigma = 10
beta = 8 / 3
rho = 28
H = np.array([1,0,0])
imQ = 10 # Variance of the model error 
R = 1 # Variance of the data error
Qinit = 100 # Variance of the initial condition
ncycle = 100
Dt = 1
nbv = 1000
which = 1
