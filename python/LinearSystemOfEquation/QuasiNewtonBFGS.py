import numpy as np
# BFGS

def fun(a):
    return a[0]*a[0]/2 + a[1]*a[1] - a[0]*a[1] - 2*a[0]

tmax = 100
nmax = 2
x = np.zeros((tmax,nmax))
Q = np.array(([1,-1],[-1,2]))
b = np.array([-2,0])
x[0,:] = 1
B = np.eye((nmax))
jac = np.zeros((tmax,nmax))
s = np.zeros((nmax))
y = np.zeros((nmax))
for t in range(0,tmax):
    # Scipy中的3-points
    for i in range(0,nmax):
        x0 = x[t,:].copy()
        x1 = x[t,:].copy()
        x0[i] += 1
        x1[i] -= 1
        jac[t,i] = (fun(x0) - fun(x1))/2
        
    if t > 0:
        s = x[t,:] - x[t-1,:]
        y = jac[t,:] - jac[t-1,:]
        norm = np.dot(np.transpose(y),s)
        if norm < 1e-8:
            break
        else:
            term1 = np.zeros((nmax))
            for i in range(0,nmax):
                for j in range(0,nmax):
                    term1[i] += B[i,j]*s[j]
            term2 = np.zeros((nmax,nmax))
            for i in range(0,nmax):
                term2[i,:] = term1[i]*s[:]
            term3 = np.zeros((nmax,nmax))
            for i in range(0,nmax):
                for j in range(0,nmax):
                    term3[i,:] += term2[i,j]*B[j,:]
            
            term4 = np.zeros((nmax))
            for i in range(0,nmax):
                term4[:] += s[i]*B[i,:]
            term5 = np.dot(term4,s)
            
            term6 = np.zeros((nmax,nmax))
            for i in range(0,nmax):
                term6[i,:] = y[i] * y[:]
            term7 = np.dot(y,s)              
            
            
            B = B -  term3 / term5 + term6 / term7
            
    Binv = np.linalg.inv(B)
    d = np.dot(Binv,-jac[t,:])
    num = np.dot(np.transpose(jac[t,:]),d)
    denum = np.dot(np.transpose(d),np.dot(Q,d))
    alpha = - num / denum
    x[t+1,:] = x[t,:] + alpha * d
    
    
    