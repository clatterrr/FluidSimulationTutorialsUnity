import numpy as np
'''
Immersed Boundary Method for Simulating
Interfacial Problems by Wanho Lee  and Seunggyu Lee 

An Immersed Boundary Method for
a Contractile Elastic Ring in a ThreeDimensional Newtonian Fluid
Seunggyu Lee, Darae Jeong, Wanho Lee
& Junseok Kim
'''
L = 1
N = 64
h = L / N
ds = 0.5 * h
Re = 100
We = 0.1
tmax = 0.1
dt = 0.0005
itermax = int(tmax / dt)

ip = np.zeros((N),dtype = int)
im = np.zeros((N),dtype = int)
for i in range(0,N):
    ip[i] = i + 1
    im[i] = i - 1
ip[N-1] = 0
im[0] = N - 1
kmax = 10000
temp = np.zeros((kmax,2))
for i in range(0,kmax):
    theta= i * 2 * np.pi / kmax
    temp[i,0] = L / 2 + L / 4 * np.cos(theta)
    temp[i,1] = L / 2 + L / 8 * np.sin(theta)
length = np.zeros((kmax))
Nb = 0
Xt = np.zeros((kmax,2))
Xt[0,0] = temp[0,0]
Xt[0,1] = temp[0,1]
# We dynamically reduce the number of Lagrangian boundary points when
# the distance between adjacent points is too small.
for i in range(0,kmax-1):
    length[i+1] = length[i] + np.sqrt((temp[i+1,0] - temp[i,0])**2 + (temp[i+1,1] - temp[i,1])**2)
    if (length[i+1] > (Nb + 1)*ds) & (length[i] < (Nb + 1)*ds):
        Nb = Nb + 1
        Xt[Nb,0] = temp[i,0]
        Xt[Nb,1] = temp[i,1]
Nb = Nb + 1
X = Xt[0:Nb,:]
kp = np.zeros((Nb),dtype = int)
km = np.zeros((Nb),dtype = int)
for i in range(0,Nb):
    kp[i] = i + 1
    km[i] = i - 1
kp[Nb-1] = 0
km[0] = Nb - 1
u = np.zeros((N,N))
v = np.zeros((N,N))
p = np.zeros((N,N))
xgrid = np.zeros((N,N))
ygrid = np.zeros((N,N))
for i in range(0,N):
    xgrid[i,:] = i * h
    ygrid[:,i] = i * h
a = np.zeros((2,N,N),dtype = complex)
b = np.zeros((2,N,N),dtype = complex)
Amat = np.zeros((N,N))
for i in range(0,N):
    for j in range(0,N):
        if (((i == 0)|(i == N//2))&((j == 0)|(j == N//2))):
            test = 1
        else:
            t1 = 2 * np.pi / N * i
            t2 = 2 * np.pi / N * j
            s1 = np.sin(t1)
            s2 = np.sin(t2)
            ss1 = s1 / (s1*s1 + s2*s2)
            ss2 = s2 / (s1*s1 + s2*s2)
            a[0,i,j] = -1j*dt/h*s1
            a[1,i,j] = -1j*dt/h*s2
            b[0,i,j] = -1j*ss1*h/dt
            b[1,i,j] = -1j*ss2*h/dt
        s1 = np.sin(np.pi / N * i)
        s2 = np.sin(np.pi / N * j)
        Amat[i,j] = 1 + (4*dt)/(Re*h*h)*(s1*s1 + s2*s2)
            
# setting the Dirac-delta function along the x-direction
def phi1(r):
    w = np.zeros((4,4))
    s = 1 - r
    w[0,:]=(5-2*abs(-2+s)-np.sqrt(-7+12*abs(-2+s)-4*(-2+s)**2))/8;
    w[1,:]=(3-2*abs(-1+s)+np.sqrt(1+4*abs(-1+s)-4*(-1+s)**2))/8;
    w[2,:]=(3-2*abs(s)+np.sqrt(1+4*abs(s)-4*s**2))/8;
    w[3,:]=(5-2*abs(1+s)-np.sqrt(-7+12*abs(1+s)-4*(1+s)**2))/8;  
    return w        

# setting the Dirac-delta function along the y-direction
def phi2(r):
    w = np.zeros((4,4))
    s = 1 - r
    w[:,0]=(5-2*abs(-2+s)-np.sqrt(-7+12*abs(-2+s)-4*(-2+s)**2))/8;
    w[:,1]=(3-2*abs(-1+s)+np.sqrt(1+4*abs(-1+s)-4*(-1+s)**2))/8;
    w[:,2]=(3-2*abs(s)+np.sqrt(1+4*abs(s)-4*s**2))/8;
    w[:,3]=(5-2*abs(1+s)-np.sqrt(-7+12*abs(1+s)-4*(1+s)**2))/8;  
    return w   

Fx = np.zeros((Nb))
Fy = np.zeros((Nb))
for t in range(0,2):
    # Next, the elastic boundary force density is evaluated from the IB structure. It can be obtained from
    # a minus gradient of the elastic energy functional with respect to X
    for i in range(0,Nb):
        Fx[i] = (X[kp[i],0] + X[km[i],0] - 2*X[i,0])/ds/ds
        Fy[i] = (X[kp[i],1] + X[km[i],1] - 2*X[i,1])/ds/ds
    
    # Spread Start
    # spreading out the boundary force density to the fluid force density using the
    # Dirac-delta function
    c = ds / h / h
    f = np.zeros((2,N,N))
    for k in range(0,Nb):
        s = X[k,:] / h
        npf = np.floor(s)
        r = s - npf
        w = phi1(r[0])*phi2(r[1])
        for i in range(0,4):
            for j in range(0,4):
                i1 = int(np.mod(npf[0] + i - 1 + N,N))
                i2 = int(np.mod(npf[1] + j - 1 + N,N))
                f[0,i1,i2] += (c*Fx[k])*w[i,j]
                f[1,i1,i2] += (c*Fy[k])*w[i,j]
    
    # Spread End
    # NSsolver Start
    uskew = np.zeros((N,N))
    vskew = np.zeros((N,N))
    for i in range(0,N):
        for j in range(0,N):
            uskew[i,j] = (u[ip[i],j] + u[i,j])*u[ip[i],j] - (u[im[i],j]+u[i,j])*u[im[i],j]
            uskew[i,j] += (v[i,ip[j]] + v[i,j])*u[i,ip[j]] - (u[i,im[j]]+u[i,j])*u[im[i],j]
            vskew[i,j] = (u[ip[i],j] + u[i,j])*v[ip[i],j] - (u[im[i],j]+u[i,j])*v[im[i],j]
            vskew[i,j] += (v[i,ip[j]] + v[i,j])*v[i,ip[j]] - (u[i,im[j]]+u[i,j])*v[im[i],j]
    w1 = u - 0.5*dt*uskew + (dt/We)*f[0,:,:]
    w2 = u - 0.5*dt*uskew + (dt/We)*f[1,:,:]
    w1 = np.fft.fft2(w1)
    w2 = np.fft.fft2(w2)
    p = b[0,:,:]*w1 + b[1,:,:]*w2
    u = (w1 + a[0,:,:]*p[:,:])/Amat
    v = (w2 + a[1,:,:]*p[:,:])/Amat
    
    # p = np.fft.ifft2(p,1)
    p = np.fft.ifft2(p).real
    u = np.fft.ifft2(u).real
    v = np.fft.ifft2(v).real
    #NS solver End
    
    U = np.zeros((Nb,2))
    for k in range(0,Nb):
        s = X[k,:] / h
        npf = np.floor(s)
        r = s - npf
        # Subsequently, the Dirac-delta function is used again to obtain
        # the velocity on the elastic IB boundary
        w = phi1(r[0])*phi2(r[1])
        for i in range(0,4):
            for j in range(0,4):
                i1 = int(np.mod(npf[0] + i - 1 + N,N))
                i2 = int(np.mod(npf[1] + j - 1 + N,N))
                U[k,0] += w[i,j]*u[i1,i2]
                U[k,1] += w[i,j]*v[i1,i2]
    X = X + dt * U
    test = 1
    
    
    
    