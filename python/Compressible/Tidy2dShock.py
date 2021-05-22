import numpy as np
# 二维欧拉可压缩MacCormack模拟斜激波
# https://github.com/vasko6d/finite-volume-solver
nmax = 82
mmax = 162
tmax = 4

bigx = np.zeros((nmax+1,mmax+1))
bigy = np.zeros((nmax+1,mmax+1))
gridres = 4
IS = 5*gridres + 1
H = 1
L = 3.2
theta = 10.94
dx = L/(mmax - IS)
dy = np.zeros((mmax))
dy2 = np.zeros((mmax + 1))
'''
     dx
---------------
\\\)10.94     \
   \\\        \
      \\\     \dh
         \\\  \
            \\\

'''
dh = np.tan(theta*np.pi/180)*dx
dy[0:IS-1] = H/(nmax-2)
for i in range(IS-1,mmax):
    dy[i] = (H-dh*(i-IS+1.5))/(nmax-2)
for j in range(0,mmax+1):
    bigx[:,j] = j*dx
for i in range(0,mmax + 1):
    if i < 21:
        dy2[i] = H/(nmax - 2)
    else:
        dy2[i] = (H - dh*(i + 1 - IS))/(nmax - 2)
for j in range(0,mmax+1):
    for i in range(0,nmax+1):
        bigy[i,j] = dy2[j]*(i-1)
dis = 0.85
dymin = bigy[1,mmax-1] - bigy[0,mmax-1]
p1 = 100000
rho1 = 1
M1 = 2.9
gamma = 1.4
c1 = np.sqrt(gamma*p1/rho1)
u1 = M1*c1
v1 = 0
e1 = p1/(gamma-1)+rho1*(u1*u1 + v1*v1)/2
U = np.zeros((4,nmax,mmax))
F = np.zeros((4,nmax,mmax))
E = np.zeros((4,nmax,mmax))
N = np.zeros((8,nmax,mmax))
S = np.zeros((4,nmax,mmax))
V = np.zeros((nmax,mmax))
newU = np.zeros((4,nmax,mmax))
newU2 = np.zeros((4,nmax,mmax))

rho = np.zeros((nmax,mmax))
u = np.zeros((nmax,mmax))
v = np.zeros((nmax,mmax))
p = np.zeros((nmax,mmax))
e = np.zeros((nmax,mmax))
c = np.zeros((nmax,mmax))


uu1 = np.zeros((nmax,mmax))
uu2 = np.zeros((nmax,mmax))
uu3 = np.zeros((nmax,mmax))
uu4 = np.zeros((nmax,mmax))
uc = np.zeros((nmax,mmax))
pterm = np.zeros((nmax,mmax))
ptermabs = np.zeros((nmax,mmax))
mindt = 0
dphi = np.zeros((nmax,mmax))
dmhi = np.zeros((nmax,mmax))
dphj = np.zeros((nmax,mmax))
dmhj = np.zeros((nmax,mmax))

Edot = np.zeros((4,nmax,mmax))
Epdot = np.zeros((4,nmax,mmax))
Emdot = np.zeros((4,nmax,mmax))
Fdot = np.zeros((4,nmax,mmax))
Fpdot = np.zeros((4,nmax,mmax))
Fmdot = np.zeros((4,nmax,mmax))
def det(a1,a2,a3,b1,b2,b3,c1,c2,c3):
    return a1*(b2*c3 - b3*c2) + a2*(b3*c1 - b1*c3) + a3*(b1*c2 - b2*c1) 


for i in range(0,nmax):
    for j in range(0,mmax):
        ndx = bigx[i,j+1] - bigx[i,j]
        ndy = bigy[i,j+1] - bigy[i,j]
        area = np.sqrt(ndx*ndx + ndy*ndy)
        S[0,i,j] = area
        N[0,i,j] = ndy/area
        N[1,i,j] = -ndx/area
        ndx = bigx[i+1,j+1] - bigx[i,j+1]
        ndy = bigy[i+1,j+1] - bigy[i,j+1]
        area = np.sqrt(ndx*ndx + ndy*ndy)
        S[1,i,j] = area
        N[2,i,j] = ndy/area
        N[3,i,j] = -ndx/area
        ndx = bigx[i+1,j] - bigx[i+1,j+1]
        ndy = bigy[i+1,j] - bigy[i+1,j+1]
        area = np.sqrt(ndx*ndx + ndy*ndy)
        S[2,i,j] = area
        N[4,i,j] = ndy/area
        N[5,i,j] = -ndx/area
        ndx = bigx[i,j] - bigx[i+1,j]
        ndy = bigy[i,j] - bigy[i+1,j]
        area = np.sqrt(ndx*ndx + ndy*ndy)
        S[3,i,j] = area
        N[6,i,j] = ndy/area
        N[7,i,j] = -ndx/area
        S1 = det(bigx[i,j],bigy[i,j],1,bigx[i+1,j+1],bigy[i+1,j+1],1,bigx[i+1,j],bigy[i+1,j],1)
        S2 = det(bigx[i,j],bigy[i,j],1,bigx[i,j+1],bigy[i,j+1],1,bigx[i+1,j+1],bigy[i+1,j+1],1)
        V[i,j] = 0.5*(abs(S1) + abs(S2))
        
U[0,:,:] = rho1
U[1,:,:] = rho1*u1
U[2,:,:] = rho1*v1
U[3,:,:] = e1
E[0,:,:] = rho1*u1
E[1,:,:] = rho1*u1*u1 + p1
E[2,:,:] = rho1*u1*v1
E[3,:,:] = (e1 + p1)*u1
F[0,:,:] = rho1*v1
F[1,:,:] = rho1*u1*v1
F[2,:,:] = rho1*v1*v1 + p1
F[3,:,:] = (e1 + p1)*v1


newU2 = U.copy()
newU = U.copy()

for t in range(0,tmax):
    #边界条件
    U[:,:,mmax - 1] = 2*U[:,:,mmax-2] - U[:,:,mmax-3]
    U[0,nmax-1,1:mmax-1] = U[0,nmax-2,1:mmax-1]
    U[1,nmax-1,1:mmax-1] = U[1,nmax-2,1:mmax-1]
    U[2,nmax-1,1:mmax-1] = -U[2,nmax-2,1:mmax-1]
    U[3,nmax-1,1:mmax-1] = U[3,nmax-2,1:mmax-1]
    theta = 10.94*np.pi/180
    cthe = np.cos(theta)
    sthe = np.sin(theta)
    vn = U[1,nmax-2,IS-1:mmax]*sthe + U[2,nmax-2,IS-1:mmax]*cthe 
    vt = U[1,nmax-2,IS-1:mmax]*cthe - U[2,nmax-2,IS-1:mmax]*sthe 
    U[1,nmax-1,IS-1:mmax] = vt*cthe - vn*sthe
    U[2,nmax-1,IS-1:mmax] = -vt*sthe - vn*cthe
    U[0,nmax-1,IS-1:mmax] = U[0,nmax-2,IS-1:mmax]
    U[3,nmax-1,IS-1:mmax] = U[3,nmax-2,IS-1:mmax]
    U[0,0,:] = U[0,1,:]
    U[1,0,:] = U[1,1,:]
    U[2,0,:] = -U[2,1,:]
    U[3,0,:] = U[3,1,:]
    
    #计算参数
    rho[:,:] = U[0,:,:]
    u[:,:] = U[1,:,:]/rho[:,:]
    v[:,:] = U[2,:,:]/rho[:,:]
    e[:,:] = U[3,:,:]
    p[:,:] = (gamma - 1)*(e[:,:] - rho[:,:]*(u[:,:]*u[:,:] + v[:,:]*v[:,:])/2)
    E[0,:,:] = rho[:,:]*u[:,:]
    E[1,:,:] = rho[:,:]*u[:,:]*u[:,:] + p[:,:]
    E[2,:,:] = rho[:,:]*u[:,:]*v[:,:]
    E[3,:,:] = (e[:,:] + p[:,:])*u[:,:]
    F[0,:,:] = rho[:,:]*v[:,:]
    F[1,:,:] = rho[:,:]*u[:,:]*v[:,:]
    F[2,:,:] = rho[:,:]*v[:,:]*v[:,:] + p[:,:]
    F[3,:,:] = (e[:,:] + p[:,:])*v[:,:]
    
    #更新时间
    mindt = 100
    c[:,:] = (gamma*p[:,:]/rho[:,:])**0.5
    for i in range(1,nmax-1):
        for j in range(1,mmax-1):
            nowdt = 1/(abs(u[i,j])/dx + abs(v[i,j])/dymin + c[i,j]*(1/dx/dx+1/dymin/dymin)**0.5)
            if mindt > nowdt:
                mindt = nowdt
                
    uu1[:,:] = abs(u[:,:]*N[0,:,:] + v[:,:]*N[1,:,:])
    uu2[:,:] = abs(u[:,:]*N[2,:,:] + v[:,:]*N[3,:,:])
    uu3[:,:] = abs(u[:,:]*N[4,:,:] + v[:,:]*N[5,:,:])
    uu4[:,:] = abs(u[:,:]*N[6,:,:] + v[:,:]*N[7,:,:])
    uc[:,:] = uu2[:,:] + c[:,:]
    for i in range(1,nmax-1):
        for j in range(1,mmax-1):
            pterm[i,j] = p[i,j+1] + 2*p[i,j] + p[i,j-1]
            ptermabs[i,j] = abs(p[i,j+1] - 2*p[i,j] + p[i,j-1])
            dphi[i,j] = dis*uc[i,j]*(ptermabs[i,j]/pterm[i,j])
    uc[:,:] = uu4[:,:] + c[:,:]
    for i in range(1,nmax-1):
        for j in range(1,mmax-1):
            pterm[i,j] = p[i,j] + 2*p[i,j-1] + p[i,j-2]
            if j == 1:
                ptermabs[i,j] = 0
            else:
                ptermabs[i,j] = abs(p[i,j] - 2*p[i,j-1] + p[i,j-2])
            dmhi[i,j] = dis*uc[i,j]*(ptermabs[i,j]/pterm[i,j])   
    uc[:,:] = uu3[:,:] + c[:,:]
    for i in range(1,nmax-1):
        for j in range(1,mmax-1):
            pterm[i,j] = p[i+1,j] + 2*p[i,j] + p[i-1,j]
            ptermabs[i,j] = abs(p[i+1,j] - 2*p[i,j] + p[i-1,j])
            dphj[i,j] = dis*uc[i,j]*(ptermabs[i,j]/pterm[i,j])
    uc[:,:] = uu1[:,:] + c[:,:]
    for i in range(1,nmax-1):
        for j in range(1,mmax-1):
            pterm[i,j] = p[i,j] + 2*p[i-1,j] + p[i-2,j]
            if i == 1:
                ptermabs[i,j] = 0
            else:
                ptermabs[i,j] = abs(p[i,j] - 2*p[i-1,j] + p[i-2,j])
            dmhj[i,j] = dis*uc[i,j]*(ptermabs[i,j]/pterm[i,j])       
            
    
    for i in range(1,nmax-1):
        for j in range(1,mmax-1):
            Epdot[:,i,j] = E[:,i,j+1]*N[2,i,j] + F[:,i,j+1]*N[3,i,j]
            Edot[:,i,j] = E[:,i,j]*N[6,i,j] + F[:,i,j]*N[7,i,j]
            Fpdot[:,i,j] = E[:,i+1,j]*N[4,i,j] + F[:,i+1,j]*N[5,i,j]
            Fdot[:,i,j] = E[:,i,j]*N[0,i,j] + F[:,i,j]*N[1,i,j]
            newU[:,i,j] = U[:,i,j] - mindt/V[i,j]*(
                  (Epdot[:,i,j] - dphi[i,j]*(U[:,i,j+1]-U[:,i,j]))*S[1,i,j]
                + (Edot[:,i,j] + dmhi[i,j]*(U[:,i,j] - U[:,i,j-1]))*S[3,i,j] 
                + (Fpdot[:,i,j] - dphj[i,j]*(U[:,i+1,j] - U[:,i,j]))*S[2,i,j]
                + (Fdot[:,i,j] + dmhj[i,j]*(U[:,i,j] - U[:,i-1,j]))*S[0,i,j])
            
    newU[:,:,mmax - 1] = 2*newU[:,:,mmax-2] - newU[:,:,mmax-3]
    newU[0,nmax-1,1:mmax-1] = newU[0,nmax-2,1:mmax-1]
    newU[1,nmax-1,1:mmax-1] = newU[1,nmax-2,1:mmax-1]
    newU[2,nmax-1,1:mmax-1] = -newU[2,nmax-2,1:mmax-1]
    newU[3,nmax-1,1:mmax-1] = newU[3,nmax-2,1:mmax-1]
    theta = 10.94*np.pi/180
    cthe = np.cos(theta)
    sthe = np.sin(theta)
    vn = newU[1,nmax-2,IS-1:mmax]*sthe + newU[2,nmax-2,IS-1:mmax]*cthe 
    vt = newU[1,nmax-2,IS-1:mmax]*cthe - newU[2,nmax-2,IS-1:mmax]*sthe 
    newU[1,nmax-1,IS-1:mmax] = vt*cthe - vn*sthe
    newU[2,nmax-1,IS-1:mmax] = -vt*sthe - vn*cthe
    newU[0,nmax-1,IS-1:mmax] = newU[0,nmax-2,IS-1:mmax]
    newU[3,nmax-1,IS-1:mmax] = newU[3,nmax-2,IS-1:mmax]
    newU[0,0,:] = newU[0,1,:]
    newU[1,0,:] = newU[1,1,:]
    newU[2,0,:] = -newU[2,1,:]
    newU[3,0,:] = newU[3,1,:]
    
    #计算参数
    rho[:,:] = newU[0,:,:]
    u[:,:] = newU[1,:,:]/rho[:,:]
    v[:,:] = newU[2,:,:]/rho[:,:]
    e[:,:] = newU[3,:,:]
    p[:,:] = (gamma - 1)*(e[:,:] - rho[:,:]*(u[:,:]*u[:,:] + v[:,:]*v[:,:])/2)
    E[0,:,:] = rho[:,:]*u[:,:]
    E[1,:,:] = rho[:,:]*u[:,:]*u[:,:] + p[:,:]
    E[2,:,:] = rho[:,:]*u[:,:]*v[:,:]
    E[3,:,:] = (e[:,:] + p[:,:])*u[:,:]
    F[0,:,:] = rho[:,:]*v[:,:]
    F[1,:,:] = rho[:,:]*u[:,:]*v[:,:]
    F[2,:,:] = rho[:,:]*v[:,:]*v[:,:] + p[:,:]
    F[3,:,:] = (e[:,:] + p[:,:])*v[:,:]

    uu1[:,:] = abs(u[:,:]*N[0,:,:] + v[:,:]*N[1,:,:])
    uu2[:,:] = abs(u[:,:]*N[2,:,:] + v[:,:]*N[3,:,:])
    uu3[:,:] = abs(u[:,:]*N[4,:,:] + v[:,:]*N[5,:,:])
    uu4[:,:] = abs(u[:,:]*N[6,:,:] + v[:,:]*N[7,:,:])
    uc[:,:] = uu2[:,:] + c[:,:]
    for i in range(1,nmax-1):
        for j in range(1,mmax-1):
            if j == 1:
                ptermabs[i,j] = 0
                pterm[i,j] = 1
            else:
                pterm[i,j] = p[i,j] + 2*p[i,j-1] + p[i,j-2]
                ptermabs[i,j] = abs(p[i,j] - 2*p[i,j-1] + p[i,j-2])
            dphi[i,j] = dis*uc[i,j]*(ptermabs[i,j]/pterm[i,j])
    uc[:,:] = uu4[:,:] + c[:,:]
    for i in range(1,nmax-1):
        for j in range(1,mmax-1):
            if j < 3:
                ptermabs[i,j] = 0
            else:
                ptermabs[i,j] = abs(p[i,j-1] - 2*p[i,j-2] + p[i,j-3])
                pterm[i,j] = p[i,j-1] + 2*p[i,j-2] + p[i,j-3]
            dmhi[i,j] = dis*uc[i,j]*(ptermabs[i,j]/pterm[i,j])   
    uc[:,:] = uu3[:,:] + c[:,:]
    for i in range(1,nmax-1):
        for j in range(1,mmax-1):
            if i == 1:
                ptermabs[i,j] = 0
            else:
                ptermabs[i,j] = abs(p[i,j] - 2*p[i-1,j] + p[i-2,j])
                pterm[i,j] = p[i,j] + 2*p[i-1,j] + p[i-2,j]
            dphj[i,j] = dis*uc[i,j]*(ptermabs[i,j]/pterm[i,j])
    uc[:,:] = uu1[:,:] + c[:,:]
    for i in range(1,nmax-1):
        for j in range(1,mmax-1):
            if i  < 3:
                ptermabs[i,j] = 0
            else:
                ptermabs[i,j] = abs(p[i-1,j] - 2*p[i-2,j] + p[i-3,j])
                pterm[i,j] = p[i-1,j] + 2*p[i-2,j] + p[i-3,j]
            dmhj[i,j] = dis*uc[i,j]*(ptermabs[i,j]/pterm[i,j])   

    for i in range(1,nmax-1):
        for j in range(1,mmax-1):
            Edot[:,i,j] = E[:,i,j]*N[2,i,j] + F[:,i,j]*N[3,i,j]
            Emdot[:,i,j] = E[:,i,j-1]*N[6,i,j] + F[:,i,j-1]*N[7,i,j]
            Fdot[:,i,j] = E[:,i,j]*N[4,i,j] + F[:,i,j]*N[5,i,j]
            Fmdot[:,i,j] = E[:,i-1,j]*N[0,i,j] + F[:,i-1,j]*N[1,i,j]
                  #          newU(a,b,i) = 0.5*(U(a,b,i) + newU_(a,b,i) - dt./V(a,b).*...
                  # ((E_dot(a,b,i)-dphi.*(U(a,3:mmax,i)-U(a,b,i))).*S(a,b,2)... 
                  # +(Em_dot(a,b,i)+dmhi.*(U(a,b,i)-U(a,1:mmax-2,i))).*S(a,b,4)...
                  # +(F_dot(a,b,i)-dphj.*(U(3:nmax,b,i)-U(a,b,i))).*S(a,b,3)...
                  # +(Fm_dot(a,b,i)+dmhj.*(U(a,b,i)-U(1:nmax-2,b,i))).*S(a,b,1)));
            newU2[:,i,j] = 0.5*(newU[:,i,j] + U[:,i,j] - mindt/V[i,j]*(
                  (Edot[:,i,j] - dphi[i,j]*(U[:,i,j+1]-U[:,i,j]))*S[1,i,j]
                + (Emdot[:,i,j] + dmhi[i,j]*(U[:,i,j] - U[:,i,j-1]))*S[3,i,j] 
                + (Fdot[:,i,j] - dphj[i,j]*(U[:,i+1,j] - U[:,i,j]))*S[2,i,j]
                + (Fmdot[:,i,j] + dmhj[i,j]*(U[:,i,j] - U[:,i-1,j]))*S[0,i,j]))
    U[:,:,:] = newU2[:,:,:]
    tes = 1
    