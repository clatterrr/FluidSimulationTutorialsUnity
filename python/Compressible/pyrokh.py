import numpy as np

nmax = 72

dens = np.zeros((nmax,nmax))
xmom = np.zeros((nmax,nmax))
ymom = np.zeros((nmax,nmax))
ener = np.zeros((nmax,nmax))
cs = np.zeros((nmax,nmax))
Q = np.zeros((nmax,nmax,4))

dens[:, :] = 1.0
xmom[:, :] = 0.0
ymom[:, :] = 0.0

rho_1 = 1
v_1 = -0.5
rho_2 = 2
v_2 = 0.5

gamma = 1.4

dy = 0.025
w0 = 0.01
vm = 0.5*(v_1 - v_2)
rhom = 0.5*(rho_1 - rho_2)

mygx2d = np.zeros((72,72))
mygy2d = np.zeros((72,72))
griddx = 0.015625
for i in range(0,72):
    for j in range(0,72):
        mygx2d[i,j] = -0.0546875 + i * griddx
        mygy2d[i,j] = -0.0546875 + j * griddx
        
        if mygy2d[i,j] < 0.25:
            dens[i,j] = rho_1 - rhom*np.exp((mygy2d[i,j] - 0.25)/dy)
            xmom[i,j] = v_1 - vm*np.exp((mygy2d[i,j] - 0.25)/dy)
        elif mygy2d[i,j] < 0.5:
            dens[i,j]  = rho_2 + rhom*np.exp((0.25 - mygy2d[i,j])/dy)
            xmom[i,j] = v_2 + vm*np.exp((0.25 - mygy2d[i,j])/dy)
        elif mygy2d[i,j] < 0.75:
            dens[i,j] = rho_2 + rhom*np.exp((mygy2d[i,j] - 0.75)/dy)
            xmom[i,j] = v_2 + vm*np.exp((mygy2d[i,j] - 0.75)/dy)
        else:
            dens[i,j] = rho_1 - rhom*np.exp((0.75 - mygy2d[i,j])/dy)
            xmom[i,j] = v_1 - vm*np.exp((0.75 - mygy2d[i,j])/dy)
            
        xmom[i,j] *= dens[i,j]
        ymom[i,j] = dens[i,j] * w0 * np.sin(4*np.pi*mygx2d[i,j])
        p = 2.5
        ener[i,j] = p/(gamma - 1.0) + 0.5*(xmom[i,j]**2 + ymom[i,j]**2)/dens[i,j]
        cs[i,j] = np.sqrt(gamma*p/dens[i,j])
        
e = np.zeros((72,72))
U = np.zeros((72,72,4))
U[:,:,0] = dens[:,:]
U[:,:,1] = ener[:,:]
U[:,:,2] = xmom[:,:]
U[:,:,3] = ymom[:,:]
for t in range(0,1):
    cfl = 0.8
    xtmp = griddx/(abs(xmom) + cs)
    ytmp = griddx/(abs(xmom) + cs)
    dt = cfl*float(min(xtmp.min(), ytmp.min()))
    
    
    
    Q[:,:,0] = U[:,:,0]
    Q[:,:,1] = xmom[:,:] / dens[:,:]
    Q[:,:,2] = ymom[:,:] / dens[:,:]
    e[:,:] = (ener[:,:] - 0.5*Q[:,:,0]*(Q[:,:,1]**2 + Q[:,:,2]**2))/dens[:,:]
    Q[:,:,3] = dens[:,:]*e[:,:]*(gamma - 1)
    
    z0 = 0.75
    z1 = 0.85
    z = np.zeros((nmax,nmax))
    xi = np.zeros((nmax,nmax))
    for i in  range(0,nmax):
        for j in range(0,nmax):
            xi[i,j] = np.minimum(1.0, np.maximum(0.0, 1.0 - (z[i,j] - z0)/(z1 - z0)))
            # if t1 > 0 & t2 > delta :
            # x[i,j] = 1
            