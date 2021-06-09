import numpy as np
# D:\FluidSim\FluidSim\Immersed\IB2d-master\matIB2d\Examples\Example_VIV_Cylinder\Single_VIV_Cylinder_128
import scipy.io as scio
dx = 0.0078125
nmax = 128
x = np.zeros((nmax))
y = np.zeros((nmax))
u = np.zeros((nmax,nmax))
v = np.zeros((nmax,nmax))
for i in range(0,nmax):
    x[i] = y[i] = i * dx
    
 
dataFile = 'xlag.mat'
xlagx = scio.loadmat(dataFile)['xLag']
xlag = xlagx[:,0]
dataFile = 'ylag.mat'
ylagy = scio.loadmat(dataFile)['yLag']
ylag = ylagy[:,0]
Lx = 1
supp = 4
for t in range(0,1):
    # please_Move_Lagrangian_Point_Positions
    # xLag_Next = xLag_Prev + dt* int( u(x,t) delta( x - xLag_n ) dX ) 
    Xindex = np.zeros((570,16),dtype = int)
    Yindex = np.zeros((570,16),dtype = int)
    Xdist = np.zeros((570,16))
    Ydist = np.zeros((570,16))
    Xdelta = np.zeros((570,16))
    Ydelta = np.zeros((570,16))
    Xmove = np.zeros((570))
    Ymove = np.zeros((570))
    Xnext = np.zeros((570))
    Ynext = np.zeros((570))
    for i in range(0,nmax):
        idx = int(xlag[i] / dx)
        Xindex[i,0] = Xindex[i,4] = Xindex[i,8] = Xindex[i,12] = idx - 1
        Xindex[i,1] = Xindex[i,5] = Xindex[i,9] = Xindex[i,13] = idx
        Xindex[i,2] = Xindex[i,6] = Xindex[i,10] = Xindex[i,14] = idx + 1
        Xindex[i,3] = Xindex[i,7] = Xindex[i,11] = Xindex[i,15] = idx + 2
        idx = int(ylag[i] / dx)
        Yindex[i,0] = Yindex[i,4] = Yindex[i,8] = Yindex[i,12] = idx - 1
        Yindex[i,1] = Yindex[i,5] = Yindex[i,9] = Yindex[i,13] = idx
        Yindex[i,2] = Yindex[i,6] = Yindex[i,10] = Yindex[i,14] =idx + 1
        Yindex[i,3] = Yindex[i,7] = Yindex[i,11] = Yindex[i,15] =idx + 2
        for j in range(0,16):
            ab = abs(x[Xindex[i,j]] - np.mod(xlag[i],Lx))
            Xdist[i,j] = min(ab,Lx - ab)
            ab = abs(y[Yindex[i,j]] - np.mod(ylag[i],1))
            Ydist[i,j] = min(ab,Lx - ab)
            
            rmat = abs(Xdist[i,j]) / dx
            if rmat < 1:
                Xdelta[i,j] = ((3 - 2*rmat + np.sqrt(1 + 4*rmat - 4*rmat*rmat))/(8*dx))
            elif (rmat >= 1)&(rmat < 2):
                Xdelta[i,j] = ((5 - 2*rmat - np.sqrt(-7 + 12*rmat - 4*rmat*rmat))/(8*dx))
            else:
                Xdelta[i,j] = 0
                
            rmat = abs(Ydist[i,j]) / dx
            if rmat < 1:
                Ydelta[i,j] = ((3 - 2*rmat + np.sqrt(1 + 4*rmat - 4*rmat*rmat))/(8*dx))
            elif (rmat >= 1)&(rmat < 2):
                Ydelta[i,j] = ((5 - 2*rmat - np.sqrt(-7 + 12*rmat - 4*rmat*rmat))/(8*dx))
            else:
                Ydelta[i,j] = 0        
            # give_Me_Perturbed_Distance
            Xmove[i] += u[Xindex[i,j],Yindex[i,j]]*Xdelta[i,j]*Ydelta[i,j]*dx*dx
            Ymove[i] += v[Xindex[i,j],Yindex[i,j]]*Xdelta[i,j]*Ydelta[i,j]*dx*dx
    
    dt = 0.0000125
    Xnext[:] = xlag[:] + dt * Xmove[:]
    Ynext[:] = ylag[:] + dt * Ymove[:]
    # please_Find_Lagrangian_Forces_On_Eulerian_grid
    
    # give_Me_Target_Lagrangian_Force_Densities
    target = np.zeros((462,3))
    target[:,0] = Xnext[0:462] # Original x-Values of x-Target Pts.
    target[:,1] = Ynext[0:462] # Original y-Values of y-Target Pts.
    target[:,2] = 25000000 # Stores Target Stiffnesses 
    
    # give_Me_Beam_Lagrangian_Force_Densities
    nbeams = 108
    beam = np.zeros((nbeams,5))
    for i in range(0,nbeams):
        beam[i,0] = i + 462 - 1 # Initialize storage for 1ST NODE for BEAM
        beam[i,1] = i + 462 # Initialize storage for 2ND NODE for BEAM
        beam[i,2] = i + 462 + 1 # Initialize storage for 3RD NODE for BEAM
        beam[i,3] = 5e9 # Stores spring stiffness associated with each spring
        beam[i,4] = (-2.233)/(1e7) # Stores spring resting length associated with each spring
    beam[0,4] = beam[107,4] = (-3.21)/(1e8)
    beam[0,0] = 569
    beam[107,2] = 462
    fx = np.zeros((570,3))
    fy = np.zeros((570,3))
    for i in range(0,nbeams):
        
    
        xpos1 = Xnext[int(beam[i,0])]
        xpos2 = Xnext[int(beam[i,1])]
        xpos3 = Xnext[int(beam[i,2])]
        ypos1 = Ynext[int(beam[i,0])]
        ypos2 = Ynext[int(beam[i,1])]
        ypos3 = Ynext[int(beam[i,2])]
        kbeam = beam[i,3]
        cbeam = beam[i,4]
        ds = 0.00390625
        # check_If_Beam_Points_Pass_Through_Boundary
        dx12 = xpos1 - xpos2
        dx23 = xpos2 - xpos3
        xpn = xpos1
        xqn = xpos2
        xrn  =xpos3
        if abs(dx12) > 5*ds:
            if dx12 < 0:
                xpn = Lx + xpos1
            else:
                xpn = - Lx + xpos1
        if abs(dx23) > 5*ds:
            if dx23 < 0:
                xrn = -Lx + xpos3
            else:
                xrn = Lx + xpos3
                
        dy12 = ypos1 - ypos2
        dy23 = ypos2 - ypos3
        ypn = ypos1
        yqn = ypos2
        yrn = ypos3
        if abs(dy12) > 5*ds:
            if dy12 < 0:
                ypn = Lx + ypos1
            else:
                ypn = - Lx + ypos1
        if abs(dy23) > 5*ds:
            if dx23 < 0:
                yrn = -Lx + ypos3
            else:
                yrn = Lx + ypos3
                
        cross_prod = (xrn - xqn)*(yqn - ypn) - (yrn - yqn)*(xqn - xpn)
        # Force for left node
        bF_x_L = - kbeam * (cross_prod - cbeam) * (yrn - yqn)
        bF_y_L = kbeam * (cross_prod - cbeam) * (xrn - xqn)
        # Force for middle node
        bF_x_M = kbeam * (cross_prod - cbeam) * (yqn - ypn + yrn - yqn)
        bF_y_M = - kbeam * (cross_prod - cbeam) * (xrn - xqn + xqn - xpn)
        # Force for right node
        bF_x_R = - kbeam * (cross_prod - cbeam) * (yqn - ypn)
        bF_y_R = kbeam * (cross_prod - cbeam) * (xqn - xpn)
        fx[i,0] -= bF_x_L
        fy[i,0] -= bF_y_L
        fx[i,1] -= bF_x_M
        fy[i,1] -= bF_y_M
        fx[i,2] -= bF_x_R
        fy[i,2] -= bF_y_R
        
    mu = 0.01
    dt = 2.5/(1e5)
    rho = 1
    # Ahat 
    Ahat = np.zeros((nmax,nmax))
    for i in range(nmax):
        for j in range(nmax):
            Ahat[i,j] = 1 + 2*mu*dt/rho*((np.sin(np.pi * i / nmax)/dx)**2 + (np.sin(np.pi * j / nmax)/dx)**2)
            
    # rhs_u = give_RHS_HALF_Step_Velocity(dt,rho,Nx,Ny,U,Ux,Uy,U_sq_x,UV_y,V,Fx,'x');
    # rhs_v = give_RHS_HALF_Step_Velocity(dt,rho,Nx,Ny,V,Vx,Vy,V_sq_y,UV_x,U,Fy,'y');
    # give_RHS_HALF_Step_Velocity(dt,rho,Nx,Ny,A,Ax,Ay,A_sq_j,AB_j,B,Fj,string)
    # rhs = A + (0.5*dt/rho)*( Fj - 0.5*rho*(A.*Ax + B.*Ay) - 0.5*rho*( A_sq_j + AB_j ) ); %RHS: u-component
    # rhs = A + (0.5*dt/rho)*( Fj - 0.5*rho*( B.*Ax + A.*Ay ) - .5*rho*( AB_j + A_sq_j ) ); %RHS: v-compoent
    
    #     num = -( 1i/dx*sin(2*pi*idX/Nx).*rhs_u_hat + 1i/dy*sin(2*pi*idY/Ny).*rhs_v_hat );
    # den = ( dt/rho*( ( sin(2*pi*idX/Nx)/dx ).^2 + ( sin(2*pi*idY/Ny)/dy ).^2 ) );
    # p_hat = num./den;

    # % Zero out modes.
    # p_hat(1,1) = 0;
    # p_hat(1,Nx/2+1) = 0;
    # p_hat(Ny/2+1,Nx/2+1) = 0;  
    # p_hat(Ny/2+1,availab1) = 0;
  
    # vel_hat = ( rhs_VEL_hat - 1i*dt/(rho*dj)*sin(2*pi*idMat/Nx).*p_hat ) ./ A_hat;
   
    # give_RHS_FULL_Step_Velocity(dt,mu,rho,Nx,Ny,A,Ax,Ay,A_sq_j,AB_j,B,Fj,Axx,Ayy,string)