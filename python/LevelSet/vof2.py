import numpy as np
# 数值模拟方法和运动界面追踪
# 刘儒勋 王志峰

# 积分平均格式

nmax = 8
c0 = np.zeros((nmax,nmax))
c1 = np.zeros((nmax,nmax))
c2 = np.zeros((nmax,nmax))
gam = np.zeros((nmax,nmax))
h = 1
fi = np.zeros((nmax,nmax))
while(1):
    for i in range(nmax):
        for j in range(nmax):
            cita = 0
            if c0[i+1,j] - c0[i,j] > 1e-6:
                cita = (c0[i,j] - c0[i-1,j])/(c0[i+1,j] - c0[i,j])
            gam[i,j] = fi[cita]*(c0[i+1,j] - c0[i,j])/h
    for i in range(1,nmax-1):
        for j in range(1,nmax-1):
            vx = 0
            dt = 0
            sigma = dt / h
            d1 = -sigma*max(0,vx)*(c0[i,j] - c0[i-1,j])
            d2 = sigma/2*h*max(0,vx)*(vx*sigma - 1)*(gam[i,j] - gam[i-1,j])
            d3 = -sigma*min(0,vx)*(c0[i+1,j] - c0[i,j])
            d4 = sigma/2*h*min(0,vx)*(vx*sigma + 1)*(gam[i+1,j] - gam[i,j])
            c1[i,j] = c0[i,j] + d1 + d2 + d3 + d4
            
# Hirt Nichols
f1 = np.zeros((nmax,nmax))
f0 = np.zeros((nmax,nmax))
yi = np.zeros((nmax,nmax))
xi = np.zeros((nmax,nmax))
ihv = np.zeros((nmax,nmax))
u = np.zeros((nmax,nmax))
while(1):
    for i in range(nmax):
        for j in range(nmax):
            f1[i,j] = f0[i,j]
    for i in range(1,nmax-1):
        for j in range(1,nmax-1):
            yi[i,j] = (f0[i,j-1] + f0[i,j] + f0[i,j+1])*h
            xi[i,j] = (f0[i-1,j] + f0[i,j] + f0[i+1,j])*h
    for i in range(2,nmax-2):
        for j in range(2,nmax-2):
            sl1 = abs(2*(yi[i+1,j] - yi[i-1,j])/(4*h))
            sl2 = abs(2*(xi[i,j+1] - xi[i,j-1])/(4*h))
        if sl1 < sl2:
            ihv[i,j] = 0 # 垂直
        else:
            ihv[i,j] = 1 # 水平
    emikro = 1
    # X 方向
    for i in range(2,nmax-2):
        for j in range(2,nmax-2):
            dv = (abs(u[i,j]*dt))
            if u[i,j] > 0:
                id = i
                ia = i + 1
            else:
                id = i + 1
                ia = i
            if ihv[id,j] == 1: # 水平的话，iad 就是即将要进入的那个
                iad = ia
            else:
                iad = id
            if f0[ia,j] < emikro:
                iad = i
            # 目标网格空气体积乘速度，减去自己网格空气体积
            cf = max((1 - f0[iad,j])*dv - (1 - f0[id,j])*h,0)
            # 再加上自己网格水的体积
            df = min(f0[iad,j]*dv+cf,f0[id,j]*h)
            f1[id,j] = f1[id,j] - df / h
            f1[ia,j] = f1[ia,j] + df / h
    for i in range(1,nmax-1):
        for j in range(1,nmax-1):
            if f1[i,j] <  emikro:
                f1[i,j] = 0
                if f1[i-1,j] > 1:
                    f1[i-1,j] = f1[i-1,j] - 1.1*emikro
                if f1[i+1,j] > 1:
                    f1[i+1,j] = f1[i+1,j] - 1.1*emikro
                if f1[i,j-1] > 1:
                    f1[i,j-1] = f1[i,j-1] - 1.1*emikro
                if f1[i,j+1] > 1:
                    f1[i,j+1] = f1[i,j+1] - 1.1*emikro
                f0[i,j] = min(1,max(f1[i,j],0))
                    
# Young`s VOF
f1 = 0
f2 = 0
f3 = 0
f4 = 0
itype = 1
if itype == 1:
    s1 = 0
    s3 = np.sqrt(2*c*cotalfa)
    s4= 0
    s2 = np.sqrt(2*c*tanalfa)
    if u1 > 0:
        if u1 * dt < (1 - s2)*dy:
            f1 = 0
        else:
            temp1 = u1*dt - (1 - s2)*dy
    