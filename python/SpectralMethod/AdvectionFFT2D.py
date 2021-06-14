import numpy as np
import matplotlib.pyplot as plt

def FFT(arr,flag):
    dftres = np.asarray(arr,dtype = complex)
    if flag == 1:
        dftres /= total
    for level0 in range(0,level):
        maxi = 2**level0
        xdis = nmax//maxi
        for i in range(0,maxi):
            x0 = int(i % nmax)
            y0 = int(i / nmax)
            even =  dftres[x0,y0] # 蝴蝶变换中，偶数项已经计算好，因此之间复制就行
            odd = 0
            # 蝴蝶变换
            for j in range(0,nmax,xdis):
                idx = j + xdis//2
                x1 = int(idx % nmax)
                y1 = int(idx / nmax)
                odd += arr[x1,y1] * np.exp(flag*2j*np.pi*x0*x1/nmax)#奇数项需要重新计算
            
            # 消去引理
            if flag == 1:
                odd /= total
            dftres[x0,y0] = even + odd
            dftres[x0+maxi,y0] = even - odd

    for level0 in range(0,level):
        maxi = 2**level0
        xdis = total//maxi # 方块之间一下跳跃的距离
        maxii = maxi * nmax # 目前已有的行数
        for i in range(0,maxii):
            x0 = int(i % nmax)
            y0 = int(i / nmax)
            even =  dftres[x0,y0] # 蝴蝶变换中，偶数项已经计算好，因此之间复制就行
            odd = 0
            # 蝴蝶变换
            for j0 in range(0,total,xdis):
                for j1 in range(0,nmax):
                    idx = j0 + xdis//2 + j1
                    x1 = int(idx % nmax)
                    y1 = int(idx / nmax)
                    odd += arr[x1,y1] * np.exp(flag*2j*np.pi*(x0*x1 + y0*y1)/nmax)#奇数项需要重新计算
            # 消去引理
            if flag == 1:
                odd /= total
            dftres[x0,y0] = even + odd
            dftres[x0,y0+maxi] = even - odd
    return dftres

def derivation(q):
    # qhat = FFT(q,-1)
    # qx = np.fft.ifft2(1j*kwave*qhat,1).real
    # qxx = np.fft.ifft2(-kwave**2*qhat,1).real
    # qxxx = np.fft.ifft2(-1j*kwave**3*qhat,1).real
    # qy = np.fft.ifft2(1j*lwave*qhat,1).real
    # qyy = np.fft.ifft2(-lwave**2*qhat,1).real
    # qyyy = np.fft.ifft2(-1j*lwave**3*qhat,1).real
    
    qhat = np.fft.fft2(q)
    qx = np.fft.ifft2(1j*kwave*qhat).real
    qxx = np.fft.ifft2(-kwave**2*qhat).real
    qxxx = np.fft.ifft2(-1j*kwave**3*qhat).real
    qy = np.fft.ifft2(1j*lwave*qhat).real
    qyy = np.fft.ifft2(-lwave**2*qhat).real
    qyyy = np.fft.ifft2(-1j*lwave**3*qhat).real
    a = 1 # 系数
    res = a*(qx + qy)# 对流方程
    # res = - a*(qxx + qyy) # 扩散方程
    return res

L = 4
level = 5
nmax = 2**level
total = nmax*nmax
tmax = 100
amp = 2
b = np.sqrt(amp/2)  # 孤子宽度的倒数
h = np.zeros((tmax,nmax,nmax))
rk = np.zeros((4,nmax,nmax))
dt = 0.01 # 步长不能太大
kwave = np.zeros((nmax,nmax))
lwave = np.zeros((nmax,nmax))
for i in range(0,nmax):
    for j in range(0,nmax):
        h[0,i,j] = amp*1/np.cosh(b*(i+j)-4)**2
        kwave[i,j] = 2*np.pi * i / L
        if i >= nmax//2:
            kwave[i,j] = 2*np.pi * (i - nmax) / L
        lwave[i,j] = 2*np.pi * j / L
        if i >= nmax//2:
            lwave[i,j] = 2*np.pi * (j - nmax) / L
            
x1 = np.linspace(0, 5, nmax)
x2 = np.linspace(0, 5, nmax)

X1, X2 = np.meshgrid(x1, x2)            

for t in range(0,tmax-1):
    # # Runge Kutta 4 阶，第一种写法最难记也不好理解
    # rk[0,:] = -dt*derivation(h[t,:])
    # rk[1,:] = -dt*derivation(h[t,:] + 0.5*rk[0,:])
    # rk[2,:] = -dt*derivation(h[t,:] + 0.5*rk[1,:])
    # rk[3,:] = -dt*derivation(h[t,:] + rk[2,:])
    # h[t+1,:] = h[t,:] + (rk[0,:] + 2*rk[1,:] + 2*rk[2,:] + rk[3,:])/6 
    
    # # 第二种写法就是泰勒公式
    # rk[0,:] = -dt*derivation(h[t,:])
    # rk[1,:] = -dt*derivation(rk[0,:])
    # rk[2,:] = -dt*derivation(rk[1,:])
    # rk[3,:] = -dt*derivation(rk[2,:])
    # h[t+1,:] = h[t,:] + rk[0,:] + rk[1,:]/2 + rk[2,:]/6 + rk[3,:]/24
    
    # 第三种更加简洁明了一些
    rk[0,:] = h[t,:] - dt*derivation(h[t,:])
    rk[1,:] = rk[0,:] - dt*derivation(rk[0,:])
    rk[2,:] = rk[1,:] - dt*derivation(rk[1,:])
    rk[3,:] = rk[2,:] -dt*derivation(rk[2,:])
    h[t+1,:] = h[t,:]*9/24 + rk[0,:]/3 + rk[1,:]/4 + rk[3,:]/24
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X1, X2, h[t+1,:,:], cmap='bwr', linewidth=0)
    fig.colorbar(surf)
    # fig.show()