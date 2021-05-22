import numpy as np

def FFT(arr):
    dftres = np.asarray(arr,dtype = complex)
    for level0 in range(0,level):
        maxi = 2**level0
        interval = nmax//maxi
        for i in range(0,maxi): 
            even =  dftres[i] 
            odd = 0
            for j in range(0,nmax,interval):
                idx = j + interval//2
                odd += arr[idx] * np.exp(-2j*np.pi*idx*i/nmax)
            dftres[i] = even + odd
            dftres[i + maxi] = even - odd
    return dftres.copy()

def IFFT(arr):
    idftres = np.asarray(arr,dtype = complex)/nmax
    for level0 in range(0,level):
        maxi = 2**level0
        interval = nmax//maxi
        for i in range(0,maxi): 
            even =  idftres[i] 
            odd = 0
            for j in range(0,nmax,interval):
                idx = j + interval//2
                odd += arr[idx] * np.exp(2j*np.pi*idx*i/nmax)/nmax
            idftres[i] = even + odd
            idftres[i + maxi] = even - odd
    return idftres.copy()


def derivation(q):
    qhat = FFT(q)
    qx = IFFT(1j*k[:]*qhat[:]).real
    qxx = IFFT(-k[:]**2*qhat[:]).real
    qxxx = IFFT(-1j*k[:]**3*qhat[:]).real
    
    # qhat = np.fft.fft(q)
    # qx = np.fft.ifft(1j*k[:]*qhat[:]).real
    # qxx = np.fft.ifft(-k[:]**2*qhat[:]).real
    # qxxx = np.fft.ifft(-1j*k[:]**3*qhat[:]).real
    
    a = 1 # 系数
    # res = a*qx # 对流方程
    # res = - a*qxx # 扩散方程
    res = q*qx # 无粘性burgers 方程
    res = q*qx + qxxx
    # res = q*qx - a*qxx # 粘性burgers 方程
    return res

L = 16
level = 5
nmax = 2**level
tmax = 200
amp = 2
b = np.sqrt(amp/2)  # 孤子宽度的倒数
h = np.zeros((tmax,nmax))
M = np.zeros((nmax))
rk = np.zeros((4,nmax))
dt = 0.01 # 步长不能太大
for i in range(0,nmax):
    h[0,i] = amp*1/np.cosh(b*i-4)**2
    M[i] = i
    if (i >= nmax/2):
        M[i] -= nmax
k = 2*np.pi*M/L#波数
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