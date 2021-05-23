import numpy as np
# 参考： http://www.netlib.org/templates/matlab/
A = np.array([[1,4,7],[2,9,7],[5,8,3]],dtype = float)
M = np.array([[1,0,0],[0,1,0],[0,0,1]],dtype = float)
Minv = np.linalg.inv(M)
b = np.array([30,41,30],dtype = float)
nmax = 3
x = np.zeros((nmax))

res = np.dot(Minv,b - np.dot(A,x)) # 残差
tmax = 100
m = 10
V = np.zeros((nmax,m))
H = np.zeros((m+1,m))

y = np.zeros((m))
cs = np.zeros((m))
sn = np.zeros((m))
e1 = np.zeros((nmax))
s = np.zeros((nmax+1))
e1[0] = 1
bnrm2 = np.linalg.norm(b)
rho1 = 1
for t in range(0,tmax):
    res = np.dot(Minv,b - np.dot(A,x)) # 残差
    normr = np.linalg.norm(res)
    s[0:nmax] = normr * e1
    
    ''' Arnoldi 迭代开始 '''
    # 用于产生krylov空间{Ab}的一组正交基{v}
    V[:,0] = res / normr # 第一个正交基，其实就是归一化的b
    for i in range(0,m):
        w = np.dot(A,V[:,i])
        for k in range(0,i+1):
            H[k,i] = np.dot(np.transpose(w),V[:,k])
            w = w - H[k,i]*V[:,k]
        H[i+1,i] = np.linalg.norm(w)
        V[:,i+1] = w / H[i+1,i] # 其余的w的正交基
        ''' Arnoldi 迭代结束 '''
        for k in range(0,i):
            temp = cs[k]*H[k,i] + sn[k]*H[k+1,i]
            H[k+1,i] = -sn[k]*H[k,i] + cs[k]*H[k+1,i]
            H[k,i] = temp
        
        paraA = H[i,i]
        paraB = H[i+1,i]
        if paraB == 0:
            cs[i] = 1
            sn[i] = 0
        elif abs(paraB) > abs(paraA):
            temp = paraA / paraB
            sn[i] = 1 / np.sqrt(1 + temp*temp)
            cs[i] = temp * sn[i]
        else:
            temp = paraB / paraA
            cs[i] = 1 / np.sqrt(1 + temp*temp)
            sn[i] = temp * cs[i]
        temp = cs[i]*s[i]
        s[i+1] = -sn[i]*s[i]
        s[i] = temp
        H[i,i] = cs[i]*H[i,i] + sn[i]*H[i+1,i]
        H[i+1,i] = 0
        error = abs(s[i+1]) / bnrm2
        if error < 1e-8:
            # 解方程 argmin || beta * e1 - Hmym ||2
            temph = H[0:i+1,0:i+1]
            hinv = np.linalg.inv(temph)
            y[0:i+1] = np.dot(hinv,s[0:i+1])
            x = x + np.dot(V[:,0:i+1],y[0:i+1])
            break
    if error < 1e-8:
        break
    # y = H[1:m,1:m] / s[1:m]
    # x = x + V[:,1:m]*y
    # r = b - np.dot(A,x)
    # normr = np.linalg.norm(r)
    # s[i+1] = normr
        
            
            
    
    