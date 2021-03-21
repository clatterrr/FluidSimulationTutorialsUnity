import numpy as np
'''
解的偏微分方程为 dU/dt + adU/dx= 0，a为常数
矩阵参数使用数值二重积分计算


'''
#可修改的初始变量
num = 4#num行num列
a = 0.5#上面公式中的a，也就是线性对流的速度
nmax = num*num
initval = np.zeros((num))
for i in range(0,num):
    initval[i] = i

def phi(n,t,x):
    if n == 0:
        return (1-t)*(1-x)
    elif n == 1:
        return (1-t)*x
    elif n == 2:
        return t*(1-x)
    elif n == 3:
        return t*x

def simpon(point0,point1,point2,h):
    return h/3*(point0 + 4*point1 + point2)

#以下变量用于计算数值积分
nodes = 4
p = np.zeros((9,nodes,nodes))
pt = np.zeros((9,nodes,nodes))
ppt = np.zeros((3))
tempt = np.zeros((3))
px = np.zeros((9,nodes,nodes))
ppx = np.zeros((3))
tempx = np.zeros((3))
k1mat = np.zeros((nodes,nodes))
k2mat = np.zeros((nodes,nodes))

#以下变量用于计算逆矩阵
kmat = np.zeros((nmax,nmax))
kinv = np.zeros((nmax,nmax))
bmat = np.zeros((nmax))
temprow = np.zeros((nmax))
res = np.zeros((num,num))

#计算权函数的积分
for i in range(0,nodes):
    for j in range(0,nodes):
        '''
        p6 p7 p8
        p3 p4 p5
        p0 p1 p2
        矩阵有16个元素
        每个元素都需要做三阶精度的数值二重积分
        因此每个元素都需要9个值来算出
        这9个值每个值都由两个数相乘，结果由ppt表示
        一个数是phi，由p表示
        另一个是phi对t求导，由pt表示
        '''
        tempt[:] = tempx[:] = 0
        for z1 in range(0,3):#j y t
            ppt[:] = ppx[:] = 0
            for z2 in range(0,3):# i x 
                p[z1*3+z2,i,j] = phi(i,z1/2,z2/2)
                pt[z1*3+z2,i,j] = phi(j,2,z2/2) - phi(j,1,z2/2)
                ppt[z2] += p[z1*3 + z2,i,j]*pt[z1*3 + z2,i,j] 
                px[z1*3+z2,i,j] = phi(j,z1/2,2) - phi(j,z1/2,1)
                ppx[z2] += p[z1*3 + z2,i,j]*px[z1*3 + z2,i,j] 
            tempt[z1] = simpon(ppt[0],ppt[1],ppt[2],0.5)
            tempx[z1] = simpon(ppx[0],ppx[1],ppx[2],0.5)
        k1mat[i,j] = simpon(tempt[0],tempt[1],tempt[2],0.5)#phi * phit
        k2mat[i,j] = simpon(tempx[0],tempx[1],tempx[2],0.5)#phi * phix

#计算刚度矩阵
idx = np.zeros((4))
for i in range(0,num-1):
    for j in range(0,num-1):
        idx[0] = j*num + i
        idx[1] = idx[0] + 1
        idx[2] = idx[0] + num
        idx[3] = idx[2] + 1
        '''
        idx2 --- idx3
        |       |
        idx0 --- idx1
        
        '''
        for z1 in range(0,4):
            for z2 in range(0,4):
                kmat[(int)(idx[z1]),(int)(idx[z2])] +=(k1mat[z1,z2] + a*k2mat[z1,z2])

for i in range(0,nmax):
    kinv[i,i] = 1
for i in range(0,num):#只给t = 0的那一排赋值
    kmat[i,:] = 0
    kmat[i,i] = 1
    bmat[i] = initval[i]

#高斯消元法求逆矩阵
for i in range(0,nmax):
    if abs(kmat[i,i]) < 1e-10:
        for j in range(i+1,nmax):
            if kmat[j,i] != 0:
                temprow[:] = kmat[i,:]
                kmat[i,:] = kmat[j,:]
                kmat[j,:] = temprow[:]
                temprow[:] = kinv[i,:]
                kinv[i,:] = kinv[j,:]
                kinv[j,:] = temprow[:]
                break
    temp = kmat[i,i]
    for j in range(0,nmax):
        kmat[i,j] = kmat[i,j]/temp
        kinv[i,j] = kinv[i,j]/temp
    for k in range(0,nmax):
        if i == k:
            continue
        if kmat[k,i] == 0:
            continue
        mul = kmat[k,i]
        for j in range(0,nmax):
            kmat[k,j] -= mul*kmat[i,j]
            kinv[k,j] -= mul*kinv[i,j]
# 求解最后结果
for i in range(0,nmax):
    temp = 0
    for j in range(0,nmax):
        temp += kinv[i,j]*bmat[j]
    res[(int)(i/num),(int)(i%num)] = temp
