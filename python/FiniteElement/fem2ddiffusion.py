import numpy as np
'''
解的偏微分方程，扩散 dU/dt = ad^2U/dx^2，a为常数
矩阵参数使用数值二重积分计算
'''
#可修改的初始变量
nnum = 10#变量个数，网格个数，至少是3
tnum = 100#时间步数
a = 1#上面公式中的a，也就是扩散的速度
nmax = tnum*nnum
initval = np.zeros((nnum))#t = 0的初始条件
for i in range(0,nnum):
    #initval[i] = -(i/2-3)**2 + 10
    initval[i] = np.sin(i*3/nnum)

#接下来就是硅基生物的事了
#权函数
def phit(n,t):
    if ((int)(n / 3) ==  0):
        return 1 - t
    elif ((int)(n / 3) == 1):
        return t
def phix(n,x):
    if ((n % 3) ==  0):
        return x*x/2 - x/2
    elif ((n % 3) ==  1):
        return 1 - x*x
    elif ((n % 3) ==  2):
        return x*x/2 + x/2

def bool5(point0,point1,point2,point3,point4,h):#五阶精度
    return 2*h/45*(7*point0 + 32*point1 + 
                   12*point2 + 32*point3 + 7*point4)

def simpon(point0,point1,point2,h):#三阶精度
    return h/3*(point0 + 4*point1 + point2)

#以下变量用于计算装配矩阵
nodes = 6#权函数数量，t有两种，x轴三种，共6种
IntegralNum = 5#五阶精度的数值积分需无个点
p = np.zeros((IntegralNum**2,nodes,nodes))
pt = np.zeros((IntegralNum**2,nodes,nodes))
px = np.zeros((IntegralNum**2,nodes,nodes))
k1mat = np.zeros((nodes,nodes))#phi * phit^T
k2mat = np.zeros((nodes,nodes))#用分部积分替换后的左半，phi * phix^T这半直接计算无需数值积分
k3mat = np.zeros((nodes,nodes))#用分部积分替换后右半，phix * phix^T 这半需要计算数值积分
k4mat = np.zeros((nodes,nodes))

#计算数值积分
ppt = np.zeros((IntegralNum))
tempt = np.zeros((IntegralNum))
ppx = np.zeros((IntegralNum))
tempx = np.zeros((IntegralNum))
pxpx = np.zeros((IntegralNum))
tempxpx = np.zeros((IntegralNum))
#以下变量用于计算逆矩阵
kmat = np.zeros((nmax,nmax))
bmat = np.zeros((nmax))
res = np.zeros((tnum,nnum))

#计算权函数的积分，当然要是你闲得慌也可以手算这3x36 = 108个元素
# 数值积分的结果与初始条件没啥关系，所以我更倾向于用一个程序算数值积分
# 然后算好的结果填到另一个程序里去，就不用每次重新算数值积分了
for i in range(0,nodes):#说明K是一个6x6的矩阵
    for j in range(0,nodes):
        '''
        矩阵有36个元素
        每个元素都需要做五阶精度的数值二重积分
        因此每个元素都需要25个值来算出
        虽然t轴只需要二阶精度，x轴只需要四阶精度
        但这里还是用了五阶精度，对称一点
        '''
        #下面计算矩阵36个元素之中的一个
        for z1 in range(0,IntegralNum):#i y t 
            # 积分区域0到，数值积分取五个点的值，即0,0.25,0.5,0.75,1，对应z1/4
            ppt[:] = ppx[:] = 0
            for z2 in range(0,IntegralNum):# j x 
            # 积分区域-1到1，数值积分三个点分别为-1,-0.5,0,0.5,1，对应z2/2 - 1
                index = z1 * IntegralNum + z2
                p[index,i,j] = phit(i,z1/4) * phix(i,z2/2-1)
                pt[index,i,j] = (phit(j,1) - phit(j,0))*phix(j,z2/2-1)
                ppt[z2] = p[index,i,j]*pt[index,i,j] 
                #接下来两处均用二阶中心差分计算斜率
                pxpx[z2] = phit(i,z1/4)*(phix(i,z2/2) - phix(i,z2/2-2)
                    )/2*phit(j,z1/4)*(phix(j,z2/2) - phix(j,z2/2-2))/2 
                
            tempt[z1] = bool5(ppt[0],ppt[1],ppt[2],ppt[3],ppt[4],0.5)
            tempx[z1] = phit(i,z1/4)*phit(j,z1/4)*(phix(i,1)*(phix(j,2) - phix(j,0))/2 
                        - phix(i,-1)*(phix(j,0) - phix(j,-2))/2)#分部积分后无需积x轴的部分
            tempxpx[z1] = bool5(pxpx[0],pxpx[1],pxpx[2],pxpx[3],pxpx[4],0.5)
            
        k1mat[i,j] = bool5(tempt[0],tempt[1],tempt[2],tempt[3],tempt[4],0.25)#phi * phidt
        k2mat[i,j] = bool5(tempx[0],tempx[1],tempx[2],tempx[3],tempx[4],0.25)#phi * phidt
        k3mat[i,j] = bool5(tempxpx[0],tempxpx[1],tempxpx[2],tempxpx[3],tempxpx[4],0.25)#phidx * phidx
        k4mat[i,j] = k1mat[i,j] - a*(k2mat[i,j] - k3mat[i,j])

#计算并装配刚度矩阵
idx = np.zeros((nodes))
for i in range(0,tnum-1):
    for j in range(0,nnum-2):
        idx[0] = i*nnum + j
        idx[1] = idx[0] + 1
        idx[2] = idx[0] + 2
        
        idx[3] = idx[0] + nnum
        idx[4] = idx[3] + 1 
        idx[5] = idx[3] + 2
        '''
        idx2 --- idx5
          \        \
        idx1 --- idx4 
        |           /        
        idx0 --- idx3
        
        '''
        for z1 in range(0,nodes):
            for z2 in range(0,nodes):
                kmat[(int)(idx[z1]),(int)(idx[z2])] += k4mat[z1,z2]

for i in range(0,nnum):#只给t = 0的那一排赋
    kmat[i,:] = 0
    kmat[i,i] = 1
    bmat[i] = initval[i]
    
# 下面强制第0个网格和最后网格的值永远为零
# 虽然按照理论来说不应该加这个，但不加的话数值会不收敛
# 肯定有控制收敛的办法，不过我现在并不知道
for i in range(0,tnum):
    idx = i * nnum
    kmat[idx,:] = 0
    kmat[idx,idx] = 1
    bmat[idx] = 0
    idx = idx + nnum - 1
    kmat[idx,:] = 0
    kmat[idx,idx] = 1
    bmat[idx] = 0.1    
        
#偷懒用库函数求逆矩阵
kinvmat = np.linalg.inv(kmat)

for i in range(0,nmax):
    temp = 0
    for j in range(0,nmax):
        temp += kinvmat[i,j]*bmat[j]
    res[(int)(i/nnum),(int)(i%nnum)] = temp
