import numpy as np
#仅仅是将正方形拆成左下和右上两个三角形，网格很规整

#初始条件
Totalcol = 8
Totalrow = 8
dx = 1
dy = 1
coff1 = 1
coff2 = 2
coff3 = -2
#c1x + c2y + c3 = 0
a = 1
b = 1 # - coff1 - coff2

#接下来是计算机的事情
trinum = Totalcol*Totalrow*2
pnum = (Totalcol + 1)*(Totalrow + 1)
tri = np.zeros((trinum,6))#每个三角形三个顶点，每个顶点有个二维坐标xyxyxy
triidx = np.zeros((trinum,3))#三角形三个顶点的编号，不同三角形共用
for i in range(0,trinum):
    '''
    0/0 -------- 2
    \            \
    \            \
    \            \
    2 -----------1/1
    '''
    idx = int(i/2)
    rev = int(i%2)
    row = int(idx%Totalrow)
    col = int(idx/Totalrow)
    tri[i,0] = row*dx
    tri[i,1] = col*dy + rev*dy
    tri[i,2] = row*dx + dx
    tri[i,3] = col*dy
    tri[i,4] = row*dx + rev*dx
    tri[i,5] = col*dy + dy
    triidx[i,0] = col*(Totalrow + 1) + row + rev*(Totalrow + 1)
    triidx[i,1] = col*(Totalrow + 1) + row + 1
    triidx[i,2] = col*(Totalrow + 1) + Totalrow + 1 + row + rev*1

kmat = np.zeros((pnum,pnum))
bmat = np.zeros((pnum))
res = np.zeros((pnum))
U = np.zeros((Totalrow + 1,Totalcol + 1))
for i in range(0,trinum):
    b2 = tri[i,5] - tri[i,1]
    b3 = tri[i,1] - tri[i,3]
    c2 = tri[i,0] - tri[i,4]
    c3 = tri[i,2] - tri[i,0]
    Area2 = abs(b2*c3 - c2*b3)
    idx0 = int(triidx[i,0])
    idx1 = int(triidx[i,1])
    idx2 = int(triidx[i,2])
    kmat[idx0,idx0] += (-1/6)*(b2 + a*c2)/Area2 + (-1/6)*(b3 + a*c3)/Area2
    kmat[idx0,idx1] += (1/6)*(b2 + a*c2)/Area2
    kmat[idx0,idx2] += (1/6)*(b3 + a*c3)/Area2
    bmat[idx0] += (-b/6)
for i in range(0,Totalrow+1):
    for j in range(0,Totalcol+1):
        idx = j*(Totalrow + 1) + i
        if(bmat[idx] == 0):
            kmat[idx,:] = 0
            kmat[idx,idx] = 1
            bmat[idx] = coff1*i + coff2*j + coff3
kinv = np.linalg.inv(kmat)
for i in range(0,pnum):
    summ = 0
    for j in range(0,pnum):
        summ += kinv[i,j]*bmat[j]
    res[i] = summ  
for i in range(0,Totalrow+1):
    for j in range(0,Totalcol+1):
        U[i,j] = res[j*(Totalrow + 1) + i]
    
    
    