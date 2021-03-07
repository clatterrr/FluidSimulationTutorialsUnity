import numpy as np
#压力泊松方程交错网格二维示例，作者：光影帽子
nmax = 4#宽X轴
mmax = 6#高Y轴
tmax = 1000
u = np.zeros((nmax+1,mmax+2))
v = np.zeros((nmax+2,mmax+1))
p = np.zeros((nmax+2,mmax+2))
div = np.zeros((nmax+2,mmax+2))
dx = 1/nmax
dy = 1/mmax
u2 = np.zeros((nmax+1,mmax+2))
v2 = np.zeros((nmax+2,mmax+1))
div2 = np.zeros((nmax+2,mmax+2))
for i in range(1,nmax):
    for j in range(1,mmax+1):
        u[i,j] = j
for i in range(1,nmax+1):
    for j in range(1,mmax):
        v[i,j] = i
for i in range(1,nmax+1):
    for j in range(1,mmax+1):
        div[i,j] = (u[i,j] - u[i-1,j])/dx + (v[i,j] - v[i,j-1])/dy
for t in range(0,tmax-2):
    for i in range(1,nmax+1):
        for j in range(1,mmax+1):
            term = (p[i+1,j] + p[i-1,j])*dy*dy + (p[i,j+1] + p[i,j-1])*dx*dx - dx*dx*dy*dy*div[i,j]
            p[i,j] = p[i,j] + (term/(2*(dx*dx + dy*dy)) - p[i,j])
for i in range(0,nmax+1):
    for j in range(0,mmax+2):
        u2[i,j] = u[i,j] + (p[i,j] - p[i+1,j])/dx
for i in range(0,nmax+2):
    for j in range(0,mmax+1):
        v2[i,j] = v[i,j] + (p[i,j] - p[i,j+1])/dy
for i in range(1,nmax+1):
    for j in range(1,mmax+1):
        div2[i,j] = (u2[i,j] - u2[i-1,j])/dx + (v2[i,j] - v2[i,j-1])/dy