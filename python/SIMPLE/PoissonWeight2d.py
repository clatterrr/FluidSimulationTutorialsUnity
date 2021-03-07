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
n = np.ones((nmax+2,mmax+2))
s = np.ones((nmax+2,mmax+2))
w = np.ones((nmax+2,mmax+2))
e = np.ones((nmax+2,mmax+2))
n[2,2] = 5
s[2,3] = 5
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
            term = (e[i,j]*p[i+1,j] + w[i,j]*p[i-1,j])*dy*dy + (n[i,j]*p[i,j+1] + s[i,j]*p[i,j-1])*dx*dx
            p[i,j] = (term - dx*dx*dy*dy*div[i,j])/((e[i,j] + w[i,j])*dy*dy + (n[i,j] + s[i,j])*dx*dx)
for i in range(0,nmax+1):
    for j in range(0,mmax+2):
        u2[i,j] = u[i,j] + (p[i,j] - p[i+1,j])/dx*e[i,j]
for i in range(0,nmax+2):
    for j in range(0,mmax+1):
        v2[i,j] = v[i,j] + (p[i,j] - p[i,j+1])/dy*n[i,j]
for i in range(1,nmax+1):
    for j in range(1,mmax+1):
        div2[i,j] = (u2[i,j] - u2[i-1,j])/dx + (v2[i,j] - v2[i,j-1])/dy