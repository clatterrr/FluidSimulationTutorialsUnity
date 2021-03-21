import numpy as np
#Runge Kutta四阶
nmax = 20
h = np.zeros((6,nmax))
k = np.zeros((nmax))

for i in range(nmax):
    h[0,i] = i*i*i*i
    
# 因为没有处理边界条件，所以不会把整列都算上
# 计算完成后，h的最后一行的第8列到第11列
# 正好就是h的第一行的第7列到第10列
for i in range(2,nmax-2):
    #k[i] = 3*i*i
    k[i] = (h[0,i-2] - 8*h[0,i-1] + 8*h[0,i+1] -h[0,i+2])/12
    h[1,i] = h[0,i] - k[i]
    
for i in range(4,nmax-4):
    k[i] = (h[1,i-2] - 8*h[1,i-1] + 8*h[1,i+1] -h[1,i+2])/12
    h[2,i] = h[1,i] - k[i]
    
for i in range(6,nmax-6):
    k[i] = (h[2,i-2] - 8*h[2,i-1] + 8*h[2,i+1] -h[2,i+2])/12
    h[3,i] = h[2,i] - k[i]
    
for i in range(8,nmax-8):
    k[i] = (h[3,i-2] - 8*h[3,i-1] + 8*h[3,i+1] -h[3,i+2])/12
    h[4,i] = h[3,i] - k[i]
    h[5,i] = h[4,i]/24 + h[2,i]/4 + h[1,i]/3 + h[0,i]*9/24
    
'''
\a  b  c  d  e\ \1\ = \1  \
\0  b 2c 3d 4e\ \1\ = \1  \
\0  0  c 3d 6e\ \1\ = \1/2\
\0  0  0 d  4e\ \1\ = \1/6\
\0  0  0 0   e\ \1\ = \1/24\

e = 1/24
d = 0
c = 1/4
b = 1/3
a = 9/24

'''
    