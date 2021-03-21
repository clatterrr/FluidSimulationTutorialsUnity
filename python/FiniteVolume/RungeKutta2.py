import numpy as np
#Runge Kutta二阶
nmax = 16
h = np.zeros((4,nmax))
k = np.zeros((nmax))

for i in range(nmax):
    h[0,i] = i*i
    
# 因为没有处理边界条件，所以不会把整列都算上
# 计算完成后，h的最后一行的第4列到11列
# 正好就是h的第一行的第3列到第10列
for i in range(2,nmax-2):
    #k[i] = 3*i*i
    k[i] = (h[0,i-2] - 8*h[0,i-1] + 8*h[0,i+1] -h[0,i+2])/12
    h[1,i] = h[0,i] - k[i]
    
for i in range(4,nmax-4):
    k[i] = (h[1,i-2] - 8*h[1,i-1] + 8*h[1,i+1] -h[1,i+2])/12
    h[2,i] = h[1,i] - k[i]
    h[3,i] = h[2,i]/2 + h[0,i]/2

'''

\a  b  c\ \1\ = \1  \
\0  b 2c\ \1\ = \1  \
\0  0  c\ \1\ = \1/2\  

c = 1/2 
a = 1/2
'''
    

    