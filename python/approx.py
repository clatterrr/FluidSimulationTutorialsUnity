import numpy as np
import matplotlib.pyplot as plt
#定义x、y散点坐标
nx = 3
x = np.array([0,1,2])
print('x is :\n',x)
num = [1,1.125,1.875]
y = np.array(num)
print('y is :\n',y)
#用3次多项式拟合
f1 = np.polyfit(x, y, 2)
print('f1 is :\n',f1)

p1 = np.poly1d(f1)
print('p1 is :\n',p1)
#也可使用yvals=np.polyval(f1, x)
yvals = p1(x) #拟合y值
print('yvals is :\n',yvals)
#绘图
plot1 = plt.plot(x, y, 's',label='original values')
plot2 = plt.plot(x, yvals, 'r',label='polyfit values')
plt.xlabel('x')
plt.ylabel('y')
plt.title('polyfitting')
plt.show()
nmax = nx * 10
res = np.zeros((nmax))
newx = np.zeros((nmax))
for i in range(0,nmax):
    newx[i] = i / 10 + x[0]
    for j in range(0,nx):
        res[i] += f1[nx-j-1]*(newx[i]**j)
plot1 = plt.plot(newx, res, 's',label='original values')   
plt.show() 