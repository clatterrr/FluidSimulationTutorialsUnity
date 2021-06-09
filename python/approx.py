###拟合年龄

import numpy as np
import matplotlib.pyplot as plt
#定义x、y散点坐标
x = [1,2,3,4,5,6]
x = np.array(x)
print('x is :\n',x)
num = [1,1,1,1,1,20]
y = np.array(num)
print('y is :\n',y)
#用3次多项式拟合
f1 = np.polyfit(x, y, 5)
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
nmax = 60
res = np.zeros((nmax))
newx = np.zeros((nmax))
for i in range(0,nmax):
    newx[i] = i / (nmax / 6) + 1
    res[i] = f1[0]*newx[i]**5 + f1[1]*newx[i]**4 + f1[2]*newx[i]**3 + f1[3]*newx[i]**2 + f1[4]*newx[i] + f1[5]
plot1 = plt.plot(newx, res, 's',label='original values')   
plt.show() 