# 简介

游戏流体力学基础及unity代码

github地址https://github.com/clatterrr/FluidSimulationTutorialsUnity

gitee码云地址https://gitee.com/clatterrr/FluidSimulationTutorialsUnity

qq模拟流体交流群1001290801，欢迎加入

代码作者：光影帽子

# 地址

【游戏流体力学基础及Unity代码（一）】热传导方程

https://zhuanlan.zhihu.com/p/263053689

【游戏流体力学基础及Unity代码（二）】有限差分法

https://zhuanlan.zhihu.com/p/264153771

【游戏流体力学基础及Unity代码（三）】用波动方程模拟三维落雨池塘，连续性方程

https://zhuanlan.zhihu.com/p/264585002

![Alt Text](images\ch03WaveEquation.PNG)

【游戏流体力学基础及Unity代码（四）】用欧拉方程模拟无粘性染料之公式推导

https://zhuanlan.zhihu.com/p/270530827

【游戏流体力学基础及Unity代码（五）】用欧拉方程模拟无粘性染料之代码实现

https://zhuanlan.zhihu.com/p/270531017

![Alt Text](images\ch05EulerEquation.png)

【游戏流体力学基础及Unity代码（六）】用NavierStokes方程模拟粘性染料流动

https://zhuanlan.zhihu.com/p/283662524

![Alt Text](images\ch06NavierStokes.png)

【游戏流体力学基础及Unity代码（七）】车流量问题，非线性水波以及burgers方程

https://zhuanlan.zhihu.com/p/309860521

【游戏流体力学基础及Unity代码（八）】有限体积法

https://zhuanlan.zhihu.com/p/331771977

【游戏流体力学基础及Unity代码（九）】用浅水波方程模拟雨落池塘和DamBreak

https://zhuanlan.zhihu.com/p/331781508

![Alt Text](images\ch09ShallowWater.png)

B站视频https://www.bilibili.com/video/BV1Ry4y167MV

【游戏流体力学基础及Unity代码（十）】漩涡和模拟二维烟雾

https://zhuanlan.zhihu.com/p/340842666

![Alt Text](images\ch10smoke2d.png)

【游戏流体力学基础及Unity代码（十一）】理想流体机翼绕流和升力原理

https://zhuanlan.zhihu.com/p/340848576

![Alt Text](images\ch11SourceSink.png)

【游戏流体力学基础及Unity代码（十二）】卡门涡街，边界层，涡方法

https://zhuanlan.zhihu.com/p/345332340

 ![Alt Text](images\ch12VorteXStreet.png)

B站视频https://www.bilibili.com/video/BV1u5411H7hr

【游戏流体力学基础及Unity代码（十三）】泊松压力方程，SIMPLE算法

https://zhuanlan.zhihu.com/p/347410166

【游戏流体力学基础及Unity代码（十四）】舌尖上的有限元Galerkin法

https://zhuanlan.zhihu.com/p/358033368

【游戏流体力学基础及Unity代码（十五）】线性有限元及弹性物体模拟

https://zhuanlan.zhihu.com/p/369505527

波前推进法网格生成https://www.bilibili.com/video/BV1ZK4y1w7R6/

![Alt Text](images/ch15meshgen.png)

弹性果冻模拟https://www.bilibili.com/video/BV1w84y1c7K2/

![Alt Text](images/ch15jelly.png)

【游戏流体力学基础及Unity代码（十六）】非线性有限元及牛顿迭代法

https://zhuanlan.zhihu.com/p/369521901

# 收集

目前代码和视频仍然会尽量更新，但文字版教程将尽量放缓，是为了确保教程质量，以及不犯低级错误。所以你可能会发现一大堆没对应教程的代码，这些宝贵的代码收集起来很不容易，我会尽量保证代码可读性，以及我是从哪里得到的。不过一般我下载到的都是matlab或c++的，我会把它转写成python以加深理解。

以下是我收集的一些代码和网址，觉得很不错就贴上来，不定时更新。不过这些网址可能随时会挂掉

### 浸入边界法Immersed Boundary Method

https://sites.google.com/view/sglee/research 进入页面，搜索code，找到“\21. Wanho Lee and Seunggyu Lee, Immersed boundary method for simulating interfacial problems, Mathematics 8(11) (2020) 1982 [..](https://drive.google.com/file/d/1aeJv_8TqKAmYVYxu9R36yhnCNwJFoxFd/view?usp=sharing) (IF2019:1.747) ([code](https://drive.google.com/file/d/1zwRWtJG8cu6lqRBnzYHYbpWHuZNfMmFU/view?usp=sharing))”字样，点击code即可下载。这个页面上还有许多文章是可免费下载的

https://github.com/nickabattista/IB2d 很棒的开源库，有代码，论文和视频

### 线性方程组的迭代解法

http://www.netlib.org/templates/matlab/ 就是一些共轭梯度，预处理的共轭梯度，最小残差GMRES，双共轭梯度

### 网格生成

http://persson.berkeley.edu/ 特别棒的二维三维网格生成的matlab代码

### 非定常流

https://www.mathworks.com/matlabcentral/fileexchange/?q=profileid:4187051

### 杂

https://people.sc.fsu.edu/~jburkardt/m_src/ 有各种各样的有限元matlab代码

