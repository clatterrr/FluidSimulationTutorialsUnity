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

之后预计很长一段时间不会更新了。是为了确保教程质量，以及不犯低级错误。所以你可能会发现一大堆没对应教程的代码，这些宝贵的代码收集起来很不容易，我会尽量保证代码可读性，以及我是从哪里得到的。不过一般我下载到的都是matlab或c++的，我会把它转写成python以加深理解。

以下是我收集的一些代码和网址，觉得很不错就贴上来，不定时更新。不过这些网址可能随时会挂掉

### 浸入边界法Immersed Boundary Method

https://sites.google.com/view/sglee/research 进入页面，搜索code，找到“\21. Wanho Lee and Seunggyu Lee, Immersed boundary method for simulating interfacial problems, Mathematics 8(11) (2020) 1982 [..](https://drive.google.com/file/d/1aeJv_8TqKAmYVYxu9R36yhnCNwJFoxFd/view?usp=sharing) (IF2019:1.747) ([code](https://drive.google.com/file/d/1zwRWtJG8cu6lqRBnzYHYbpWHuZNfMmFU/view?usp=sharing))”字样，点击code即可下载。这个页面上还有许多文章是可免费下载的

https://github.com/nickabattista/IB2d 很棒的开源库，有代码，论文和视频

https://www.math.nyu.edu/~peskin/ib_lecture_notes/index.html

https://github.com/shurikkuzmin/ImmersedBoundary

### 半拉格朗日Semi-Lagrange

https://github.com/iCFD/SemiLagrangian

D:\FluidSim\FluidSim\semilagrange\SemiLagrangian-master

https://github.com/abarret/SemiLagrangian

### Fluid Solid Interaction

An Introduction to Fluid-Structure Interaction: Application to the Piston Problem  

项目地址 ： http://www.utc.fr/~elefra02/ifs/

代码地址：http://www.utc.fr/~elefra02/ifs/archive_FSI.tar.gz

本地代码：D:\FluidSim\FluidSim\FluidSolidInteraction\archive_FSI

https://github.com/WhiteTshirtXI/IBFS_M

### 有限元

Efficient implementation of adaptive P1-FEM in Matlab

https://www.pplusplus.lima-city.de/femfluid.html Pressure Solve with Finite Elements  很好的matlab库

代码地址：https://www.pplusplus.lima-city.de/lib/data/femfluid/FEM%20Fluid%20Source.zip

本地地址：D:\FluidSim\FluidSim\FEMNEW\FEM Fluid Source\FEM Fluid

https://www.math.hu-berlin.de/~cc/cc_homepage/software/software.shtml 

Computationally Solving Nonlinear Membranes with Plane Stress Condition  

https://github.com/vasko6d/finite-element-solver

D:\FluidSim\FluidSim\FEMGOOD\finite-element-solver-master\finite-element-solver-master

H^1-Stability of the L^2-Projection onto Finite Element Spaces on Adaptively Refined Quadrilateral Meshes

https://github.com/aschmidtuulm/h1-stability

https://github.com/Vinay5SVeerapur/Finite-element-analysis/blob/master/BEAM%20equation.ipynb

https://github.com/tobyvg/Fluid-codes 方强流

https://github.com/Milad-Rakhsha/FEM_PDE 有限元解势流

https://github.com/emarinhoss/FEM_PETSC

https://github.com/RnkSngh/Double-Slit-Experiment 双缝

https://github.com/shardoolk/FEM

https://github.com/jborggaard/ns2d

https://github.com/Satchit4/Navier-stokes

https://github.com/coltonjconroy/DG_2d_lava_flows

https://github.com/nileshjchoudhary/Flow-through-driven-cavity-Finite-element-analysis-CFD

https://github.com/michelrobijns/pyBurgersFEM

https://github.com/Hahany/Finite-element-method

### Galerkin

An Introduction to Element-based Galerkin Methods on Tensor-Product Bases: Analysis, Algorithms, and Applications

https://github.com/fxgiraldo/Element-based-Galerkin-Methods

D:\FluidSim\FluidSim\Galerkin\Element-based-Galerkin-Methods-master\Element-based-Galerkin-Methods-master

https://github.com/tuhouwang?tab=repositories

Nodal Based Galerkin

https://github.com/Achyut2404/nodalDG

D:\FluidSim\FluidSim\Galerkin\nodalDG-master\src

https://github.com/Jacklswalsh/DGM-Advection-AD 自适应快速一维Galerkin

D:\FluidSim\FluidSim\Galerkin\DGM-Advection-AD-main

https://github.com/asdf123101/HDPG1D

D:\FluidSim\FluidSim\Galerkin\HDPG1D-master

https://github.com/hanveiga/higher-order-methods/blob/master/dg1d.py

https://github.com/wme7/cprlinearexamples 画正方形

https://github.com/AndrewWang996/Discontinuous-Galerkin

https://github.com/Chang-Liu-0520/1D_advec_DG

https://github.com/pinkieli/Interpolation-Nodes-for-High-order-Lagrange-Finite-Elements.Nodal Discontinuous Galerkin Methods: Algorithms, Analysis, and Applications", Jan S Hesthaven and Tim Warburton.

### 谱方法

SPECTRAL METHOD FOR TIME DEPENDENT NAVIER-STOKES
EQUATIONS  

http://cpraveen.github.io/teaching/chebpy.html

### 数值积分

https://github.com/jgressier/flowdyn/blob/master/flowdyn/integration.py RK4 low storage

### 边界元

https://team-pancho.github.io/deltaBEM/download.html

### 线性方程组的迭代解法

http://www.netlib.org/templates/matlab/ 就是一些共轭梯度，预处理的共轭梯度，最小残差GMRES，双共轭梯度

https://github.com/JuliaLinearAlgebra/IterativeSolvers.jl Julia语言实现的各种方程组解法。

*Optimization in Practice with MATLAB®: For Engineering Students and Professionals* 最优化书籍

https://github.com/Manchery/numerical-analysis-practice

https://github.com/JordanFisher/Paper-Implicit-IBM-2D/blob/master/CODE%20FREEZE/NewHeartValveSim.py 有预处理共轭梯度，快速雅可比等

### 网格生成

http://persson.berkeley.edu/ 特别棒的二维三维网格生成的matlab代码

https://ifsnumericaltools.weebly.com/ 也是很棒的代码 D:\FluidSim\MathsWorkMisc\mesh2d_v24\Mesh2d_v24

https://github.com/aschmidtuulm/ameshref Adaptive Mesh Refinement in 2D–An Efficient Implementation in Matlab论文对应的代码

TetGen http://wias-berlin.de/software/tetgen/formAction12.jsp

D:\FluidSim\OpenSource\tetgen1.5.1\tetgen1.5.1

https://doc.cgal.org/4.13/Manual/tutorials.html



### 非定常流

https://www.mathworks.com/matlabcentral/fileexchange/?q=profileid:4187051

### 声学

Physically Based Sound for Computer Animation and Virtual Environments

[ACM SIGGRAPH 2016 Course](http://s2016.siggraph.org/courses/events/physically-based-sound-computer-animation-and-virtual-environments)

### 水平集

教授Osher Stanley,

*Geometric Level Set Methods in Imaging,Vision & Graphics*

https://www.cs.ubc.ca/~mitchell/ToolboxLS/ matlab示例，非常棒的成系统的代码。

https://github.com/scikit-image/scikit-image/blob/main/skimage/segmentation/_chan_vese.py 有一个有名的python库叫scikit-image，里面实现了chanvese算法

*A discrete level-set topology optimization code written in Matlab*

### Closet Point Method

https://www.math.ubc.ca/~cbm/cpm/

### 多孔介质

An Introduction to the Numerics of Flow in Porous Media using Matlab  

https://github.com/pmgbergen/porepy

https://github.com/jjhidalgo/HGCchem2 有分层现象，大佬主页https://jjhidalgo.wordpress.com/codes/

### 浅水波

项目地址：https://web.cse.ohio-state.edu/~wang.3602/courses/cse3541-2019-fall/index.html

unity 包：https://web.cse.ohio-state.edu/~wang.3602/courses/cse3541-2019-fall/lab4/wave_example.unitypackage

### 多重网格

http://pages.cs.wisc.edu/~sifakis/project_pages/mgpcg.html 并行多重网格泊松求解器附代码

https://github.com/JuliaLinearAlgebra/AlgebraicMultigrid.jl/tree/master/src 代数多重网格

https://github.com/danfortunato 直接关注这个人就行了

https://github.com/pymatting/pymatting/blob/master/pymatting/preconditioner/vcycle.py

https://github.com/lyc102/ifem

https://amgcl.readthedocs.io/en/latest/examples.html

https://github.com/evstigneevnm/GMG_2D_tests

https://github.com/gnitish18/FEM_Multigrid

### 多相流

2D Cartesian Quadtree Adaptive Mesh Refinement (AMR) for multiphase Five Equations Model.https://github.com/dattv/2D_CARFIVE

https://github.com/Spoonacular/LBM_python

https://github.com/mirsandiharyo/multiphase_flows_front_tracking_python 模拟泡泡和水滴

https://github.com/rarbarim/multiphase_flow_simulator 附带报告，作者还有一些别的代码

### 数学

http://pages.cs.wisc.edu/~sifakis/project_pages/svd.html Computing the Singular Value Decomposition of 3x3 matrices with minimal branching and elementary floating point operations附代码

### 泡泡

Role of the Dynamic Contact Angle on Splashing  

### 湍流

DNSLABhttps://md-datasets-cache-zipfiles-prod.s3.eu-west-1.amazonaws.com/6gtnjwwg8j-1.zip

LESCODE https://cfd.engr.uconn.edu/ 不过代码需要填个申请表格才能获取

本地代码地址：D:\FluidSim\FluidSim\LES\les.r123\les

Turbulent Fluids – SIGGRAPH Course https://ge.in.tum.de/research/turbulent-fluids-siggraph-course/

https://github.com/thijsbon/CMF_project_thijs_victor 有湍流和墙函数

https://github.com/Maikuelet/Turbulence_Modelling_Burgulence

### 边界元

https://github.com/Timmmdavis/CutAndDisplace

https://github.com/Timmmdavis/CutAndDisplace

### 碰撞检测

Fast Continuous Collision Detection using Deforming Non-Penetration Filters

项目地址：http://gamma.cs.unc.edu/DNF/

代码地址：http://gamma.cs.unc.edu/DNF/request.html

### 形状优化

https://github.com/jorgensd/MultiMeshShapeOpt_code

### 电磁学

A generalized polynomial chaos based ensemble Kalman filter with high accuracy  

https://github.com/Andrewpensoneault/Lorenz_63_Stochastic_Galerkin_EnKF

https://github.com/joebling/graduate_essay

https://github.com/Mjjnuu/DoublePendulum 似乎是宇宙学

D:\FluidSim\Electron\DoublePendulum-master\python

https://github.com/keileg/fvbiot

https://github.com/tarcisiofischer/helmholtz-solver/tree/master/src/python

https://github.com/ep2lab 一些磁流体

https://github.com/mgoycoolea/twofluid/blob/master/twofluid.py

https://github.com/trevorcrupi/EM-MG 电磁学的多重网格

https://github.com/rasalkumar/FEM

https://github.com/ocramz/lib_FEM_py

https://github.com/eduardobehr/pyjoule

### 声学

http://www.k-wave.org/download.php

https://github.com/pvanberg/DGFEM-Acoustic

https://github.com/1ceaham/AcousticFVTD_GeneralImpedance

https://github.com/ivanmartinezsuarez/Matlab_FVM

### 粘弹性

https://github.com/labsin-unesp/Viscoel-stico-Kelvin-Voigt

### 杂

https://people.sc.fsu.edu/~jburkardt/m_src/ 有各种各样的有限元matlab代码

https://github.com/weymouth/WaterLily.jl 漂亮的NS方程模拟

Extraction of Distinguished Hyperbolic Trajectories for 2D Time-Dependent Vector Field Topology

介绍页面：https://vcg.iwr.uni-heidelberg.de/people/sadlo/

代码地址：https://github.com/lhofmann/eurovis2020_hyperbolic_trajectories

A PArallel Robust Interface Simulator that combines VOF and Front-Tracking

介绍页面：http://www.ida.upmc.fr/~zaleski/paris/index.html

代码地址：http://www.ida.upmc.fr/~zaleski/paris/paris-stable.tar.gz

A Hyperbolic Geometric Flow for Evolving Films and Foams

项目地址：https://ryichando.graphics/

代码地址：https://github.com/sdsgisd/HGF

Interpolation Nodes for High-order Lagrange Finite Elements

https://github.com/pinkieli/Interpolation-Nodes-for-High-order-Lagrange-Finite-Elements.

Semi-Riemannian Manifold Optimization

https://github.com/trgao10/SemiRiem

https://github.com/noamaig/hyperbolic_orbifolds

An entropy-stable hybrid scheme for simulations of transcritical real-fluid flow[JCP的]

https://github.com/peterma123456789/DoubleFlux-1D

Ice sheet flow with thermally activated sliding

https://github.com/elisamantelli/subtemperate_sliding_rspa_2019

Wavelet-Fourier CORSING techniques for multi-dimensional advection-diffusion-reaction equations

https://github.com/simone-brugiapaglia/corsing-wavelet-fourier-adr

River Profile

https://github.com/sfgallen/ChiProfiler

https://github.com/ISSI2015/M4 Real-Time Deformation

### 超声速可压缩

https://github.com/Fanxiaotsing/One-dimensional-aero-heating-code 平板对流换热

https://github.com/holdmygithub/ASOInviscidSupersonicFlow 机翼设计

D:\FluidSim\FluidSim\NavierStokes\ASOInviscidSupersonicFlow-master\ASOInviscidSupersonicFlow-master

https://github.com/amikkonen/lidDrivenCavityCompressibleFlowPython 可压缩顶盖驱动

https://github.com/GerardBoberg/CompressiblePipeFlow

D:\FluidSim\FluidSim\CompressibeNewgood\CompressiblePipeFlow-master\CompressiblePipeFlow-master

https://home.cscamm.umd.edu/centpack/examples/euler2d.htm#press

https://github.com/silentmovie/RTmodel RT不稳定

PYRO2 开源库https://python-hydro.github.io/pyro2/compressible_basics.html

超级好https://github.com/jingchangshi/NumericalMethodsForConservationLawsDG

### 空气动力

https://github.com/Maikuelet/FEM_Airplane

### SIMPLE/PISO

https://github.com/mehrdadyo/LS-IBM

### 混合网格粒子法

PolyPIC: the Polymorphic-Particle-in-Cell Method for Fluid-Kinetic Coupling

https://github.com/smarkidis/fluid-kinetic-PIC

### 大佬主页

下面的主页全部是附有开源代码的

https://cs.uwaterloo.ca/~c2batty/

https://zhxx1987.github.io/#cod 猜猜这是谁？

http://gamma.cs.unc.edu/software/ 这是个项目主页，开源代码很多

https://www.cc.gatech.edu/~turk/

https://people.llnl.gov/lindstrom2 偏向几何数据处理

http://ntoken.com/pubs.html#Thuerey_2016_ofblend

http://www.tkim.graphics/

http://www.cmap.polytechnique.fr/~allaire/

https://sites.google.com/view/valentinresseguier/projects

https://www.konrad-simon.eu/wordpress/?page_id=91

### 开源项目

pyro2

netgen

scipy这玩意有一些矩阵迭代求解法

http://granoo.52083.n8.nabble.com/

### 不错的论文

内容很棒的论文以及讲义

Lecture notes Introduction to numerical methods for interfacial flows  

### 很有个性的论文标题

内容不管，但是标题值得写上一万字来吐槽

A massive fractal in days, not years

