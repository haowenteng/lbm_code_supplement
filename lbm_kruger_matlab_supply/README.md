# 中文简介
本项目基于MATLAB，对Krüger的经典LBM书籍《The Lattice Boltzmann Method: Principles and Practice》进行代码补充。原作者项目中包含了基于c++的各章节代码，本项目主要出于入门学习与比较的目的，补充了第八章的matlab形式代码，此外，补充了动画绘制部分，可直接生成结果动画，仅供参考。

原cpp代码中存在部分简化，如碰撞后未进行传播等，阅读时应注意。

注意！部分代码错误待修正。

## 1. cylinder气缸扩散
无流气缸扩散，一个由圆柱形空腔组成的系统，该空腔填充有静止液体。腔内化学物质的初始浓度，圆柱体表面的浓度保持恒定。随着时间的推移，物质从圆柱体表面扩散到液体中。
D2Q9模型
## 2. film平板扩散
- antibb
两个平行板之间流的浓度层发展。板表面的浓度是恒定的。
- inamuro
半无限板问题。底部边界保持在恒定的板浓度。假设背景流是均匀的。
- uniform
半无限板问题。
## 3. Gaussian hill扩散
具有初始高斯分布浓度的物种在均匀速度场的存在下的发展。观察到由于非零速度引起的平流和由于不均匀浓度引起的扩散。
- 1d_bgk
- 1d_magic12
- 2d_bgk
- 2d_trt

---
# Introduction
This repository is based on Matlab language and supplements the modeling of the case in Chapter 8 of the book

*The Lattice Boltzmann Method: Principles and Practice*  
T. Krüger, H. Kusumaatmaja, A. Kuzmin, O. Shardt, G. Silva, E.M. Viggen  
authors@lbmbook.com  
ISBN 978-3-319-44649-3 (Electronic) 978-3-319-44647-9 (Print)  
http://www.springer.com/978-3-319-44647-9  

All code is shared for beginner learning LBM.