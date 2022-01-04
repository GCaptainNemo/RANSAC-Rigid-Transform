# 刚体变换 RANSAC

## 一、介绍

随机一致性采样(Random Sample Consensus, RANSAC)算法是一种在存在噪声条件下通过采样提取内点(inliers)，并计算模型的算法。

传统求解噪声下的模型常把它转换为一个最小二乘问题，构造关于模型参数的损失函数，在**二范数**下最小化损失函数，求解模型参数。从这个角度上说，RANSAC算法追求的是**零范数**下最小化损失函数(即符合模型的内点个数最大)。

RANSAC具体算法如下：

1. 从包含噪声的集合**X**中以均匀分布采放回采m个样品(m>=模型求解所需最小样品个数)。
2. 根据m个样品计算模型(二范数下最小)。
3. 设置阈值e，计算集合**X**中与模型的偏差小于阈值e的点数(内点)。
4. 重复1，2，3。找到内点数最大的模型。

令点集中内点所占的比例(概率)为d，得到正确模型的概率为p。则采样次数至少为 

<p align="center">N = log(1 - p) / log(1 - d <sup>M</sup>)</p>

本仓库实现3d-3d匹配点下求解刚体变化模型下的RANSAC提取算法，模型求解方法见[ICP算法](https://github.com/GCaptainNemo/ICP-registration-PCL/blob/main/docs/ICP.ipynb)。



## 二、依赖和使用

* 依赖

`Eigen3`

* 使用方法

```
mkdir build & cd build
cmake ..
make
./RANSAC_RIGID_TRANSFORM
```







