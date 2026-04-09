# PyImageJ-Script-Examples

本仓库包含我所实现的一些代码示范，它们主要的目的是通过PyImageJ调用Fiji以分析图像。

本仓库大概率不会为每个脚本提供示例数据，除非这个脚本是基于某篇公开文献的实验复现。

本项目所有脚本都遵循如下风格：

- 使用PyImageJ调用Fiji以载入图像，并遵循ImageJ2的风格
- 核心分析步骤通常采用Numpy或SkImage库实现，而非全部依赖Fiji
