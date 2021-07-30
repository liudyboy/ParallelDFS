# ParallelDFS
Deep first search (DFS) parallel implementatoin based on MPI.

local_n :   进程内所包含的边的数目

1. local path search
对于每条边都作为起始边进行路径搜索     --------->   得到local_n 条路径

搜索所得数据存储格式


|前local_n| 边信息|
|----|----|
|前local_n 个整数保存以该边为起始边可形成的路径个数|	<end_edge, local_path_len>|

2. global path search
拥有全局的起始边和终点边的进程开始进行全局的搜索

全局搜索路径表示
```
[global_start_edge, end_edge, local_path, end_edge, local_path, .....,  global_end_edge, local_path]
```

3. path sequence search
各个进程基于`<start_edge, end_edge, len>`三元组，进行路径序列的搜索。


![https://github.com/liudyboy/ParallelDFS/tree/main/imgs/example3.PNG]()

![](https://github.com/liudyboy/ParallelDFS/tree/main/imgs/example4.PNG)
