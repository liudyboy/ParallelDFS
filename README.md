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
[global_start_edge, end_edge, path_length, end_edge, path_length, .....,  global_end_edge, path_length]
```
其中的path_length 表示的是当前路径已有的长度。

3. path sequence search
各个进程基于`<start_edge, end_edge, len>`三元组，进行路径序列的搜索。

## Example 3: circle exit among processes

![example3](./imgs/example3.png)

+ P0: 0 1 2
+ P1: 3 4 5
+ P2: 6 7 8 9


## Example 4: circle exit in a process

![example4](./imgs/example4.png)

+ P0: 0 1 2 3 4
+ P1: 5 6 7 8 9


## Memory-constrain implementation
从拥有全局起始点的进程出发，一条一条边的进行延长，在各个进程本地会保存一个信息列表格式如下
```c++
<path_id, start_edge, count>
```
其中`path_id`的数据格式为前32位保存`process_id`后32位保存`local_path_count`。

针对每个start_edge 进行搜索是还需要两个变量进行传递`<cur_fork, next_fork>`分别表示当前搜索在那个岔路可以进行变道和下一个可以进行变道的岔路口（都是edge的编号）。
