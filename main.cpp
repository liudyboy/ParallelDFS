/*
 * @Brief: file description
 * @Author: liudy
 * @Email: deyin.liu@nscc-gz.cn
 * @Date: 2021-07-27 15:14:36
 * @LastEditors: liudy
 * @LastEditTime: 2021-07-27 16:04:55
 */
#include <iostream>
#include "parallel_dfs.h"

using namespace std;


int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Barrier(MPI_COMM_WORLD);
  ParallelDFS parallel_dfs;
  parallel_dfs.run();
  MPI_Finalize();
  return 0;
}