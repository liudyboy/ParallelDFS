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
  // MPI_Init(&argc, &argv);
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  ParallelDFS parallel_dfs;
  // parallel_dfs.SetOutDir("path2/");
  parallel_dfs.run();
  MPI_Finalize();
  return 0;
}
