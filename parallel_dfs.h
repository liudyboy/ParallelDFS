/*
 * @Brief: file description
 * @Author: liudy
 * @Email: deyin.liu@nscc-gz.cn
 * @Date: 2021-07-27 15:14:09
 * @LastEditors: liudy
 * @LastEditTime: 2021-07-28 09:28:40
 */
#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "mpi.h"

#define FILE_CACHE_LIMIT 2
#define LOCAL_LENGTH_LIMIT 262144  // 2 MB (unsigned long)

using std::ofstream;
using std::string;
using std::vector;
#define mix(h)                    \
  ({                              \
    (h) ^= (h) >> 23;             \
    (h) *= 0x2127599bf4325c37ULL; \
    (h) ^= (h) >> 47;             \
  })

uint64_t fasthash(const void* buf, size_t len, uint64_t seed);

class PathFile {
 public:
  PathFile();
  ~PathFile();
  void Output(unsigned long data);
  void SavePath();
  void RemovePath();
  void Flush();
  void Clear();
  void SetOutDir(string address);
  void CreateDir();
 private:
  void CreateFile();
  string file_name_;
  ofstream fout;
  unsigned long* data_;
  unsigned long count_;
  unsigned short path_count_;
};

class PathMsg {
 public:
  PathMsg();
  PathMsg(unsigned long* data);
  PathMsg(unsigned long path_id, unsigned long start_edge);
  unsigned long* Data();
  unsigned long Data(unsigned long index);
  void Set(unsigned long path_id, unsigned long start_edge);
  void Set(unsigned long path_id);

 private:
  unsigned long data_[2];  // {path_id, start_edge}
};

class LocalMsg {
 public:
  LocalMsg();
  LocalMsg(const unsigned long path_id, const unsigned long start_edge, const unsigned long hash_value);
  void Link(LocalMsg* msg);
  unsigned long Query(const unsigned long path_id, const unsigned start_edge,
                      bool activate_fork, vector<unsigned long>& next_edges, const unsigned long hash_value);

 private:
  unsigned long start_edge_;
  unsigned long path_id_;
  unsigned long hash_value_;
  unsigned short count_;  // represent already return n next edges
  LocalMsg* next_;
};

class ParallelDFS {
 public:
  ParallelDFS();
  void run();
  void Server();
  void SetOutDir(string address);
  static void* ServerWrapper(void* object);

 private:
  void GetNextEdges(unsigned long edge, vector<unsigned long>& next_edges);

  /**
   * @brief: return single next edge
   */
  unsigned long SearchNextEdge(const unsigned long edge,
                               const unsigned long path_id,
                               const unsigned start_edge,
                               const unsigned long cur_fork,
                               unsigned long& next_fork,
                               unsigned long& hash_value);

  unsigned long QueryNextEdge(const unsigned long edge, PathMsg& path_msg,
                              unsigned long* fork_edges,
                              unsigned long& hash_value);

  /**
   * @brief:
   * @param start [out] The start edge in local
   * @param end [out]  The end edge in local
   * @return The number of edges in local
   */
  void GetLocalEdges(vector<unsigned long>& local_edges);

  /**
   * @brief: Check the edge is a global start edge (in degree is 0) or not.
   * @param  edge The edge need to be checked.
   * @return {*}
   */
  bool CheckGlobalStartEdge(unsigned long edge);
  int GetEdgeOwner(unsigned long edge);

  void Client();
  void StopServer();

  vector<unsigned long> local_edges_;
  unsigned long local_edges_size_;
  unsigned long local_path_count_;
  LocalMsg* local_msg_;
  vector<unsigned long>* next_edges_;
  string output_dir_;
};
