/*
 * @Brief: file description
 * @Author: liudy
 * @Email: deyin.liu@nscc-gz.cn
 * @Date: 2021-07-27 17:06:04
 * @LastEditors: liudy
 * @LastEditTime: 2021-07-30 16:00:02
 */
/*
 * @Brief: file description
 * @Author: liudy
 * @Email: deyin.liu@nscc-gz.cn
 * @Date: 2021-07-27 15:14:09
 * @LastEditors: liudy
 * @LastEditTime: 2021-07-28 09:28:40
 */
#include <iostream>
#include <vector>
#include <string>

#include "mpi.h"


using std::vector;
using std::string;


class ParallelDFS {
 public:
  void run();

 private:
  void GetNextEdge(unsigned long edge, vector<unsigned long>& next_edges);

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
  unsigned long GetRemoteLocalId(unsigned long edge);
  unsigned long GetRemoteLocalN(int target_process);


  bool CheckExitLocal(unsigned long edge);
  void LocalDFS(unsigned long start_edge, unsigned long cur_edge,
                vector<unsigned long>& local_paths, unsigned long path_len);
  void LocalPathSearch(vector<unsigned long>& local_paths);
  void GetPathFromRemote(int target_process, unsigned long target_edge,
                         vector<unsigned long>& remote_path, MPI_Win& win);

  /**
   * @brief: Collect the global paths  
   * @param  master_process The global paths is collected in master porcess.
   * @param global_paths The global paths search in each process. 
   * @return {*}
   */  
  void CollectGlobalPath(int master_process, vector<vector<unsigned long> >& global_paths);

  void LocalPathSeqSearch(vector<unsigned long>& query_msg,
                          unsigned long cur_edge, unsigned long cur_len,
                          vector<unsigned long>& path_seq);

  void SinglePathSeqSearch(int master_process,
                           vector<unsigned long>& global_path,
                           vector<unsigned long>& path_seq);

  /**
   * @brief: According to the global path to get the path sequence. 
   */
  void PathSeqSearch(int master_process,
                     vector<vector<unsigned long> >& global_paths,
                     vector<vector<unsigned long> >& global_path_seq);

  /**
   * @brief: Search in local paths to find paths that start from the local
   * start edge
   * @param local_paths [in]
   * @param  local_start_edge [in]
   * @param  local_end_edges [out] The local end edges of paths
   * @param  paths_len [out] The length of paths
   * @return
   */
  void QueryLocalPath(vector<unsigned long>& local_paths,
                      unsigned long local_start_edge,
                      vector<unsigned long>& local_end_edges,
                      vector<unsigned long>& paths_len);

  void GlobalDFS(MPI_Win& win, unsigned long cur_end_edge,
                 vector<unsigned long>& single_global_path,
                 vector<vector<unsigned long> >& global_paths);
  void GlobalPathSearch(MPI_Win& win,
                        vector<vector<unsigned long> >& global_paths,
                        vector<unsigned long>& local_paths);

  /**
   * @brief: Send a signal to other processes, tell them job is finish, stop to
   * wait for receiving message.
   * @param master_process
   */
  void SendStopSignal(int master_process);

  void SendQueryToRemote(int master_process,
                         vector<unsigned long>& global_path);

  void GetAnswerFromRemote(int master_process,
                           vector<unsigned long>& global_path,
                           vector<unsigned long>& path_seq);

  void WaitForQuerying(int master_process);
  void ResetLocalVisitedMap();
  void ResetGlobalVisitedMap();
  void SaveToFile(vector<vector<unsigned long> >&  path_seq, string file_name = "Path.csv");
  vector<unsigned long> local_edges_;
  bool local_path_seq_complete_;  // Used in LocalPathSeqSearch function to
                                 // indicate the lcoal path sequence search is
                                 // completed
  vector<bool> local_visited_map_;
  vector<bool> global_visited_map_;
};