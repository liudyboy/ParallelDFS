/*
 * @Brief: file description
 * @Author: liudy
 * @Email: deyin.liu@nscc-gz.cn
 * @Date: 2021-07-25 16:39:41
 * @LastEditors: liudy
 * @LastEditTime: 2021-07-29 09:54:15
 */
#include "parallel_dfs.h"
#include <algorithm>
#include <map>
#include <fstream>
#include "mpi.h"

using std::map;
using std::cout;
using std::endl;
using std::ofstream;
using std::ios;

// caution: actually the values is 18446744073709551615 (Maximum of unsigned long)
const unsigned long kGlobalEndEdge = -1;  
const unsigned long kStopSignal = -2;



// Test example 1

// void ParallelDFS::GetNextEdge(unsigned long edge, vector<unsigned long>& next_edges) {
//   // cout << "query edge =" << edge << endl;
//   map<unsigned long, vector<unsigned long> > graph;
//   vector<unsigned long> children0{1, 2};
//   vector<unsigned long> children1{10};
//   vector<unsigned long> children2{5, 6};
//   vector<unsigned long> children3{5, 6};
//   vector<unsigned long> children4{10};
//   vector<unsigned long> children5{7};
//   vector<unsigned long> children6{8};
//   vector<unsigned long> children7;
//   vector<unsigned long> children8{9};
//   vector<unsigned long> children9;
//   vector<unsigned long> children10{11};
//   vector<unsigned long> children11;
//   vector<unsigned long> children12{10};
//   vector<unsigned long> children13{5, 6};
//   graph[0] = children0;
//   graph[1] = children1;
//   graph[2] = children2;
//   graph[3] = children3;
//   graph[4] = children4;
//   graph[5] = children5;
//   graph[6] = children6;
//   graph[7] = children7;
//   graph[8] = children8;
//   graph[9] = children9;
//   graph[10] = children10;
//   graph[11] = children11;
//   graph[12] = children12;
//   graph[13] = children13;

//   next_edges = graph[edge];
// }

// void ParallelDFS::GetLocalEdges(vector<unsigned long>& local_edges) {
//   int my_id, process_num;
//   MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
//   MPI_Comm_size(MPI_COMM_WORLD, &process_num);
//   unsigned long start, end;
//   if (my_id == 0) {
//     start = 0;
//     end = 6;
//   }
//   if (my_id == 1) {
//     start = 6;
//     end = 14;
//   }
//   for (size_t i = start; i < end; i++) {
//     local_edges.push_back(i);
//   }
// }

// bool ParallelDFS::CheckGlobalStartEdge(unsigned long edge) {
//   if (edge == 0 || edge == 3 || edge == 4 || edge == 13 || edge == 12) return true;
//   return false;
// }

// int ParallelDFS::GetEdgeOwner(unsigned long edge) {
//   if (edge >= 0 && edge < 6)
//     return 0;
//   else
//     return 1;
// }

// unsigned long ParallelDFS::GetRemoteLocalId(unsigned long edge) {
//   if (edge > 5) return edge - 6;
//   return edge;
// }

// unsigned long ParallelDFS::GetRemoteLocalN(int target_process) {
//   if (target_process == 0) return 6;
//   return 8;

// }

// test example 2

void ParallelDFS::GetNextEdge(unsigned long edge, vector<unsigned long>& next_edges) {
  map<unsigned long, vector<unsigned long> > graph;
  vector<unsigned long> children0{10};
  vector<unsigned long> children1{4};
  vector<unsigned long> children2{5};
  vector<unsigned long> children3{7, 8};
  vector<unsigned long> children4{11};
  vector<unsigned long> children5{7, 8};
  vector<unsigned long> children6{9};
  vector<unsigned long> children7 {1};
  vector<unsigned long> children8{12};
  vector<unsigned long> children9{1};
  vector<unsigned long> children10{2};
  vector<unsigned long> children11;
  vector<unsigned long> children12{6};

  graph[0] = children0;
  graph[1] = children1;
  graph[2] = children2;
  graph[3] = children3;
  graph[4] = children4;
  graph[5] = children5;
  graph[6] = children6;
  graph[7] = children7;
  graph[8] = children8;
  graph[9] = children9;
  graph[10] = children10;
  graph[11] = children11;
  graph[12] = children12;

  next_edges = graph[edge];
}

void ParallelDFS::GetLocalEdges(vector<unsigned long>& local_edges) {
  int my_id, process_num;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Comm_size(MPI_COMM_WORLD, &process_num);
  unsigned long start, end;
  if (my_id == 0) {
    start = 0;
    end = 3;
  }
  if (my_id == 1) {
    start = 3;
    end = 8;
  }

  if (my_id == 2) {
    start = 8;
    end = 13;
  }
  for (size_t i = start; i < end; i++) {
    local_edges.push_back(i);
  }
}

bool ParallelDFS::CheckGlobalStartEdge(unsigned long edge) {
  if (edge == 0 || edge == 3) return true;
  return false;
}

int ParallelDFS::GetEdgeOwner(unsigned long edge) {
  if (edge >= 0 && edge < 3)
    return 0;
  else if (edge >= 3 && edge < 8)
    return 1;
  else
    return 2;
}

unsigned long ParallelDFS::GetRemoteLocalId(unsigned long edge) {
  if (edge >= 0 && edge < 3)
    return edge;
  else if (edge >= 3 && edge < 8)
    return edge - 3;
  else
    return edge - 8;
}

unsigned long ParallelDFS::GetRemoteLocalN(int target_process) {
  if (target_process == 0) return 3;
  if (target_process == 1) return 5;
  return 5;
}

bool ParallelDFS::CheckExitLocal(unsigned long edge) {
  if (local_edges_.size() != 0) {
    auto it = find(local_edges_.begin(), local_edges_.end(), edge);
    if (it != local_edges_.end()) return true;
  }
  return false;
}

void ParallelDFS::LocalDFS(unsigned long start_edge,
                           unsigned long cur_edge,
                           vector<unsigned long>& local_paths,
                           unsigned long path_len) {
  vector<unsigned long> next_edges;
  GetNextEdge(cur_edge, next_edges);
  unsigned long next_edges_size = next_edges.size();
  path_len++;
  if (next_edges_size == 0) {
    local_paths.push_back(kGlobalEndEdge);
    local_paths.push_back(path_len);
    local_paths[start_edge - local_edges_[0]]++;
    return;
  }
  for (size_t i = 0; i < next_edges_size; i++) {
    if (CheckExitLocal(next_edges[i]))
      LocalDFS(start_edge, next_edges[i], local_paths, path_len);
    else {
      local_paths.push_back(next_edges[i]);
      local_paths.push_back(path_len);
      local_paths[start_edge - local_edges_[0]]++;
    }
  }
}

void ParallelDFS::LocalPathSearch(vector<unsigned long>& local_paths) {
  GetLocalEdges(local_edges_);
  unsigned long local_n = local_edges_.size(); 
  for (size_t i = 0; i < local_n; i++) local_paths.push_back(0);
  for (size_t i = 0; i < local_n; i++) {
    LocalDFS(local_edges_[i], local_edges_[i], local_paths, 0);
  }
}

void ParallelDFS::QueryLocalPath(vector<unsigned long>& local_paths,
                     unsigned long local_start_edge,
                     vector<unsigned long>& local_end_edges,
                     vector<unsigned long>& paths_len) {
  unsigned long local_n = local_edges_.size();
  unsigned long num_paths = local_paths[local_start_edge - local_edges_[0]];
  unsigned long start_index = 0;
  for (size_t i = 0; i < local_start_edge - local_edges_[0]; i++) {
    start_index += local_paths[i];
  }
  start_index = 2 * start_index + local_n;
  for (size_t i = 0; i < num_paths; i++) {
    local_end_edges.push_back(local_paths[start_index]);
    paths_len.push_back(local_paths[start_index + 1]);
    start_index += 2;
  }
}

void ParallelDFS::GetPathFromRemote(int target_process, unsigned long target_edge,
                       vector<unsigned long>& remote_path, MPI_Win& win) {
  MPI_Win_lock(MPI_LOCK_SHARED, target_process, 0, win);
  unsigned long remote_local_id = GetRemoteLocalId(target_edge);
  unsigned long* num_buf = new unsigned long[remote_local_id + 1];
  MPI_Get(num_buf, remote_local_id + 1, MPI_UNSIGNED_LONG, target_process, 0, remote_local_id + 1, MPI_UNSIGNED_LONG, win);
  MPI_Win_flush(target_process, win);
  unsigned long count = 0;
  for (size_t i = 0; i < remote_local_id; i++) {
    count += num_buf[i];
  } 
  unsigned long remote_local_n = GetRemoteLocalN(target_process); 
  unsigned long offset = count * 2 + remote_local_n;
  remote_path.resize(2 * num_buf[remote_local_id]);
  MPI_Get(remote_path.data(), 2 * num_buf[remote_local_id], MPI_UNSIGNED_LONG,
          target_process, offset, 2 * num_buf[remote_local_id],
          MPI_UNSIGNED_LONG, win);
  MPI_Win_unlock(target_process, win);
}

void ParallelDFS::GlobalDFS(MPI_Win& win, unsigned long cur_end_edge,
                            vector<unsigned long>& single_global_path,
                            vector<vector<unsigned long> >& global_paths) {
  if (cur_end_edge == kGlobalEndEdge) {
    global_paths.push_back(single_global_path);
    return;
  }
  int target_process = GetEdgeOwner(cur_end_edge);
  vector<unsigned long> remote_path;
  GetPathFromRemote(target_process, cur_end_edge, remote_path, win);

  for (size_t i = 0; i < remote_path.size(); i+=2) {
    unsigned long cur_len = single_global_path.back(); 
    single_global_path.push_back(remote_path[i]);
    single_global_path.push_back(remote_path[i+1] + cur_len);
    GlobalDFS(win, remote_path[i], single_global_path, global_paths);
    single_global_path.pop_back();
    single_global_path.pop_back();
  }
}

void ParallelDFS::GlobalPathSearch(MPI_Win& win,
                                   vector<vector<unsigned long> >& global_paths,
                                   vector<unsigned long>& local_paths) {
  for (size_t i = 0; i < local_edges_.size(); i++) {
    unsigned long edge = local_edges_[i];
    vector<unsigned long> single_global_path;
    if (CheckGlobalStartEdge(edge)) {
      single_global_path.push_back(edge);
      vector<unsigned long> local_end_edges;
      vector<unsigned long> paths_len;
      QueryLocalPath(local_paths, edge, local_end_edges, paths_len);

      for (size_t i = 0; i < local_end_edges.size(); i++) {
        single_global_path.push_back(local_end_edges[i]);
        single_global_path.push_back(paths_len[i]);
        GlobalDFS(win, local_end_edges[i], single_global_path, global_paths);
        single_global_path.pop_back();
        single_global_path.pop_back();
      }
    }
  }
}

void ParallelDFS::CollectGlobalPath(int master_process, vector<vector<unsigned long> >& global_paths) {
  int my_id, num_process;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Comm_size(MPI_COMM_WORLD, &num_process);
  unsigned long num_local_paths = global_paths.size();
  vector<unsigned long> num_paths;
  if (my_id == master_process) num_paths.resize(num_process);
  MPI_Gather(&num_local_paths, 1, MPI_UNSIGNED_LONG, num_paths.data(), 1,
             MPI_UNSIGNED_LONG, master_process, MPI_COMM_WORLD);

  if (my_id == master_process) {
    MPI_Status status;
    int msg_size;
    for (int send_process = 0; send_process < num_process; send_process++) {
      vector<unsigned long> single_path;
      if (send_process != master_process) {
        for (size_t i = 0; i < num_paths[send_process]; i++) {
          MPI_Probe(send_process, send_process, MPI_COMM_WORLD, &status);
          MPI_Get_count(&status, MPI_UNSIGNED_LONG, &msg_size);
          single_path.resize(msg_size);
          MPI_Recv(single_path.data(), msg_size, MPI_UNSIGNED_LONG,
                   send_process, send_process, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
          global_paths.push_back(single_path);
        }
      }
    }
  } else {
    for (size_t i = 0; i < global_paths.size(); i++) {
      MPI_Send(global_paths[i].data(), global_paths[i].size(),
               MPI_UNSIGNED_LONG, master_process, my_id, MPI_COMM_WORLD);
    }
  }
}

void ParallelDFS::LocalPathSeqSearch(vector<unsigned long>& query_msg,
                                     unsigned long cur_edge,
                                     unsigned long cur_len,
                                     vector<unsigned long>& path_seq) {
  if ((cur_len == query_msg[2] && query_msg[1] == cur_edge) || local_path_seq_complete_) {
    path_seq.pop_back();
    local_path_seq_complete_ = true;
    return;
  }
  if (cur_edge == kGlobalEndEdge) return;
  vector<unsigned long> next_edges;
  GetNextEdge(cur_edge, next_edges);
  for (size_t i = 0; i < next_edges.size(); i++) {
    path_seq.push_back(next_edges[i]);
    LocalPathSeqSearch(query_msg, next_edges[i], cur_len + 1, path_seq);
    if (local_path_seq_complete_) break;
    path_seq.pop_back();
  }
  if (next_edges.size() == 0) {
    path_seq.push_back(kGlobalEndEdge);
    LocalPathSeqSearch(query_msg, kGlobalEndEdge, cur_len + 1, path_seq);
    if (!local_path_seq_complete_) path_seq.pop_back();
  }
}

void ParallelDFS::SendStopSignal(int master_process) {
  int num_process;
  MPI_Comm_size(MPI_COMM_WORLD, &num_process);
    for (int recv_processs = 0; recv_processs < num_process; recv_processs++) {
      if (recv_processs != master_process) {
        MPI_Request request;
        vector<unsigned long> query_msg{kStopSignal, 0, 0};
        MPI_Isend(query_msg.data(), 3, MPI_UNSIGNED_LONG, recv_processs,
                  recv_processs, MPI_COMM_WORLD, &request);
      }
    }
}

void ParallelDFS::SendQueryToRemote(int master_process, vector<unsigned long>& global_path) {
    unsigned long start_edge = global_path[0];
    unsigned long pre_len = 0;
    for (size_t i = 1; i < global_path.size(); i += 2) {
      MPI_Request request;
      vector<unsigned long> query_msg{start_edge, global_path[i], global_path[i+1] - pre_len};
      int remote_process = GetEdgeOwner(start_edge);
      if (remote_process != master_process) {
        MPI_Isend(query_msg.data(), 3, MPI_UNSIGNED_LONG, remote_process,
                  remote_process, MPI_COMM_WORLD, &request);
      }
      start_edge = global_path[i];
      pre_len = global_path[i + 1];
    }
}

void ParallelDFS::GetAnswerFromRemote(int master_process, vector<unsigned long>& global_path, vector<unsigned long>& path_seq) {
    unsigned long start_edge = global_path[0];
    unsigned long pre_len = 0;
    for (size_t i = 1; i < global_path.size(); i += 2) {
      int remote_process = GetEdgeOwner(start_edge);
      vector<unsigned long> path_seq_segment;
      MPI_Status status;
      int seq_size;
      if (remote_process != master_process) {
          MPI_Probe(remote_process, remote_process, MPI_COMM_WORLD, &status);
          MPI_Get_count(&status, MPI_UNSIGNED_LONG, &seq_size);
          path_seq_segment.resize(seq_size);
          MPI_Recv(path_seq_segment.data(), seq_size, MPI_UNSIGNED_LONG,
                   remote_process, remote_process, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
      } else {
        vector<unsigned long> query_msg{start_edge, global_path[i],
                                        global_path[i + 1] - pre_len};
        local_path_seq_complete_ = false;
        path_seq_segment.push_back(query_msg[0]);
        LocalPathSeqSearch(query_msg, query_msg[0], 0, path_seq_segment);
      }

      for (size_t j = 0; j < path_seq_segment.size(); j++) {
        path_seq.push_back(path_seq_segment[j]);
      }
      start_edge = global_path[i];
      pre_len = global_path[i + 1];
    }
}

void ParallelDFS::WaitForQuerying(int master_process) {
  int my_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  vector<unsigned long> query_msg;
  vector<unsigned long> path_seq;
  query_msg.resize(3);
  while (true) {
    MPI_Recv(query_msg.data(), 3, MPI_UNSIGNED_LONG, master_process, my_id,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (query_msg[0] != kStopSignal) {
      MPI_Request request;
      path_seq.clear();

      local_path_seq_complete_ = false;
      path_seq.push_back(query_msg[0]);
      LocalPathSeqSearch(query_msg, query_msg[0], 0, path_seq);

      MPI_Isend(path_seq.data(), path_seq.size(), MPI_UNSIGNED_LONG,
                master_process, my_id, MPI_COMM_WORLD, &request);
    } else {
      break;
    }
  }
}

void ParallelDFS::SinglePathSeqSearch(int master_process, vector<unsigned long>& global_path, vector<unsigned long>& path_seq) {
  int my_id, num_process;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Comm_size(MPI_COMM_WORLD, &num_process);

  if (my_id == master_process) {
    SendQueryToRemote(master_process, global_path);
    GetAnswerFromRemote(master_process, global_path, path_seq);
    SendStopSignal(master_process);
  } else {
    WaitForQuerying(master_process);
  } 
  MPI_Barrier(MPI_COMM_WORLD);
}

void ParallelDFS::SaveToFile(vector<vector<unsigned long> >& path_seq, string file_name) {
  ofstream fout;
  fout.open(file_name, ios::out | ios::app); 
  for (size_t i = 0; i < path_seq.size(); i++) {
    fout << path_seq[i][0];
    for (size_t j = 1; j < path_seq[i].size(); j++) {
      fout << ", " << path_seq[i][j]; 
    }
    fout << endl;
  }
  fout.close();
}
void ParallelDFS::PathSeqSearch(
    int master_process, vector<vector<unsigned long> >& global_paths,
    vector<vector<unsigned long> >& global_path_seq) {
  int my_id, process_num;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Comm_size(MPI_COMM_WORLD, &process_num);
  int num_paths = 0;
  if (my_id == master_process) num_paths = global_paths.size();
  MPI_Bcast(&num_paths, 1, MPI_INT, master_process, MPI_COMM_WORLD);
  for (size_t i = 0; i < num_paths; i++) {
    vector<unsigned long> single_path_seq;
    SinglePathSeqSearch(master_process, global_paths[i], single_path_seq);
    global_path_seq.push_back(single_path_seq);
    MPI_Barrier(MPI_COMM_WORLD);
  }
}
void ParallelDFS::run() {
  int my_id, process_num;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Comm_size(MPI_COMM_WORLD, &process_num);
  int master_process = 0;

  vector<unsigned long> local_paths;
  vector<vector<unsigned long> > global_paths;
  LocalPathSearch(local_paths);

  MPI_Win win;
  MPI_Win_create(local_paths.data(), sizeof(unsigned long) * local_paths.size(),
                 sizeof(unsigned long), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
  GlobalPathSearch(win, global_paths, local_paths);
  CollectGlobalPath(master_process, global_paths);

  vector<vector<unsigned long> > global_path_seq;
  PathSeqSearch(master_process, global_paths, global_path_seq);
  if (my_id == master_process) SaveToFile(global_path_seq);

// print the paths to screen
  if (my_id == master_process) {
    cout << "==========================" << global_paths.size()
         << " path======================" << endl;
    for (size_t i = 0; i < global_path_seq.size(); i++) {
      for (size_t j = 0; j < global_path_seq[i].size(); j++) {
        cout << global_path_seq[i][j] << " ";
      }
      cout << endl << "------------------------------------" << endl;
      ;
    }
  }
  MPI_Win_free(&win);  // free window
}