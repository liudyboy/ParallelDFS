/*
 * @Brief: file description
 * @Author: liudy
 * @Email: deyin.liu@nscc-gz.cn
 * @Date: 2021-07-25 16:39:41
 * @LastEditors: liudy
 * @LastEditTime: 2021-08-05 10:27:03
 */
#include "parallel_dfs.h"
#include <algorithm>
#include <map>
#include <fstream>
#include <pthread.h>
#include "mpi.h"


using std::map;
using std::cout;
using std::endl;
using std::ofstream;
using std::ios;
using std::to_string;

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

// void ParallelDFS::GetNextEdge(unsigned long edge, vector<unsigned long>& next_edges) {
//   map<unsigned long, vector<unsigned long> > graph;
//   vector<unsigned long> children0{10};
//   vector<unsigned long> children1{4};
//   vector<unsigned long> children2{5};
//   vector<unsigned long> children3{7, 8};
//   vector<unsigned long> children4{11};
//   vector<unsigned long> children5{7, 8};
//   vector<unsigned long> children6{9};
//   vector<unsigned long> children7 {1};
//   vector<unsigned long> children8{12};
//   vector<unsigned long> children9{1};
//   vector<unsigned long> children10{2};
//   vector<unsigned long> children11;
//   vector<unsigned long> children12{6};

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

//   next_edges = graph[edge];
// }

// void ParallelDFS::GetLocalEdges(vector<unsigned long>& local_edges) {
//   int my_id, process_num;
//   MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
//   MPI_Comm_size(MPI_COMM_WORLD, &process_num);
//   unsigned long start, end;
//   if (my_id == 0) {
//     start = 0;
//     end = 3;
//   }
//   if (my_id == 1) {
//     start = 3;
//     end = 8;
//   }

//   if (my_id == 2) {
//     start = 8;
//     end = 13;
//   }
//   for (size_t i = start; i < end; i++) {
//     local_edges.push_back(i);
//   }
// }

// bool ParallelDFS::CheckGlobalStartEdge(unsigned long edge) {
//   if (edge == 0 || edge == 3) return true;
//   return false;
// }

// int ParallelDFS::GetEdgeOwner(unsigned long edge) {
//   if (edge >= 0 && edge < 3)
//     return 0;
//   else if (edge >= 3 && edge < 8)
//     return 1;
//   else
//     return 2;
// }

// unsigned long ParallelDFS::GetRemoteLocalId(unsigned long edge) {
//   if (edge >= 0 && edge < 3)
//     return edge;
//   else if (edge >= 3 && edge < 8)
//     return edge - 3;
//   else
//     return edge - 8;
// }

// unsigned long ParallelDFS::GetRemoteLocalN(int target_process) {
//   if (target_process == 0) return 3;
//   if (target_process == 1) return 5;
//   return 5;
// }

// example 3  (circle exit among processes)
// void ParallelDFS::GetNextEdge(unsigned long edge, vector<unsigned long>& next_edges) {
//   map<unsigned long, vector<unsigned long> > graph;
//   vector<unsigned long> children0{6, 7};
//   vector<unsigned long> children1{4, 9};
//   vector<unsigned long> children2{1};
//   vector<unsigned long> children3{1};
//   vector<unsigned long> children4;
//   vector<unsigned long> children5;
//   vector<unsigned long> children6{1};
//   vector<unsigned long> children7 {3, 5};
//   vector<unsigned long> children8{3, 5};
//   vector<unsigned long> children9{3, 5};

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

//   next_edges = graph[edge];
// }

// void ParallelDFS::GetLocalEdges(vector<unsigned long>& local_edges) {
//   int my_id, process_num;
//   MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
//   MPI_Comm_size(MPI_COMM_WORLD, &process_num);
//   unsigned long start, end;
//   if (my_id == 0) {
//     start = 0;
//     end = 2;
//   }
//   if (my_id == 1) {
//     start = 3;
//     end = 5;
//   }

//   if (my_id == 2) {
//     start = 6;
//     end = 9;
//   }
//   local_edges.push_back(start);
//   local_edges.push_back(end);
// }

// bool ParallelDFS::CheckGlobalStartEdge(unsigned long edge) {
//   if (edge == 0 || edge == 2 || edge == 8) return true;
//   return false;
// }

// int ParallelDFS::GetEdgeOwner(unsigned long edge) {
//   if (edge >= 0 && edge < 3)
//     return 0;
//   else if (edge >= 3 && edge < 6)
//     return 1;
//   else
//     return 2;
// }

// unsigned long ParallelDFS::GetRemoteLocalId(unsigned long edge) {
//   if (edge >= 0 && edge < 3)
//     return edge;
//   else if (edge >= 3 && edge < 6)
//     return edge - 3;
//   else
//     return edge - 6;
// }

// unsigned long ParallelDFS::GetRemoteLocalN(int target_process) {
//   if (target_process == 0) return 3;
//   if (target_process == 1) return 3;
//   return 4;
// }


// example 4 (circle exit in one porcess)
void ParallelDFS::GetNextEdge(unsigned long edge, vector<unsigned long>& next_edges) {
  map<unsigned long, vector<unsigned long> > graph;
  vector<unsigned long> children0{1};
  vector<unsigned long> children1{2, 4};
  vector<unsigned long> children2{3, 9};
  vector<unsigned long> children3{1};
  vector<unsigned long> children4;
  vector<unsigned long> children5{1};
  vector<unsigned long> children6{5, 7};
  vector<unsigned long> children7 {3, 9};
  vector<unsigned long> children8{3, 9};
  vector<unsigned long> children9;

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

  next_edges = graph[edge];
}

void ParallelDFS::GetLocalEdges(vector<unsigned long>& local_edges) {
  int my_id, process_num;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Comm_size(MPI_COMM_WORLD, &process_num);
  unsigned long start, end;
  if (my_id == 0) {
    start = 0;
    end = 4;
  }
  if (my_id == 1) {
    start = 5;
    end = 9;
  }
  local_edges.push_back(start);
  local_edges.push_back(end);
}

bool ParallelDFS::CheckGlobalStartEdge(unsigned long edge) {
  if (edge == 0 || edge == 6 || edge == 8) return true;
  return false;
}

int ParallelDFS::GetEdgeOwner(unsigned long edge) {
  if (edge >= 0 && edge < 5)
    return 0;
  return 1;
}

unsigned long ParallelDFS::GetRemoteLocalId(unsigned long edge) {
  if (edge >= 0 && edge < 5)
    return edge;
  return edge - 5;
}

unsigned long ParallelDFS::GetRemoteLocalN(int target_process) {
  return 5;
}



bool ParallelDFS::CheckExitLocal(unsigned long edge) {
  if (local_edges_.size() != 0) {
    if (edge >= local_edges_[0] && edge <= local_edges_[1]) return true;
  }
  return false;
}

void ParallelDFS::ResetLocalVisitedMap() {
  if (local_visited_map_.size() == 0) {
    local_visited_map_.resize(local_edges_[1] - local_edges_[0] + 1);
  }
  for (size_t i = 0; i < local_visited_map_.size(); i++) {
    local_visited_map_[i] = false;
  }
}

void ParallelDFS::ResetGlobalVisitedMap() {
  if (global_visited_map_.size() == 0) {
    unsigned long global_size = 0;
    MPI_Allreduce(&local_edges_size_, &global_size, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    global_visited_map_.resize(global_size);
  }

  for (size_t i = 0; i < global_visited_map_.size(); i++) {
    global_visited_map_[i] = false;
  }
}

void ParallelDFS::LocalDFS(unsigned long start_edge, unsigned long cur_edge,
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
    if (CheckExitLocal(next_edges[i])) {
      if (local_visited_map_[next_edges[i] - local_edges_[0]])
        continue;
      local_visited_map_[next_edges[i] - local_edges_[0]] = true;
      LocalDFS(start_edge, next_edges[i], local_paths, path_len);
      local_visited_map_[next_edges[i] - local_edges_[0]] = false;
    }
    else {
      local_paths.push_back(next_edges[i]);
      local_paths.push_back(path_len);
      local_paths[start_edge - local_edges_[0]]++;
    }
  }
}

void ParallelDFS::LocalPathSearch(vector<unsigned long>& local_paths) {
  GetLocalEdges(local_edges_);
  local_edges_size_ = local_edges_[1] - local_edges_[0] + 1;
  for (size_t i = 0; i < local_edges_size_; i++) local_paths.push_back(0);
  for (size_t i = 0; i < local_edges_size_; i++) {
    unsigned long edge = local_edges_[0] + i;
    ResetLocalVisitedMap();
    local_visited_map_[i] = true;
    LocalDFS(edge, edge, local_paths, 0);
  }
}

void ParallelDFS::QueryLocalPath(vector<unsigned long>& local_paths,
                     unsigned long local_start_edge,
                     vector<unsigned long>& local_end_edges,
                     vector<unsigned long>& paths_len) {
  unsigned long local_n = local_edges_size_;
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
  delete []num_buf;
}

void ParallelDFS::GlobalDFS(MPI_Win& win, unsigned long cur_end_edge,
                            vector<unsigned long>& single_global_path,
                            vector<vector<unsigned long> >& global_paths) {

  if (single_global_path.back() > max_path_len_) // 防止路径长度异常一直在增加,确保能退出
    return;
  if (cur_end_edge == kGlobalEndEdge) {
    global_paths.push_back(single_global_path);
    return;
  }
  if (global_visited_map_[cur_end_edge]) return;
  global_visited_map_[cur_end_edge] = true;

  int target_process = GetEdgeOwner(cur_end_edge);

  vector<unsigned long> remote_path;
  GetPathFromRemote(target_process, cur_end_edge, remote_path, win);

  for (size_t i = 0; i < remote_path.size(); i += 2) {
    unsigned long cur_len = single_global_path.back();
    single_global_path.push_back(remote_path[i]);
    single_global_path.push_back(remote_path[i + 1] + cur_len);
    GlobalDFS(win, remote_path[i], single_global_path, global_paths);
    single_global_path.pop_back();
    single_global_path.pop_back();
  }
  global_visited_map_[cur_end_edge] = false;
}

void ParallelDFS::GlobalPathSearch(MPI_Win& win,
                                   vector<vector<unsigned long> >& global_paths,
                                   vector<unsigned long>& local_paths) {
  for (size_t i = 0; i < local_edges_size_; i++) {
    unsigned long edge = i + local_edges_[0];
    vector<unsigned long> single_global_path;
    if (CheckGlobalStartEdge(edge)) {
      ResetGlobalVisitedMap();
      single_global_path.push_back(edge);
      vector<unsigned long> local_end_edges;
      vector<unsigned long> paths_len;
      QueryLocalPath(local_paths, edge, local_end_edges, paths_len);

      for (size_t j = 0; j < local_end_edges.size(); j++) {
        single_global_path.push_back(local_end_edges[j]);
        single_global_path.push_back(paths_len[j]);
        GlobalDFS(win, local_end_edges[j], single_global_path, global_paths);
        single_global_path.pop_back();
        single_global_path.pop_back();
      }
    }
  }
}


void ParallelDFS::LocalPathSeqSearch(vector<unsigned long>& query_msg,
                                     unsigned long cur_edge,
                                     unsigned long cur_len,
                                     vector<unsigned long>& path_seq) {

  // int my_id, num_process;
  // MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  // MPI_Comm_size(MPI_COMM_WORLD, &num_process);
  // cout << "my_id = " << my_id  << "cur_edge = " << cur_edge << endl;
  if ((cur_len == query_msg[2] && query_msg[1] == cur_edge) || local_path_seq_complete_) {
    path_seq.pop_back();
    local_path_seq_complete_ = true;
    return;
  }
  if (!CheckExitLocal(cur_edge)) {
    return;
  }
  if (cur_edge == kGlobalEndEdge) return;
  vector<unsigned long> next_edges;
  GetNextEdge(cur_edge, next_edges);
  for (size_t i = 0; i < next_edges.size(); i++) {
    if (CheckExitLocal(next_edges[i])) {
      if (local_visited_map_[next_edges[i] - local_edges_[0]]) continue;
      local_visited_map_[next_edges[i] - local_edges_[0]] = true;
    }

    path_seq.push_back(next_edges[i]);
    LocalPathSeqSearch(query_msg, next_edges[i], cur_len + 1, path_seq);

    if (CheckExitLocal(next_edges[i]))
      local_visited_map_[next_edges[i] - local_edges_[0]] = false;

    if (local_path_seq_complete_) break;
    path_seq.pop_back();
  }
  if (next_edges.size() == 0) {
    path_seq.push_back(kGlobalEndEdge);
    LocalPathSeqSearch(query_msg, kGlobalEndEdge, cur_len + 1, path_seq);
    if (!local_path_seq_complete_) path_seq.pop_back();
  }
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

void ParallelDFS::SinglePathSeqSearch(unsigned long* global_path, unsigned long length, unsigned long path_id) {
  int my_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  unsigned long start_edge = global_path[0];
  unsigned long pre_len = 0;
  unsigned long query_times = 0;
  for (size_t i = 1; i < length; i += 2) {
    int remote_process = GetEdgeOwner(start_edge);
    vector<unsigned long> path_seq_segment;
    if (remote_process == my_id) {
      vector<unsigned long> query_msg{start_edge, global_path[i],
                                      global_path[i + 1] - pre_len};
      local_path_seq_complete_ = false;
      path_seq_segment.push_back(query_msg[0]);
      ResetLocalVisitedMap();
      local_visited_map_[query_msg[0] - local_edges_[0]] = true;
      LocalPathSeqSearch(query_msg, query_msg[0], 0, path_seq_segment);
      path_output_.ParallelWrite(path_seq_segment, pre_len, path_id);
    }
    start_edge = global_path[i];
    pre_len = global_path[i + 1];
  }
}



void ParallelDFS::PathSeqSearch(int master_process,
                                vector<vector<unsigned long> >& global_paths) {
  int my_id, process_num;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Comm_size(MPI_COMM_WORLD, &process_num);
  unsigned long num_paths = 0;
  if (my_id == master_process) num_paths = global_paths.size();
  MPI_Bcast(&num_paths, 1, MPI_UNSIGNED_LONG, master_process, MPI_COMM_WORLD);

  unsigned long num_files = 0;
  unsigned long lengths[NUM_FILE];

  for (size_t i = 0; i < num_paths; i += NUM_FILE) {
    num_files = num_paths - i > NUM_FILE ?  NUM_FILE : num_paths - i;

    if (my_id == master_process) {
      for (size_t j = 0; j < num_files; j++) {
        lengths[j] = global_paths[i + j].size();
      }
    }
    MPI_Bcast(lengths, num_files, MPI_UNSIGNED_LONG, master_process, MPI_COMM_WORLD);

    unsigned long total_size = 0;
    for (size_t j = 0; j < num_files; j++)
      total_size += lengths[j];
    unsigned long* batch_path = new unsigned long[total_size];

    if (my_id == master_process) {
      unsigned long index = 0;
      for (size_t p = 0; p < num_files; p++) {
        for (size_t q = 0; q < global_paths[i + p].size(); q++) {
          batch_path[index] = global_paths[i + p][q];
          index++;
        }
      }
    }
    MPI_Bcast(batch_path, total_size, MPI_UNSIGNED_LONG, master_process, MPI_COMM_WORLD);

    path_output_.CreateFiles(path_id_, num_files);
    unsigned long offset = 0;
    for (size_t j = 0; j < num_files; j++) {
      SinglePathSeqSearch(batch_path + offset, lengths[j], path_id_ + j);
      offset += lengths[j];
    }
    path_output_.CloseFiles();
    path_id_ += num_files;
    delete[] batch_path;
  }
}

void ParallelDFS::SetMaxLength(unsigned long max_len) {
  max_path_len_ = max_len;
}

void ParallelDFS::run() {
  // initialize variables
  path_id_ = 0;
  top_k_ = TOPK;

  int my_id, process_num;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Comm_size(MPI_COMM_WORLD, &process_num);
  int master_process = 0;

  vector<unsigned long> local_paths;
  vector<vector<unsigned long> > global_paths;
  LocalPathSearch(local_paths);

  ResetGlobalVisitedMap();

  MPI_Win win;
  MPI_Win_create(local_paths.data(), sizeof(unsigned long) * local_paths.size(),
                 sizeof(unsigned long), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
  GlobalPathSearch(win, global_paths, local_paths);
  MPI_Win_free(&win);  // free window

  vector<unsigned long>().swap(local_paths); // free memory

  for (int token = 0; token < process_num; token++) {
    PathSeqSearch(token, global_paths);
  }


  // print path result
  for (size_t i = 0; i < path_id_; i++) {
    path_output_.OpenFile(i);
    if (my_id == master_process) path_output_.SequenceRead(i);
    path_output_.CloseFile();
  }
}



ParallelStream::ParallelStream(string file_name_prefix) {
  file_handle_ = nullptr;
  file_name_prefix_ = file_name_prefix;
}
ParallelStream::~ParallelStream() {
}

void ParallelStream::CloseFiles() {
  for (size_t i = 0; i < num_files_; i++) {
    MPI_File_close(&file_handle_[i]);
  }
  delete[] file_handle_;
  file_handle_ = nullptr;
}

void ParallelStream::CloseFile() {
  MPI_File_close(&fh_);
}
void ParallelStream::ParallelWrite(vector<unsigned long>& path_seq,
                                   unsigned long offset, unsigned long path_id) {
  MPI_Status status;
  MPI_File_write_at(file_handle_[path_id - base_path_id_],
                    offset * sizeof(unsigned long), path_seq.data(),
                    path_seq.size(), MPI_UNSIGNED_LONG, &status);
}

void ParallelStream::SequenceRead(unsigned long path_id) {
  MPI_Status status;
  long long file_size;
  MPI_File_get_size(fh_, &file_size);
  long long size = file_size / sizeof(unsigned long);
  unsigned long edge;
  long long offset = 0;
  cout << "==================== path "<< path_id << ": size = " << size << "==================================" << endl;
  while(size--) {
    MPI_File_read_at(fh_, offset, &edge, 1, MPI_UNSIGNED_LONG, &status);
    offset += sizeof(unsigned long);
    cout << edge << " ";
  }
  cout << endl;
}

void ParallelStream::CreateFiles(unsigned long base_path_id, unsigned long num_files) {
  num_files_ = num_files;
  if (file_handle_ != nullptr)
    delete[] file_handle_;
  file_handle_ = new MPI_File[num_files];
  string file_name;
  base_path_id_ = base_path_id;
  for (size_t i = 0; i < num_files; i++) {
    file_name = FILE_ADDRESS + file_name_prefix_ + to_string(base_path_id_ + i);
    MPI_File_open(MPI_COMM_WORLD, file_name.c_str(),
                  MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &file_handle_[i]);
  }
}
void ParallelStream::OpenFile(unsigned long path_id) {
  string file_name;
  file_name = FILE_ADDRESS + file_name_prefix_ + to_string(path_id);
  MPI_File_open(MPI_COMM_WORLD, file_name.c_str(), MPI_MODE_RDWR, MPI_INFO_NULL,
                &fh_);
}
