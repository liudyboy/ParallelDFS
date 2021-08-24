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
void ParallelDFS::GetNextEdge(unsigned long edge, vector<unsigned long>& next_edges) {
  map<unsigned long, vector<unsigned long> > graph;
  vector<unsigned long> children0{6, 7};
  vector<unsigned long> children1{4, 9};
  vector<unsigned long> children2{1};
  vector<unsigned long> children3{1};
  vector<unsigned long> children4;
  vector<unsigned long> children5;
  vector<unsigned long> children6{1};
  vector<unsigned long> children7 {3, 5};
  vector<unsigned long> children8{3, 5};
  vector<unsigned long> children9{3, 5};

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
    end = 3;
  }
  if (my_id == 1) {
    start = 3;
    end = 6;
  }

  if (my_id == 2) {
    start = 6;
    end = 10;
  }
  for (size_t i = start; i < end; i++) {
    local_edges.push_back(i);
  }
}

bool ParallelDFS::CheckGlobalStartEdge(unsigned long edge) {
  if (edge == 0 || edge == 2 || edge == 8) return true;
  return false;
}

int ParallelDFS::GetEdgeOwner(unsigned long edge) {
  if (edge >= 0 && edge < 3)
    return 0;
  else if (edge >= 3 && edge < 6)
    return 1;
  else
    return 2;
}

unsigned long ParallelDFS::GetRemoteLocalId(unsigned long edge) {
  if (edge >= 0 && edge < 3)
    return edge;
  else if (edge >= 3 && edge < 6)
    return edge - 3;
  else
    return edge - 6;
}

unsigned long ParallelDFS::GetRemoteLocalN(int target_process) {
  if (target_process == 0) return 3;
  if (target_process == 1) return 3;
  return 4;
}


// example 4 (circle exit in one porcess)
// void ParallelDFS::GetNextEdge(unsigned long edge, vector<unsigned long>& next_edges) {
//   map<unsigned long, vector<unsigned long> > graph;
//   vector<unsigned long> children0{1};
//   vector<unsigned long> children1{2, 4};
//   vector<unsigned long> children2{3, 9};
//   vector<unsigned long> children3{1};
//   vector<unsigned long> children4;
//   vector<unsigned long> children5{1};
//   vector<unsigned long> children6{5, 7};
//   vector<unsigned long> children7 {3, 9};
//   vector<unsigned long> children8{3, 9};
//   vector<unsigned long> children9;

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
//     end = 5;
//   }
//   if (my_id == 1) {
//     start = 5;
//     end = 10;
//   }

//   for (size_t i = start; i < end; i++) {
//     local_edges.push_back(i);
//   }
// }

// bool ParallelDFS::CheckGlobalStartEdge(unsigned long edge) {
//   if (edge == 0 || edge == 6 || edge == 8) return true;
//   return false;
// }

// int ParallelDFS::GetEdgeOwner(unsigned long edge) {
//   if (edge >= 0 && edge < 5)
//     return 0;
//   return 1;
// }

// unsigned long ParallelDFS::GetRemoteLocalId(unsigned long edge) {
//   if (edge >= 0 && edge < 5)
//     return edge;
//   return edge - 5;
// }

// unsigned long ParallelDFS::GetRemoteLocalN(int target_process) {
//   return 5;
// }



bool ParallelDFS::CheckExitLocal(unsigned long edge) {
  if (local_edges_.size() != 0) {
    auto it = find(local_edges_.begin(), local_edges_.end(), edge);
    if (it != local_edges_.end()) return true;
  }
  return false;
}

void ParallelDFS::ResetLocalVisitedMap() {
  if (local_visited_map_.size() == 0) {
    local_visited_map_.resize(local_edges_.size());
  }
  for (size_t i = 0; i < local_visited_map_.size(); i++) {
    local_visited_map_[i] = false;
  }
}

void ParallelDFS::ResetGlobalVisitedMap() {
  if (global_visited_map_.size() == 0) {
    unsigned long local_size = local_edges_.size();
    unsigned long global_size = 0;
    MPI_Allreduce(&local_size, &global_size, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    global_visited_map_.resize(global_size);
  }

  for (size_t i = 0; i < global_visited_map_.size(); i++) {
    global_visited_map_[i] = false;
  }

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
  unsigned long local_n = local_edges_.size(); 
  for (size_t i = 0; i < local_n; i++) local_paths.push_back(0);
  for (size_t i = 0; i < local_n; i++) {
    ResetLocalVisitedMap();
    local_visited_map_[local_edges_[i] - local_edges_[0]] = true;
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
  delete []num_buf;
}

void ParallelDFS::GlobalDFS(MPI_Win& win, unsigned long cur_end_edge,
                            vector<unsigned long>& single_global_path,
                            vector<vector<unsigned long> >& global_paths) {
  
  // int my_id, process_num;
  // MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  // MPI_Comm_size(MPI_COMM_WORLD, &process_num);

  // if (my_id == 0) {
  //   cout << "#0 cur_end_edge = " << cur_end_edge << endl;
  // }

  if (cur_end_edge == kGlobalEndEdge) {
    global_paths.push_back(single_global_path);
    return;
  }
  if (global_visited_map_[cur_end_edge]) return;
  global_visited_map_[cur_end_edge] = true;

  int target_process = GetEdgeOwner(cur_end_edge);
  vector<unsigned long> remote_path;
  GetPathFromRemote(target_process, cur_end_edge, remote_path, win);

  for (size_t i = 0; i < remote_path.size(); i+=2) {
    if (remote_path[i] != kGlobalEndEdge) {
      // if (global_visited_map_[remote_path[i]]) continue;
      // global_visited_map_[remote_path[i]] = true;
      unsigned long cur_len = single_global_path.back();
      single_global_path.push_back(remote_path[i]);
      single_global_path.push_back(remote_path[i + 1] + cur_len);
      GlobalDFS(win, remote_path[i], single_global_path, global_paths);
      single_global_path.pop_back();
      single_global_path.pop_back();
      // global_visited_map_[remote_path[i]] = false;
    } else {
      unsigned long cur_len = single_global_path.back();
      single_global_path.push_back(remote_path[i]);
      single_global_path.push_back(remote_path[i + 1] + cur_len);
      GlobalDFS(win, remote_path[i], single_global_path, global_paths);
      single_global_path.pop_back();
      single_global_path.pop_back();
    }
  }

  global_visited_map_[cur_end_edge] = false;
}

void ParallelDFS::GlobalPathSearch(MPI_Win& win,
                                   vector<vector<unsigned long> >& global_paths,
                                   vector<unsigned long>& local_paths) {
  for (size_t i = 0; i < local_edges_.size(); i++) {
    unsigned long edge = local_edges_[i];
    vector<unsigned long> single_global_path;
    if (CheckGlobalStartEdge(edge)) {
      ResetGlobalVisitedMap();
      single_global_path.push_back(edge);
      // global_visited_map_[edge] = true;
      vector<unsigned long> local_end_edges;
      vector<unsigned long> paths_len;
      QueryLocalPath(local_paths, edge, local_end_edges, paths_len);

      for (size_t j = 0; j < local_end_edges.size(); j++) {
        // cout << "global path search " << edge << " " << local_end_edges[j] << endl;
        single_global_path.push_back(local_end_edges[j]);
        single_global_path.push_back(paths_len[j]);
        GlobalDFS(win, local_end_edges[j], single_global_path, global_paths);
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

  // master collect path size message
  vector<vector<unsigned long> > paths_size;
  MPI_Request* request = new MPI_Request[num_process];
  MPI_Status* status = new MPI_Status[num_process];

  for (size_t i = 0; i < num_process; i++)
    request[i] = MPI_REQUEST_NULL;

  if (my_id == master_process) {
    for (int send_process = 0; send_process < num_process; send_process++) {
      vector<unsigned long> tmp;
      paths_size.push_back(tmp);
      paths_size[send_process].resize(num_paths[send_process]);
      if (send_process != master_process) {
        MPI_Irecv(paths_size[send_process].data(), num_paths[send_process],
                  MPI_UNSIGNED_LONG, send_process, send_process, MPI_COMM_WORLD,
                  &request[send_process]);
      }
    }
  } else {
      vector<unsigned long> paths_msg_len;
      paths_msg_len.resize(num_local_paths);
      for (size_t i = 0; i < num_local_paths; i++) {
        paths_msg_len[i] = global_paths[i].size();
      }
      MPI_Send(paths_msg_len.data(), num_local_paths,
               MPI_UNSIGNED_LONG, master_process, my_id, MPI_COMM_WORLD);
  }

  MPI_Waitall(num_process, request, status);
  MPI_Barrier(MPI_COMM_WORLD);
  delete []request;
  delete []status;

  // if (my_id == master_process) {
  //   cout << "============ paths size =================" << endl;
  //   for (size_t i = 0; i < paths_size.size(); i++) {
  //     cout << "process # " << i << endl;
  //     for (size_t j = 0; j < paths_size[i].size(); j++) {
  //       cout << paths_size[i][j] << " ";
  //     }
  //     cout << endl;
  //   }
  // }

  unsigned long comm_times = 0; // master communication times
  if (my_id == master_process) {
    for (size_t i = 0; i < num_process; i++) {
      if (i != master_process) comm_times += num_paths[i];
    }
  } 
  MPI_Bcast(&comm_times, 1, MPI_UNSIGNED_LONG, master_process, MPI_COMM_WORLD);
  // cout << "comm_times = " << comm_times << endl;
  request = new MPI_Request[comm_times];
  status = new MPI_Status[comm_times];
  for (size_t i = 0; i < num_process; i++)
    request[i] = MPI_REQUEST_NULL;

  // master collect path message
  if (my_id == master_process) {
    int msg_size;
    unsigned long comm_count = 0;
    for (int send_process = 0; send_process < num_process; send_process++) {
      if (send_process != master_process) {
        for (size_t i = 0; i < num_paths[send_process]; i++) {
          vector<unsigned long> single_path;
          global_paths.push_back(single_path);
          msg_size = paths_size[send_process][i];
          global_paths.back().resize(msg_size);
          MPI_Irecv(global_paths.back().data(), msg_size, MPI_UNSIGNED_LONG,
                   send_process, i, MPI_COMM_WORLD,
                   &request[comm_count]);
          comm_count++;
        }
      }
    }
  } else {
    for (size_t i = 0; i < global_paths.size(); i++) {
      MPI_Send(global_paths[i].data(), global_paths[i].size(),
               MPI_UNSIGNED_LONG, master_process, i, MPI_COMM_WORLD);
    }
  }
  MPI_Waitall(num_process, request, status);
  MPI_Barrier(MPI_COMM_WORLD);
  delete []request;
  delete []status;

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

void ParallelDFS::SinglePathSeqSearch(vector<unsigned long>& global_path) {
  int my_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  unsigned long start_edge = global_path[0];
  unsigned long pre_len = 0;
  unsigned long query_times = 0;
  for (size_t i = 1; i < global_path.size(); i += 2) {
    int remote_process = GetEdgeOwner(start_edge);
    vector<unsigned long> path_seq_segment;
    if (remote_process == my_id) {
      vector<unsigned long> query_msg{start_edge, global_path[i],
                                      global_path[i + 1] - pre_len};
      // cout << "master start local search: " << query_msg[0] << " " <<
      // query_msg[1] << " " << query_msg[2] << endl;
      local_path_seq_complete_ = false;
      path_seq_segment.push_back(query_msg[0]);
      ResetLocalVisitedMap();
      local_visited_map_[query_msg[0] - local_edges_[0]] = true;
      LocalPathSeqSearch(query_msg, query_msg[0], 0, path_seq_segment);

      // cout << "path_seq_segment size = " << path_seq_segment.size() << endl;
      // for (size_t j = 0; j < path_seq_segment.size(); j++) {
      //   cout << "path_segament = " << path_seq_segment[j] << endl;
      // }
      path_output_.ParallelWrite(path_seq_segment, pre_len);
      // cout << "parallelwirte done" << endl;
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
  int num_paths = 0;
  if (my_id == master_process) num_paths = global_paths.size();
  MPI_Bcast(&num_paths, 1, MPI_INT, master_process, MPI_COMM_WORLD);
  
  for (size_t i = 0; i < num_paths; i++) {
    // if (my_id == master_process) {
    //   cout << "============single path search ===============" << endl;
    //   for (size_t j = 0; j < global_paths[i].size(); j++) {
    //     cout << global_paths[i][j] << " ";
    //   }
    //   cout << endl;
    // }
    vector<unsigned long> single_path;
    unsigned long path_len;
    if (my_id == master_process) {
      single_path = global_paths[i];
      path_len = single_path.size();
    }
    MPI_Bcast(&path_len, 1, MPI_UNSIGNED_LONG, master_process, MPI_COMM_WORLD);
    if (my_id != master_process) {
      single_path.resize(path_len);
    }
    MPI_Bcast(single_path.data(), path_len, MPI_UNSIGNED_LONG, master_process, MPI_COMM_WORLD);
    path_output_.CreateFile(i);
    SinglePathSeqSearch(single_path);
    path_output_.CloseFile();
    MPI_Barrier(MPI_COMM_WORLD);
  }
}
void ParallelDFS::run() {
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
  CollectGlobalPath(master_process, global_paths);

  PathSeqSearch(master_process, global_paths);


  unsigned long num_paths;
  if (my_id == master_process) num_paths = global_paths.size();
  MPI_Bcast(&num_paths, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  for (size_t i = 0; i < num_paths; i++) {
    path_output_.OpenFile(i);
    if (my_id == master_process) path_output_.SequenceRead(i);
    path_output_.CloseFile();
  }
  MPI_Win_free(&win);  // free window
}



ParallelStream::ParallelStream(string file_name_prefix) {
  file_name_prefix_ = file_name_prefix;
}
ParallelStream::~ParallelStream() {
}

void ParallelStream::CloseFile() { MPI_File_close(&fh_); }

void ParallelStream::SetPrefix(string prefix) {
  file_name_prefix_ = prefix + "path_";
}
void ParallelStream::ParallelWrite(vector<unsigned long>& path_seq,
                                   unsigned long offset) {

  MPI_Status status;
  MPI_File_write_at(fh_, offset * sizeof(unsigned long), path_seq.data(), path_seq.size(), MPI_UNSIGNED_LONG,
                    &status);
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

void ParallelStream::CreateFile(unsigned long path_id) {
  string file_name;
  file_name = file_name_prefix_ + to_string(path_id);
  MPI_File_open(MPI_COMM_WORLD, file_name.c_str(),
                MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh_);
}
void ParallelStream::OpenFile(unsigned long path_id) {
  string file_name;
  file_name = file_name_prefix_ + to_string(path_id);
  MPI_File_open(MPI_COMM_WORLD, file_name.c_str(), MPI_MODE_RDWR, MPI_INFO_NULL,
                &fh_);
}