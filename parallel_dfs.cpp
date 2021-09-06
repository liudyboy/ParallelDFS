/*
 * @Brief: file description
 * @Author: liudy
 * @Email: deyin.liu@nscc-gz.cn
 * @Date: 2021-07-25 16:39:41
 * @LastEditors: liudy
 * @LastEditTime: 2021-08-05 10:27:03
 */
#include "parallel_dfs.h"

#include <pthread.h>
#include <sys/stat.h>

#include <algorithm>
#include <cstdio>
#include <map>

using std::cout;
using std::endl;
using std::ios;
using std::map;
using std::to_string;

// caution: actually the values is 18446744073709551615 (Maximum of unsigned
// long)
const unsigned long kGlobalEndEdge = -1;  // 18446744073709551615
const unsigned long kStopSignal = -2;     // 18446744073709551614
const unsigned long kInvalidEdge = -3;    // 18446744073709551613
const unsigned long kSearchEnd = -4;      // 18446744073709551612
const unsigned long Seed = 2021;

// Test example 1

// void ParallelDFS::GetNextEdge(unsigned long edge, vector<unsigned long>&
// next_edges) {
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
//   if (edge == 0 || edge == 3 || edge == 4 || edge == 13 || edge == 12) return
//   true; return false;
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

// void ParallelDFS::GetNextEdge(unsigned long edge, vector<unsigned long>&
// next_edges) {
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
// void ParallelDFS::GetNextEdges(unsigned long edge, vector<unsigned long>&
// next_edges) {
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
void ParallelDFS::GetNextEdges(unsigned long edge,
                               vector<unsigned long>& next_edges) {
  map<unsigned long, vector<unsigned long> > graph;
  vector<unsigned long> children0{1};
  vector<unsigned long> children1{2, 4};
  vector<unsigned long> children2{3, 9};
  vector<unsigned long> children3{1};
  vector<unsigned long> children4;
  vector<unsigned long> children5{1};
  vector<unsigned long> children6{5, 7};
  vector<unsigned long> children7{3, 9};
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
  if (edge >= 0 && edge < 5) return 0;
  return 1;
}

ParallelDFS::ParallelDFS() {
  output_dir_ = "paths/";
}


void ParallelDFS::SetOutDir(string address) {
  output_dir_ = address;
}

unsigned long ParallelDFS::SearchNextEdge(const unsigned long edge,
                                          const unsigned long path_id,
                                          const unsigned start_edge,
                                          const unsigned long cur_fork,
                                          unsigned long& next_fork,
                                          unsigned long& hash_value) {
  unsigned long local_id = edge - local_edges_[0];
  if (next_edges_[local_id].size() == 0) return kGlobalEndEdge;
  bool activate_fork = cur_fork == edge ? true : false;
  unsigned long next_edge = local_msg_[local_id].Query(
      path_id, start_edge, activate_fork, next_edges_[local_id], hash_value);
  if (next_edge != next_edges_[local_id].back() && next_edge != kInvalidEdge)
    next_fork = edge;
  if (next_edge == kInvalidEdge && edge == start_edge) next_edge = kSearchEnd;

  unsigned long array[2]{hash_value, next_edge};
  hash_value = fasthash(array, 2 * sizeof(unsigned long), Seed);
  return next_edge;
}

unsigned long ParallelDFS::QueryNextEdge(const unsigned long edge,
                                         PathMsg& path_msg,
                                         unsigned long* fork_edge,
                                         unsigned long& hash_value) {
  unsigned long query_msg[6]{edge,         path_msg.Data(0), path_msg.Data(1),
                             fork_edge[0], fork_edge[1],     hash_value};
  int owner = GetEdgeOwner(edge);
  MPI_Send(query_msg, 6, MPI_UNSIGNED_LONG, owner, 0, MPI_COMM_WORLD);

  unsigned long answer[3];  // {next_edge, next_fork, hash_value}
  MPI_Recv(answer, 3, MPI_UNSIGNED_LONG, owner, 1, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  fork_edge[1] = answer[1];
  hash_value = answer[2];
  return answer[0];
}

void ParallelDFS::Server() {
  unsigned long msg[6];
  unsigned long answer[3];
  MPI_Status status;
  while (true) {
    MPI_Recv(msg, 6, MPI_UNSIGNED_LONG, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,
             &status);

    if (msg[0] == kStopSignal) return;
    unsigned long local_id = msg[0] - local_edges_[0];
    if (next_edges_[local_id].size() == 0) {
      GetNextEdges(msg[0], next_edges_[local_id]);
    }
    unsigned long next_edge =
        SearchNextEdge(msg[0], msg[1], msg[2], msg[3], msg[4], msg[5]);
    answer[0] = next_edge;
    answer[1] = msg[4];  // next_fork
    answer[2] = msg[5];  // hash_value
    MPI_Send(answer, 3, MPI_UNSIGNED_LONG, status.MPI_SOURCE, 1,
             MPI_COMM_WORLD);
  }
}

void* ParallelDFS::ServerWrapper(void* object) {
  reinterpret_cast<ParallelDFS*>(object)->Server();
  return nullptr;
}

void ParallelDFS::Client() {
  local_path_count_ = 0;
  int my_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

  unsigned long fork_edges[2]{kInvalidEdge,
                              kInvalidEdge};  // <cur_fork, next_fork>
  PathFile file_output;
  file_output.SetOutDir(output_dir_);
  file_output.CreateDir();
  PathMsg path_msg;
  unsigned long cur_edge;
  unsigned long hash_value;
  for (size_t i = 0; i < local_edges_size_; i++) {
    unsigned long edge = i + local_edges_[0];
    if (CheckGlobalStartEdge(edge)) {
      cur_edge = edge;
      fork_edges[1] = edge;
      local_path_count_++;
      path_msg.Set(local_path_count_, edge);
      hash_value = fasthash(&edge, sizeof(unsigned long), Seed);

      while (true) {
        file_output.Output(cur_edge);
        cur_edge = QueryNextEdge(cur_edge, path_msg, fork_edges, hash_value);
        if (cur_edge == kGlobalEndEdge) {
          file_output.SavePath();
          local_path_count_++;
          path_msg.Set(local_path_count_, edge);
          cur_edge = edge;
          fork_edges[0] = fork_edges[1];
          fork_edges[1] = edge;
          hash_value = fasthash(&edge, sizeof(unsigned long), Seed);
        } else if (cur_edge == kInvalidEdge) {
          file_output.RemovePath();
          local_path_count_++;
          path_msg.Set(local_path_count_, edge);
          cur_edge = edge;
          fork_edges[0] = fork_edges[1];
          fork_edges[1] = edge;
          hash_value = fasthash(&edge, sizeof(unsigned long), Seed);
        } else if (cur_edge == kSearchEnd) {
          file_output.Clear();
          break;
        }
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  StopServer();
}

void ParallelDFS::StopServer() {
  int my_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  unsigned long msg[6]{kStopSignal, 0, 0, 0, 0, 0};
  MPI_Send(msg, 6, MPI_UNSIGNED_LONG, my_id, 0, MPI_COMM_WORLD);
}

void ParallelDFS::run() {
  int my_id, process_num;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Comm_size(MPI_COMM_WORLD, &process_num);

  GetLocalEdges(local_edges_);
  local_edges_size_ = local_edges_[1] - local_edges_[0] + 1;
  local_msg_ = new LocalMsg[local_edges_size_];
  next_edges_ = new vector<unsigned long>[local_edges_size_];

  pthread_t server_thread;
  pthread_create(&server_thread, NULL, &ParallelDFS::ServerWrapper, this);
  Client();
  pthread_join(server_thread, NULL);

  delete[] local_msg_;
  delete[] next_edges_;
  // PathFile file_output;
  // for (unsigned long i = 1; i < 50; i++) {
  //   if (i % 10 == 0)
  //     file_output.NewPath();
  //   if (i == 20)
  //     file_output.RemovePath();
  //   file_output.Output(i);
  // }
}

/****************** PathMsg class *******************************/

PathMsg::PathMsg() {}
PathMsg::PathMsg(unsigned long* data) {
  data_[0] = data[0];
  data_[1] = data[1];
}

PathMsg::PathMsg(unsigned long path_id, unsigned long start_edge) {
  int process_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
  unsigned long tmp = process_id;
  data_[0] = tmp << 32;
  data_[0] += path_id;
  data_[1] = start_edge;
}

void PathMsg::Set(unsigned long path_id, unsigned long start_edge) {
  int process_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
  unsigned long tmp = process_id;
  data_[0] = tmp << 32;
  data_[0] += path_id;
  data_[1] = start_edge;
}

void PathMsg::Set(unsigned long path_id) {
  int process_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
  unsigned long tmp = process_id;
  data_[0] = tmp << 32;
  data_[0] += path_id;
}

unsigned long* PathMsg::Data() { return data_; }

unsigned long PathMsg::Data(unsigned long index) { return data_[index]; }

/****************** PathFile class *******************************/
PathFile::PathFile() {
  data_ = new unsigned long[FILE_CACHE_LIMIT];
  count_ = 0;
  path_count_ = 0;
  file_name_ = "paths/";
}

PathFile::~PathFile() {
  if (fout.is_open()) SavePath();
  delete[] data_;
}

void PathFile::CreateDir() {
  cout << file_name_ << endl;
  int my_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  string dir = file_name_ + to_string(my_id);
  mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
}

void PathFile::CreateFile() {
  int my_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  path_count_++;
  string file_name = file_name_ + to_string(my_id) + "/path_" +
                     to_string(path_count_) + ".csv";
  fout.open(file_name, ios::out | ios::app);
}

void PathFile::Output(unsigned long data) {
  data_[count_] = data;
  count_++;
  if (count_ >= FILE_CACHE_LIMIT) {
    Flush();
  }
}

void PathFile::Flush() {
  if (count_ == 0) return;
  if (!fout.is_open()) {
    CreateFile();
  }
  for (int i = 0; i < count_; i++) {
    fout << data_[i] << ",";
  }
  count_ = 0;
  fout.flush();
}

void PathFile::SavePath() {
  Flush();
  fout << endl;
  fout.close();
  count_ = 0;
}

void PathFile::RemovePath() {
  int my_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  count_ = 0;
  fout.close();
  string file_name = file_name_ + to_string(my_id) + "/path_" +
                     to_string(path_count_) + ".csv";
  remove(file_name.c_str());
  path_count_--;
}

void PathFile::Clear() { count_ = 0; }

void PathFile::SetOutDir(string address) {
  file_name_ = address;
}

/****************** LocalMsg class *******************************/
LocalMsg::LocalMsg() : count_(0), next_(nullptr) {}

LocalMsg::LocalMsg(const unsigned long path_id, const unsigned long start_edge,
                   const unsigned long hash_value)
    : start_edge_(start_edge),
      path_id_(path_id),
      hash_value_(hash_value),
      count_(0),
      next_(nullptr) {}

void LocalMsg::Link(LocalMsg* msg) {
  msg->next_ = next_;
  next_ = msg;
}

unsigned long LocalMsg::Query(const unsigned long path_id,
                              const unsigned start_edge, bool activate_fork,
                              vector<unsigned long>& next_edges,
                              const unsigned long hash_value) {
  unsigned long next_edge = kInvalidEdge;
  bool first_query = true, circle_exit = false;
  unsigned short count = 0;
  LocalMsg* cur_msg = next_;
  while (cur_msg) {
    if (cur_msg->path_id_ == path_id) {
      circle_exit = true;
      break;
    }

    if (cur_msg->start_edge_ == start_edge) {
      first_query = false;
      if (activate_fork && cur_msg->hash_value_ == hash_value) {
        cur_msg->count_++;
      }
      if (hash_value != cur_msg->hash_value_) cur_msg->count_ = 0;

      count = cur_msg->count_;
      cur_msg->path_id_ = path_id;
      cur_msg->hash_value_ = hash_value;
      break;
    }
    cur_msg = cur_msg->next_;
  }

  if (circle_exit || count >= next_edges.size()) {
    return kInvalidEdge;
  }

  if (first_query) {
    LocalMsg* msg_block = new LocalMsg(path_id, start_edge, hash_value);
    Link(msg_block);
    next_edge = next_edges[0];
    return next_edge;
  }

  return next_edges[count];
}

uint64_t fasthash(const void* buf, size_t len, uint64_t seed) {
  const uint64_t m = 0x880355f21e6d1965ULL;
  const uint64_t* pos = (const uint64_t*)buf;
  const uint64_t* end = pos + (len / 8);
  const unsigned char* pos2;
  uint64_t h = seed ^ (len * m);
  uint64_t v;

  while (pos != end) {
    v = *pos++;
    h ^= mix(v);
    h *= m;
  }

  pos2 = (const unsigned char*)pos;
  v = 0;

  switch (len & 7) {
    case 7:
      v ^= (uint64_t)pos2[6] << 48;
    case 6:
      v ^= (uint64_t)pos2[5] << 40;
    case 5:
      v ^= (uint64_t)pos2[4] << 32;
    case 4:
      v ^= (uint64_t)pos2[3] << 24;
    case 3:
      v ^= (uint64_t)pos2[2] << 16;
    case 2:
      v ^= (uint64_t)pos2[1] << 8;
    case 1:
      v ^= (uint64_t)pos2[0];
      h ^= mix(v);
      h *= m;
  }
  return mix(h);
}
