#include <iostream>
#include <iomanip> 
#include <chrono>
#include <fstream>

#ifndef COMMA
#define COMMA ,
#endif

namespace macros{
  std::ostream *out = &(std::cout);
}
using namespace std::chrono;

#define strm_out(what) (*(macros::out)) << what << std::endl;
#define matrix_to_ptr(matrix,ptr,size) for(int i = 0; i < size; ++i) ptr[i] = matrix[i];
#define print_test( actual, expected ) \
    (*(macros::out)) << "actual: " <<  actual << " :: "; \
    (*(macros::out)) << "expected: " <<  expected << std::endl;
#define ANNOUNCE_TEST(msg) (*(macros::out)) << "TEST:" << msg << std::endl;
#define measure_time(what) \
{\
  auto start = high_resolution_clock::now();\
  what;\
  auto end = high_resolution_clock::now();\
  auto duration = duration_cast<microseconds>(end - start);\
  (*(macros::out)) << "time: " << duration.count() << "micro s" << std::endl;\
}
#define out_vector(vector,filename) \
{\
  auto file = std::ofstream(filename); \
  file << vector.size << " " << 1 << std::endl;\
  for(auto i = vector.begin(); i != vector.end(); ++i)\
    file << *i << std::endl;\
  file.close();\
}

#define out_matrix(matrix,filename) \
{\
  auto file = std::ofstream(filename); \
  file << matrix.shape_y << " " << matrix.shape_x << std::endl;\
  for(auto j = 0; j < matrix.shape_y; ++j){\
    for(auto i = matrix[j].begin(); i != matrix[j].end(); ++i)\
      file << *i << " ";\
    file << std::endl;\
  }\
  file.close();\
}
