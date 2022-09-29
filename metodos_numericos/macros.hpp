#ifndef MACROS_HPP
#define MACROS_HPP
#include <iostream>
#include <iomanip> 
#include <chrono>
#include <fstream>
#include "matrix_like/real_matrix.hpp"

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
#define ANNOUNCE_TEST(msg) \
(*(macros::out)) << "#*********************************#" << std::endl;\
(*(macros::out)) << "# * TEST:" << msg << std::endl; \
(*(macros::out)) << "#*********************************#" << std::endl;
#define start_timer(start) start = high_resolution_clock::now()
#define end_timer(end) end = high_resolution_clock::now()
#define print_timer(start,end) (*(macros::out)) << "time: " << duration_cast<microseconds>(end - start).count() << "micro s" << std::endl
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
template <typename T>
void print_array(T &array, size_t size)
{
	(*(macros::out)) << std::fixed;
	for (size_t i = 0; i < size; ++i)
	{
		(*(macros::out)) << std::setfill(' ') << std::setw(7) << std::setprecision(2) << array[i] << " ";
	}
	(*(macros::out)) << std::endl;
}
template <typename T>
void print_matrix(T &matrix, size_t size)
{
	for (size_t i = 0; i < size; ++i)
	{
		auto a = matrix[i];
		print_array(a, size);
	}
}

#endif /* MACROS_HPP */
