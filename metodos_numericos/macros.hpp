#include <iostream>
#include <iomanip> 
#include <chrono>

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
