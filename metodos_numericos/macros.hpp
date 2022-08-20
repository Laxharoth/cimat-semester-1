#include <iostream>
#include <iomanip> 

#define matrix_to_ptr(matrix,ptr,size) for(int i = 0; i < size; ++i) ptr[i] = matrix[i];
#define print_test( actual, expected ) \
    std::cout << "actual: " <<  actual << " :: "; \
    std::cout << "expected: " <<  expected << std::endl;
#define ANNOUNCE_TEST(msg) std::cout << "TEST:" << msg << std::endl;
#define measure_time(what) \
{\
  auto start = high_resolution_clock::now();\
  what;\
  auto end = high_resolution_clock::now();\
  auto duration = duration_cast<nanoseconds>(end - start);\
  std::cout << "time: " << duration.count() << "nanos" << std::endl;\
}
