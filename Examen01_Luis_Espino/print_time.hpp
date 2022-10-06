#include <chrono>
#include <iostream>

using namespace std::chrono;

namespace macros{
  std::ostream *out = &(std::cout);
}

#define measure_time(what) \
{\
  auto start = high_resolution_clock::now();\
  what;\
  auto end = high_resolution_clock::now();\
  auto duration = duration_cast<microseconds>(end - start);\
  (*(macros::out)) << "time: " << duration.count() << "micro s" << std::endl;\
}