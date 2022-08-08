#ifndef SIMPLE_FUNCTION_ANALIZER_HPP
#define SIMPLE_FUNCTION_ANALIZER_HPP
#ifndef PARTITIONS_NUM
#define PARTITIONS_NUM 1000
#endif

#include "point.hpp"

#include <cstdlib>
#include <array>
#include <float.h>
#include <exception>

void analize_single_var_function( double (*func)(double x), const double left, const double right ,double &min_y, double &max_y, std::array<point, PARTITIONS_NUM> &points );


#endif 
