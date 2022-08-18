#include "simple_function_analizer.hpp"

void analize_single_var_function( double (*func)(double), const double left, const double right ,double &min_y, double &max_y, std::array<point, PARTITIONS_NUM+1> &points ){
    min_y = DBL_MAX;
    max_y = -min_y;
    const double increment = (right - left) / PARTITIONS_NUM;
    for(unsigned int i = 0; i <= PARTITIONS_NUM; ++i){
        point &current_point = points[i];
        current_point.is_nan = false;
        const double current_x = left + increment * i;
        current_point.x = current_x;
        try{
            const double current_y = func(current_x);
            current_point.y = current_y;
            if(std::isinf(current_y) || std::isnan(current_y)){
                current_point.is_nan = true;
                continue;
            }
        }catch(...){
            current_point.is_nan = true;
            continue;
        }
        if(min_y > current_point.y)min_y = current_point.y;
        if(max_y < current_point.y)max_y = current_point.y;
    }
}