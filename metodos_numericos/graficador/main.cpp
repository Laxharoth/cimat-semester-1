#include "point.hpp"
#include "simple_function_analizer.hpp"

#include <array>
#include <iostream>

int main(){
    double min_y = 0.0;
    double max_y = 0.0;
    std::array<point, PARTITIONS_NUM> points;

    analize_single_var_function([](double x){ return x*x; }, -10, 10, min_y, max_y, points);

    std::cout << "min_y:"<<min_y << std::endl;
    std::cout << "max_y:"<<max_y << std::endl;
    std::cout << "Points:" << std::endl;
    for(size_t i=0; i<points.size(); ++i){
        std::cout << "x["<<i<<"] = "<< points[i].x << std::endl;
        std::cout << "y["<<i<<"] = "<< points[i].y << std::endl;
        std::cout << "-------------" << std::endl;
    }

    return 0;
}