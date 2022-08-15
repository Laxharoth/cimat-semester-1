#ifndef CANVAS_FACADE_CPP
#define CANVAS_FACADE_CPP
#include "canvas_facade.hpp"

namespace canvas_facade{
point choose_origin_helper( const double &min_x, const double &min_y,
                                const double &scale_x, const double &scale_y,
                                const double &canvas_heigh){
        return point{ min_x * scale_x, min_y * scale_y + canvas_heigh , false };
}
double choose_scale_helper(const double &min, const double &max, const double &space){
    return space / (max - min);
}
}

#endif /* CANVAS_FACADE_CPP */
