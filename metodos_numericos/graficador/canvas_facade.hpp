#ifndef CANVAS_FACADE_HPP
#define CANVAS_FACADE_HPP

#include "point.hpp"

#include <gtk/gtk.h>

namespace canvas_facade{

class CanvasFacade{
    public:
    point origin;
    point min_value;
    point max_value;
    double scale_x;
    double scale_y;
    void set_origin(const point &new_origin);
    void draw_line(void *drawing_context, const point &point_1, const point &point_2) const;
    void draw_multiline(void *drawing_context, point *points, size_t lenght) const;
};
point choose_origin_helper( const double &min_x, const double &min_y,
                            const double &scale_x, const double &scale_y,
                            const double &canvas_heigh);
double choose_scale_helper(const double &min, const double &max, const double &space);
point choose_origin_helper( const double &min_x, const double &min_y,
                                const double &scale_x, const double &scale_y,
                                const double &canvas_heigh){
        return point{ min_x * scale_x, min_y * scale_y + canvas_heigh , false };
    }
double choose_scale_helper(const double &min, const double &max, const double &space){
    return space / (max - min);
}
}

#endif /* CANVAS_FACADE_HPP */
