#ifndef CANVAS_FACADE_HPP
#define CANVAS_FACADE_HPP

#include "point.hpp"

#include <gtk/gtk.h>

namespace canvas_facade{

class CanvasFacade{
    public:
    point origin;
    double min_x;
    double min_y;
    void set_origin(const point &new_origin);
    void draw_line(cairo_t *cr, const point &point_1, const point &point_2);
};
point choose_origin_helper( const double &min_x, const double &min_y,
                            const double &canvas_heigh);
}


#endif /* CANVAS_FACADE_HPP */
