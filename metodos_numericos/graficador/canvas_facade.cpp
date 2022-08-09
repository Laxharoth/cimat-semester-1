#include "canvas_facade.hpp"
namespace canvas_facade{
    void CanvasFacade::set_origin(const point &new_origin){
        this->origin.x = new_origin.x;
        this->origin.y = new_origin.y;
    }
    void CanvasFacade::draw_line(cairo_t *cr, const point &point_1, const point &point_2){
    cairo_set_line_width(cr, 1.0);
        cairo_set_source_rgb(cr,0.0,0.0,0.0);
        cairo_move_to(cr, point_1.x - this->origin.x, this->origin.y - point_1.y);
        cairo_line_to(cr, point_2.x - this->origin.x, this->origin.y - point_2.y);
        //g_printerr("(%f,%f)->(%f,%f)\n",this->origin.x + point_1.x * x_scale, this->origin.y - point_1.y * y_scale, this->origin.x + point_2.x * x_scale, this->origin.y - point_2.y * y_scale);
        cairo_stroke(cr);
    }
    point choose_origin_helper( const double &min_x, const double &min_y,
                                const double &canvas_heigh){
        return point{ min_x, min_y + canvas_heigh , false };
    }
    double choose_scale_helper(const double &min, const double &max, const double &space){
        return space / (max - min);
    }
}