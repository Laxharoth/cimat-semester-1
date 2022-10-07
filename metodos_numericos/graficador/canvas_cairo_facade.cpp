#include "canvas_facade.cpp"
#include "canvas_facade.hpp"

namespace canvas_facade {
void CanvasFacade::set_origin(const point &new_origin) {
  this->origin.x = new_origin.x;
  this->origin.y = new_origin.y;
}
void CanvasFacade::draw_line(void *drawing_context, const point &point_1,
                             const point &point_2) const {
  cairo_t *cr = (cairo_t *)drawing_context;
  if (point_1.is_nan || point_2.is_nan)
    return;
  cairo_set_line_width(cr, 1.0);
  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
  cairo_move_to(cr, point_1.x * scale_x - this->origin.x,
                this->origin.y - point_1.y * scale_y);
  cairo_line_to(cr, point_2.x * scale_x - this->origin.x,
                this->origin.y - point_2.y * scale_y);
  cairo_stroke(cr);
}
void CanvasFacade::draw_multiline(void *drawing_context, point *points,
                                  size_t lenght) const {
  cairo_t *cr = (cairo_t *)drawing_context;
  cairo_set_line_width(cr, 1.0);
  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
  size_t start = 0;
  while (start + 1 < lenght && points[start].is_nan) {
    ++start;
  }
  cairo_move_to(cr, points[start].x * scale_x - this->origin.x,
                this->origin.y - points[start].y * scale_y);
  bool skiped = false;
  for (size_t i = start + 1; i < lenght; i++) {
    point &current = points[i];
    if (current.is_nan) {
      skiped = true;
      continue;
    }
    if (skiped) {
      skiped = false;
      cairo_move_to(cr, current.x * scale_x - this->origin.x,
                    this->origin.y - current.y * scale_y);
      continue;
    }
    cairo_line_to(cr, current.x * scale_x - this->origin.x,
                  this->origin.y - current.y * scale_y);
  }
  cairo_stroke(cr);
}
} // namespace canvas_facade