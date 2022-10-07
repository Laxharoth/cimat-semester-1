#ifndef CANVAS_FACADE_HPP
#define CANVAS_FACADE_HPP

#include "point.hpp"

#include <gtk/gtk.h>

namespace canvas_facade {

class CanvasFacade {
  // El valor en pixeles donde se encuentra el origen despues de aplicar las
  // escalas
  point origin;

public:
  // Los valores minimos de x y f(x)
  point min_value;
  // Los valores maximos de x y f(x)
  point max_value;
  // La escala aplicada a los valores de x
  double scale_x;
  // La escala aplicada a los valores de f(x)
  double scale_y;
  /**
   * @brief Asigna el valor de origin
   * @param new_origin el nuevo valor para origin
   */
  void set_origin(const point &new_origin);
  /**
   * @brief Dibuja una linea entre el point_1 y point_2 solo si point_1 y
   * point_2 no son nan
   *
   * @param drawing_context el contexto de dibujo
   * @param point_1 el punto de origen de la linea
   * @param point_2 el punto de destino de la linea
   */
  void draw_line(void *drawing_context, const point &point_1,
                 const point &point_2) const;
  /**
   * @brief Dibuja lineas entre cada par de puntos consecutivos.
   * Solo dibuja la linea si ninguno de los dos puntos es nan.
   *
   * @param drawing_context el contexto de dibujo
   * @param points la lista de puntos a dibujar
   * @param lenght la cantidad de puntos a dibujar
   */
  void draw_multiline(void *drawing_context, point *points,
                      size_t lenght) const;
};
/**
 * @brief Obtiene el punto "origin" para un objeto CanvasFacade
 *
 * @param min_x El valor minimo de x
 * @param min_y El valor minimo de y
 * @param scale_x La escala en x utilizada por el objeto CanvasFacade
 * @param scale_y  La escala en x utilizada por el objeto CanvasFacade
 * @param canvas_heigh La altura en pixeles del contexto de dibujo
 * @return point El punto "origin" para el objeto CanvasFacade
 */
point choose_origin_helper(const double &min_x, const double &min_y,
                           const double &scale_x, const double &scale_y,
                           const double &canvas_heigh);
/**
 * @brief Obtiene la escala para que una linea que va desde un valor a otro
 * quepa en un espacio especifico
 *
 * @param min El valor minimo
 * @param max El valor maximo
 * @param space El espacio en pixeles
 * @return double
 */
double choose_scale_helper(const double &min, const double &max,
                           const double &space);
} // namespace canvas_facade

#endif /* CANVAS_FACADE_HPP */
