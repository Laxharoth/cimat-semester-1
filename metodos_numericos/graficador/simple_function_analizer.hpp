#ifndef SIMPLE_FUNCTION_ANALIZER_HPP
#define SIMPLE_FUNCTION_ANALIZER_HPP
#ifndef PARTITIONS_NUM
#define PARTITIONS_NUM 100
#endif

#include "point.hpp"

#include <array>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <float.h>

/**
 * @brief Obtiene un conjunto de puntos de una funcion.
 * El numero de puntos esta definido por la macro PARTITIONS_NUM
 * Tambien obtiene el minimo y maximo valor de la funcion en el intervalo
 * proporcionado
 *
 * @param func la funcion a evaluar
 * @param left el valor minimo de "x"
 * @param right el valor maximo de "x"
 * @param min_y Devuelve el valor minimo de "y"
 * @param max_y Devuelve el valor maximo de "y"
 * @param points Devuelve los puntos evaluados.
 */
void analize_single_var_function(double (*func)(double x), const double left,
                                 const double right, double &min_y,
                                 double &max_y,
                                 std::array<point, PARTITIONS_NUM + 1> &points);

#endif
