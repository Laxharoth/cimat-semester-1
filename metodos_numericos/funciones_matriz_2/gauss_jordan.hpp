#ifndef GAUSS_JORDAN_HPP
#define GAUSS_JORDAN_HPP

#include "funcion_matriz.hpp"
#include <string.h>

void gauss_jordan( matrix_like<double> &matriz, array_like<double> &variables, array_like<double> &resultados, const int &size );

#endif /* GAUSS_JORDAN_HPP */
