#ifndef FUNCION_MATRIZ_HPP
#define FUNCION_MATRIZ_HPP

#include "matrix_like/matrix_like.tcc"
#include "factorizacion.hpp"

#include <cstdlib>

void solucion_diagonal(matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result);
double determinante_diagonal(matrix_like<double> &matriz_diagonal);
void inversa_diagonal(matrix_like<double> &matriz_diagonal,array_like<double> &inversa);
double determinante_triangular(matrix_like<double> &matriz_triangular);
void solucion_triangular_inf( matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result);
void solucion_triangular_sup( matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result);
void gauss( matrix_like<double> &matriz, array_like<double> &variables, array_like<double> &resultados);
void solucion_LDU( matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result);
void solucion_crout( matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result);
void solucion_doolittle( matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result);

#endif /* FUNCION_MATRIZ_HPP */
