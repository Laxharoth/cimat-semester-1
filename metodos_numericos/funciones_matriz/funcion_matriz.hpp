#ifndef FUNCION_MATRIZ_HPP
#define FUNCION_MATRIZ_HPP

#include "matrix_like/matrix_like.tcc"

#include <cstdlib>

void solucion_diagonal(matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result, size_t size);
double determinante_diagonal(matrix_like<double> &matriz_diagonal, const size_t &size);
void inversa_diagonal(matrix_like<double> &matriz_diagonal,array_like<double> &inversa, const size_t &size);
double determinante_triangular(matrix_like<double> &matriz_triangular, const size_t &size);
void solucion_triangular_inf( matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result, const size_t &size);
void solucion_triangular_sup( matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result, const size_t &size);

#endif /* FUNCION_MATRIZ_HPP */
