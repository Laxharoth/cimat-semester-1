#ifndef FUNCION_MATRIZ_HPP
#define FUNCION_MATRIZ_HPP

#include "matrix_like/matrix_like.tcc"
#include "matrix_like/real_matrix.hpp"
#include "factorizacion.hpp"

#include <cstdlib>

using mymtx::RealMatrix;
using mymtx::RealVector;

void solucion_diagonal(RealMatrix &matriz, RealVector &incognitas, RealVector &result);
double determinante_diagonal(RealMatrix &matriz_diagonal);
void inversa_diagonal(RealMatrix &matriz_diagonal,RealVector &inversa);
double determinante_triangular(RealMatrix &matriz_triangular);
void solucion_triangular_inf( RealMatrix &matriz, RealVector &incognitas, RealVector &result,bool);
void solucion_triangular_inf( RealMatrix &matriz, RealVector &incognitas, RealVector &result);
void solucion_triangular_sup( RealMatrix &matriz, RealVector &incognitas, RealVector &result,bool);
void solucion_triangular_sup( RealMatrix &matriz, RealVector &incognitas, RealVector &result);
void gauss( RealMatrix &matriz, RealVector &variables, RealVector &resultados);
void solucion_LDU( RealMatrix &matriz, RealVector &incognitas, RealVector &result);
void solucion_crout( RealMatrix &matriz, RealVector &incognitas, RealVector &result);
void solucion_doolittle( RealMatrix &matriz, RealVector &incognitas, RealVector &result);

#endif /* FUNCION_MATRIZ_HPP */
