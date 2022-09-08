#ifndef FUNCION_MATRIZ_HPP
#define FUNCION_MATRIZ_HPP

#include "matrix_like/real_matrix.hpp"
#include "factorizacion.hpp"

#include <cstdlib>
#include <cmath>

using mymtx::RealMatrix;
using mymtx::RealVector;

struct eigen{
    RealVector vect;
    double val;
};

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
void solve_gauss_seidel(mymtx::RealMatrix &matix, mymtx::RealVector &variables, mymtx::RealVector &solutions, double *error);
void solve_cholesky(mymtx::RealMatrix &cholesky_factored,mymtx::RealVector &variables, mymtx::RealVector &solutions);
double normalize(RealVector &vec);
void power_iteration(const RealMatrix &A, RealVector &V0, RealVector &V1, const double tolerance, double &value, size_t n_values, RealMatrix *_vec_holder, RealVector *_val_holder);
void inverse_power_iteration(const RealMatrix &A, RealVector &V0, RealVector &V1, const double tolerance, double &value, size_t n_values, RealMatrix *_vec_holder, RealVector *_val_holder);
#endif /* FUNCION_MATRIZ_HPP */
