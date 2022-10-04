#ifndef FUNCION_MATRIZ_HPP
#define FUNCION_MATRIZ_HPP

#include "matrix_like/real_matrix.hpp"

#include <cstdlib>
#include <cmath>
#include <random>
#include <vector>
#include <future>

using mymtx::RealMatrix;
using mymtx::RealVector;

struct eigen{
    RealVector vect;
    double val;
};

void solucion_diagonal(const RealMatrix &matriz, RealVector &incognitas, RealVector &result);
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
void randomize(RealVector& vec);
void power_iteration(const RealMatrix &A, RealVector &V1, const double tolerance, double &value, size_t n_values, RealMatrix *_vec_holder, RealVector *_val_holder,const size_t max_iter);
void inverse_power_iteration(const RealMatrix &A, RealVector &V1, const double tolerance, double &value, size_t n_values, RealMatrix *_vec_holder, RealVector *_val_holder, const size_t);
void jacobi_eigen(mymtx::RealMatrix &A, mymtx::RealVector &e, mymtx::RealMatrix  &U, const unsigned max_iter);
void subspace_pow(const mymtx::RealMatrix &A,mymtx::RealMatrix &I,mymtx::RealVector &eig);
void subspace_ipow(const mymtx::RealMatrix &A,mymtx::RealMatrix &I_t,mymtx::RealVector &eig);
void rayleigh_method(const mymtx::RealMatrix &A,mymtx::RealVector &V1, double &val);
void solve_qr(const mymtx::RealMatrix &Q, mymtx::RealMatrix &R, mymtx::RealVector &var, const mymtx::RealVector &res);
void conjugate_gradient(const mymtx::RealMatrix &A, mymtx::RealVector &var, mymtx::RealVector &result);
void conjugate_gradient_jacobi(const mymtx::RealMatrix &A, mymtx::RealVector &var, mymtx::RealVector &result);

using mymtx::RealMatrix;

void crout(RealMatrix &matriz, RealMatrix &matriz_inferior, RealMatrix &matriz_superior);
void crout_as_band(RealMatrix &A_mtx, RealMatrix &L_mtx, RealMatrix &U_mtx,int height, int width);
void crout_tridiagonal(RealMatrix &matriz, RealMatrix &matriz_inferior, RealMatrix &matriz_superior);
void doolittle(RealMatrix &matriz, RealMatrix &matriz_inferior, RealMatrix &matriz_superior);
void doolittle_as_band(RealMatrix &A_mtx, RealMatrix &L_mtx, RealMatrix &U_mtx,int height, int width);
void doolittle_tridiagonal(RealMatrix &matriz, RealMatrix &matriz_inferior, RealMatrix &matriz_superior);
void LDU_factor(RealMatrix &matriz, RealMatrix &matriz_inferior, RealMatrix &matriz_diagonal, RealMatrix &matriz_superior);
void LDU_factor_tridiagonal(RealMatrix &matriz, RealMatrix &matriz_inferior, RealMatrix &matriz_diagonal, RealMatrix &matriz_superior);
void factor_cholesky(const mymtx::RealMatrix &matix, mymtx::RealMatrix &triangular);
void factor_cholesky_as_band(const mymtx::RealMatrix &matix, mymtx::RealMatrix &triangular, int heigh);
void factor_cholesky_tridiag(mymtx::RealMatrix &matix, mymtx::RealMatrix &triangular);
void qr_decomposition(const mymtx::RealMatrix& A, mymtx::RealMatrix&Q, mymtx::RealMatrix&R);

class cant_factor_exception : public std::exception{
    using std::exception::exception;
    public:
    virtual const char* 
        what() 
        const throw();
};
#endif /* FUNCION_MATRIZ_HPP */
