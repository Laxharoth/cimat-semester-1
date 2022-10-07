#ifndef FUNCION_MATRIZ_HPP
#define FUNCION_MATRIZ_HPP

#include "matrix.hpp"

#include <cmath>
#include <cstdlib>
#include <future>
#include <random>
#include <vector>

struct eigen {
  mymtx::vector vect;
  double val;
};

void solucion_diagonal(const mymtx::matrix &matriz, mymtx::vector &incognitas,
                       mymtx::vector &result);
double determinante_diagonal(mymtx::matrix &matriz_diagonal);
void inversa_diagonal(mymtx::matrix &matriz_diagonal, mymtx::vector &inversa);
double determinante_triangular(mymtx::matrix &matriz_triangular);
void solucion_triangular_inf(mymtx::matrix &matriz, mymtx::vector &incognitas,
                             mymtx::vector &result, bool);
void solucion_triangular_inf(mymtx::matrix &matriz, mymtx::vector &incognitas,
                             mymtx::vector &result);
void solucion_triangular_sup(mymtx::matrix &matriz, mymtx::vector &incognitas,
                             mymtx::vector &result, bool);
void solucion_triangular_sup(mymtx::matrix &matriz, mymtx::vector &incognitas,
                             mymtx::vector &result);
void gauss(mymtx::matrix &matriz, mymtx::vector &variables,
           mymtx::vector &resultados);
void solucion_LDU(mymtx::matrix &matriz, mymtx::vector &incognitas,
                  mymtx::vector &result);
void solucion_crout(mymtx::matrix &matriz, mymtx::vector &incognitas,
                    mymtx::vector &result);
void solucion_doolittle(mymtx::matrix &matriz, mymtx::vector &incognitas,
                        mymtx::vector &result);
void solve_gauss_seidel(mymtx::matrix &matix, mymtx::vector &variables,
                        mymtx::vector &solutions, double *error);
void solve_cholesky(mymtx::matrix &cholesky_factored, mymtx::vector &variables,
                    mymtx::vector &solutions);
double normalize(mymtx::vector &vec);
void randomize(mymtx::vector &vec);
void power_iteration(const mymtx::matrix &A, mymtx::vector &V1,
                     const double tolerance, double &value, size_t n_values,
                     mymtx::matrix *_vec_holder, mymtx::vector *_val_holder,
                     const size_t max_iter);
void inverse_power_iteration(const mymtx::matrix &A, mymtx::vector &V1,
                             const double tolerance, double &value,
                             size_t n_values, mymtx::matrix *_vec_holder,
                             mymtx::vector *_val_holder, const size_t);
void jacobi_eigen(mymtx::matrix &A, mymtx::vector &e, mymtx::matrix &U,
                  const unsigned max_iter);
void subspace_pow(const mymtx::matrix &A, mymtx::matrix &I, mymtx::vector &eig);
void subspace_ipow(const mymtx::matrix &A, mymtx::matrix &I_t,
                   mymtx::vector &eig);
void rayleigh_method(const mymtx::matrix &A, mymtx::vector &V1, double &val);
void solve_qr(const mymtx::matrix &Q, mymtx::matrix &R, mymtx::vector &var,
              const mymtx::vector &res);
void conjugate_gradient(const mymtx::matrix &A, mymtx::vector &var,
                        mymtx::vector &result);
void conjugate_gradient_jacobi(const mymtx::matrix &A, mymtx::vector &var,
                               mymtx::vector &result);

// factorizacion

void crout(mymtx::matrix &matriz, mymtx::matrix &matriz_inferior,
           mymtx::matrix &matriz_superior);
void crout_as_band(mymtx::matrix &A_mtx, mymtx::matrix &L_mtx,
                   mymtx::matrix &U_mtx, int height, int width);
void crout_tridiagonal(mymtx::matrix &matriz, mymtx::matrix &matriz_inferior,
                       mymtx::matrix &matriz_superior);
void doolittle(mymtx::matrix &matriz, mymtx::matrix &matriz_inferior,
               mymtx::matrix &matriz_superior);
void doolittle_as_band(mymtx::matrix &A_mtx, mymtx::matrix &L_mtx,
                       mymtx::matrix &U_mtx, int height, int width);
void doolittle_tridiagonal(mymtx::matrix &matriz,
                           mymtx::matrix &matriz_inferior,
                           mymtx::matrix &matriz_superior);
void LDU_factor(mymtx::matrix &matriz, mymtx::matrix &matriz_inferior,
                mymtx::matrix &matriz_diagonal, mymtx::matrix &matriz_superior);
void LDU_factor_tridiagonal(mymtx::matrix &matriz,
                            mymtx::matrix &matriz_inferior,
                            mymtx::matrix &matriz_diagonal,
                            mymtx::matrix &matriz_superior);
void factor_cholesky(const mymtx::matrix &matix, mymtx::matrix &triangular);
void factor_cholesky_as_band(const mymtx::matrix &matix,
                             mymtx::matrix &triangular, int heigh);
void factor_cholesky_tridiag(mymtx::matrix &matix, mymtx::matrix &triangular);
void qr_decomposition(const mymtx::matrix &A, mymtx::matrix &Q,
                      mymtx::matrix &R);

class cant_factor_exception : public std::exception {
  using std::exception::exception;

public:
  virtual const char *what() const throw();
};
#endif /* FUNCION_MATRIZ_HPP */
