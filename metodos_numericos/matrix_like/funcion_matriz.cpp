#include "funcion_matriz.hpp"
#include "factorizacion.cpp"
#include "matrix.hpp"
#include <cstdio>
#define PI 3.141592653589793
#define JACOBI_UMBRAL 10e-4
#define ZERO_UMBRAL 10e-5

class randgen {
  std::random_device *rd;
  std::mt19937 *gen;
  std::uniform_real_distribution<double> *dis;
  randgen() {
    rd = new std::random_device();
    gen = new std::mt19937((*rd)());
    dis = new std::uniform_real_distribution<double>(0, 50);
  }
  ~randgen() {
    delete rd;
    delete gen;
    delete dis;
  }

public:
  double generate() { return (*dis)(*gen); }
  static randgen &get_randgen() {
    static randgen rng;
    return rng;
  }
};

void solucion_diagonal(const mymtx::matrix &matriz, mymtx::vector &incognitas,
                       mymtx::vector &result) {
  const size_t size = matriz.shape_x;
  for (int i = 0; i < size; ++i) {
    incognitas[i] = result[i] / matriz[i][i];
  }
}
double determinante_diagonal(mymtx::matrix &matriz_diagonal) {
  const size_t size = matriz_diagonal.shape_x;
  double result{1};
  for (int i = 0; i < size; ++i)
    result *= matriz_diagonal[i][i];
  return result;
}
void inversa_diagonal(mymtx::matrix &matriz_diagonal, mymtx::matrix &inversa) {
  const size_t size = matriz_diagonal.shape_x;
  for (int i = 0; i < size; ++i) {
    inversa[i][i] = 1 / matriz_diagonal[i][i];
  }
}
double determinante_triangular(mymtx::matrix &matriz_triangular) {
  const size_t size = matriz_triangular.shape_x;
  double result{1};
  for (int i = 0; i < size; ++i)
    result *= matriz_triangular[i][i];
  return result;
}
void solucion_triangular_inf(const mymtx::matrix &matriz,
                             mymtx::vector &incognitas,
                             const mymtx::vector &result) {
  solucion_triangular_inf(matriz, incognitas, result, true);
}
void solucion_triangular_inf(const mymtx::matrix &matriz,
                             mymtx::vector &incognitas,
                             const mymtx::vector &result,
                             bool compute_diagonale) {
  const size_t size = matriz.shape_x;
  for (int i = 0; i < size; ++i) {
    incognitas[i] = result[i];
    auto iter = matriz[i].begin();
    for (int j = 0; j <= i - 1; ++j) {
      incognitas[i] -= matriz[i][j] * incognitas[j];
    }
    if (compute_diagonale)
      incognitas[i] /= matriz[i][i];
  }
}
void solucion_triangular_inf_as_band(mymtx::matrix &matriz,
                                     mymtx::vector &incognitas,
                                     mymtx::vector &result,
                                     bool compute_diagonale, size_t heigh) {
  const size_t size = matriz.shape_x;
  for (int i = 0; i < size; ++i) {
    incognitas[i] = result[i];
    auto iter = matriz[i].begin();
    size_t start = 0;
    if (i > heigh)
      start = i - heigh;
    for (int j = start; j <= i - 1; ++j) {
      incognitas[i] -= matriz[i][j] * incognitas[j];
    }
    if (compute_diagonale)
      incognitas[i] /= matriz[i][i];
  }
}
void solucion_triangular_sup(const mymtx::matrix &matriz,
                             mymtx::vector &incognitas,
                             const mymtx::vector &result) {
  solucion_triangular_sup(matriz, incognitas, result, true);
}
void solucion_triangular_sup(const mymtx::matrix &matriz,
                             mymtx::vector &incognitas,
                             const mymtx::vector &result,
                             bool compute_diagonale) {
  const size_t size = matriz.shape_x;
  for (int i = size - 1; i >= 0; --i) {
    incognitas[i] = result[i];
    for (int j = i + 1; j < size; ++j) {
      incognitas[i] -= matriz[i][j] * incognitas[j];
    }
    if (compute_diagonale)
      incognitas[i] /= matriz[i][i];
  }
}
void solucion_triangular_sup_as_band(mymtx::matrix &matriz,
                                     mymtx::vector &incognitas,
                                     mymtx::vector &result,
                                     bool compute_diagonale, size_t width) {
  const size_t size = matriz.shape_x;
  for (int i = size - 1; i >= 0; --i) {
    incognitas[i] = result[i];
    size_t end = size;
    if (i < width)
      end = width - i;
    for (int j = i + 1; j < size; ++j) {
      incognitas[i] -= matriz[i][j] * incognitas[j];
    }
    if (compute_diagonale)
      incognitas[i] /= matriz[i][i];
  }
}
void gauss(mymtx::matrix &matriz, mymtx::vector &variables,
           mymtx::vector &resultados) {
  const size_t size = matriz.shape_x;
  for (int i = 0; i < size; ++i) {
    const double divide_privote = matriz[i][i];
    if (std::abs(divide_privote) < ZERO_UMBRAL)
      throw std::runtime_error("zero division");
    for (auto j = matriz[i].begin() + i; j < matriz[i].end(); ++j) {
      *j /= divide_privote;
    }
    resultados[i] /= divide_privote;
    for (int j = i + 1; j < size; ++j) {
      auto iter_pivote = matriz[i].begin() + (i + 1);
      auto iter = matriz[j].begin() + i;
      auto end = matriz[j].end();
      const double coeficiente_eliminar = *iter;
      *iter = 0;
      for (++iter; iter < end; ++iter, ++iter_pivote) {
        *iter -= coeficiente_eliminar * (*iter_pivote);
      }
      resultados[j] -= coeficiente_eliminar * resultados[i] / matriz[i][i];
    }
  }
  solucion_triangular_sup(matriz, variables, resultados);
}
void solucion_LDU(mymtx::matrix &matriz, mymtx::vector &incognitas,
                  mymtx::vector &result) {
  mymtx::vector aux1(incognitas.size);
  mymtx::vector aux2(incognitas.size);
  solucion_triangular_inf(matriz, aux1, result, false);
  solucion_diagonal(matriz, aux2, aux1);
  solucion_triangular_sup(matriz, incognitas, aux2, false);
}
void solucion_crout(const mymtx::matrix &matriz, mymtx::vector &incognitas,
                    const mymtx::vector &result) {
  mymtx::vector aux1(incognitas.size);
  solucion_triangular_inf(matriz, aux1, result);
  solucion_triangular_sup(matriz, incognitas, aux1, false);
}
void solve_crout_as_band(mymtx::matrix &matriz, mymtx::vector &incognitas,
                         mymtx::vector &result, const size_t heigh,
                         const size_t width) {
  mymtx::vector aux1(incognitas.size);
  solucion_triangular_inf_as_band(matriz, aux1, result, true, heigh);
  solucion_triangular_sup_as_band(matriz, incognitas, aux1, false, width);
}
void solucion_doolittle(mymtx::matrix &matriz, mymtx::vector &incognitas,
                        mymtx::vector &result) {
  mymtx::vector aux1(incognitas.size);
  solucion_triangular_inf(matriz, aux1, result, false);
  solucion_triangular_sup(matriz, incognitas, aux1);
}
double normalize(mymtx::vector &vec) {
  double sum{0};
  for (auto i = vec.begin(); i != vec.end(); ++i) {
    sum += *i * *i;
  }
  sum = std::sqrt(sum);
  for (auto i = vec.begin(); i != vec.end(); ++i) {
    *i /= sum;
  }
  return sum;
}
void normalize(mymtx::matrix &mtx) {
  double norm1{0}, norm2{0};
  for (size_t i = 0; i < mtx.shape_y; ++i) {
    norm1 += mtx[i][0] * mtx[i][0];
  }
  norm1 = std::sqrt(norm1);
  for (size_t k = 1; k < mtx.shape_x; ++k) {
    for (size_t i = 0; i < mtx.shape_y; ++i) {
      norm2 += mtx[i][k] * mtx[i][k];
      mtx[i][k - 1] /= norm1;
    }
    norm1 = std::sqrt(norm2);
    norm2 = 0;
  }
  for (size_t i = 0; i < mtx.shape_y; ++i) {
    norm1 += mtx[i][mtx.shape_x - 1] * mtx[i][mtx.shape_x - 1];
  }
}
void randomize(mymtx::vector &vec) {
  randgen &rng = randgen::get_randgen();
  for (auto i = vec.begin(); i != vec.end(); ++i) {
    *i = rng.generate();
  }
}
void randomize(mymtx::matrix &mtx) {
  randgen &rng = randgen::get_randgen();
  for (size_t i = 0; i < mtx.shape_y; ++i) {
    for (auto j = mtx.begin(i); j < mtx.end(i); ++j) {
      *j = rng.generate();
    }
  }
}
void power_iteration(const mymtx::matrix &A, mymtx::vector &V1,
                     const double tolerance, double &value, size_t n_values,
                     mymtx::matrix *_vec_holder, mymtx::vector *_val_holder,
                     const size_t max_iter) {
  mymtx::vector V0(V1.size);
  auto &vec_holder = *_vec_holder;
  auto &val_holder = *_val_holder;
  double error = 1E20;
  double old_val{0};
  size_t found{0}, iter{0};
  for (size_t k = 0; k < n_values; ++k) {
    iter = 0;
    randomize(V0);
    normalize(V0);
    while (error > tolerance && iter < max_iter) {
      iter++;
      for (size_t i = 0; i < found; i++) {
        V0 -= vec_holder[i] * (V0 * vec_holder[i]);
      }
      V1 = A.prod_as_band(V0, 2, 2);
      value = V0 * V1;
      error = std::abs(old_val - value);
      old_val = value;
      normalize(V1);
      V0 = V1;
    }
    if (_vec_holder == nullptr || _val_holder == nullptr || k >= n_values ||
        iter >= max_iter) {
      return;
    }
    vec_holder[found] = V1;
    val_holder[found] = value;
    found++;
    error = 1;
  }
}
void inverse_power_iteration(const mymtx::matrix &A, mymtx::vector &V1,
                             const double tolerance, double &value,
                             size_t n_values, mymtx::matrix *_vec_holder,
                             mymtx::vector *_val_holder,
                             const size_t max_iter) {
  mymtx::vector V0(V1.size);
  auto &vec_holder = *_vec_holder;
  auto &val_holder = *_val_holder;
  double error = 1E20;
  double old_val{0};
  size_t found{0};
  unsigned iter = 0;
  mymtx::matrix working_m = A;
  // crout(working_m,working_m,working_m);
  factor_cholesky_as_band(A, working_m, 2);
  for (size_t k = 0; k < n_values; ++k) {
    iter = 0;
    randomize(V0);
    normalize(V0);
    while (error > tolerance && iter < max_iter) {
      iter++;
      for (size_t i = 0; i < found; i++) {
        V0 -= vec_holder[i] * (V0 * vec_holder[i]);
      }
      solve_cholesky(working_m, V1, V0);
      value = 1 / (V0 * V1);
      error = std::abs(old_val - value);
      old_val = value;
      normalize(V1);
      V0 = V1;
    }
    if (_vec_holder == nullptr || _val_holder == nullptr || k >= n_values ||
        iter >= max_iter) {
      return;
    }
    vec_holder[found] = V1;
    val_holder[found] = value;
    found++;
    error = 1;
  }
}
void solve_cholesky(const mymtx::matrix &cholesky_factored,
                    mymtx::vector &variables, const mymtx::vector &solutions) {
  mymtx::vector tmp(solutions.size);
  solucion_triangular_inf(cholesky_factored, tmp, solutions);
  solucion_triangular_sup(cholesky_factored, variables, tmp);
}
void jacobi_eigen(mymtx::matrix &A, mymtx::vector &e, mymtx::matrix &U,
                  const unsigned max_iter) {
  size_t col, row;
  const size_t n = A.shape_y;
  unsigned iter = 0;
  mymtx::matrix B = A;
  const mymtx::matrix I = mymtx::matrix::identity(n);
  mymtx::matrix R(n, n);
  mymtx::vector new_row_i(n);
  mymtx::vector new_row_j(n);
  auto IndexOfMax = [](const mymtx::matrix &inA, size_t &incol, size_t &inrow) {
    double max = inrow = 0;
    incol = 1;
    for (size_t k = 0; k < inA.shape_y; ++k)
      for (auto i = inA.begin(k); i < inA.end(k); ++i) {
        if (i.get_col() == i.get_row())
          continue;
        if (std::abs(*i) > std::abs(max)) {
          max = *i;
          incol = i.get_col();
          inrow = i.get_row();
        }
      }
  };
  auto rotate = [&](mymtx::matrix &U_mtx, const size_t row_lmb,
                    const size_t col_lmb) {
    double tan_2, tan_, cos_, sin_, theta;
    const size_t n = B.shape_y;
    if (row_lmb == col_lmb)
      return;
    if (B[row_lmb][row_lmb] != B[col_lmb][col_lmb]) {
      tan_2 = (2 * B[row_lmb][col_lmb]) /
              (B[row_lmb][row_lmb] - B[col_lmb][col_lmb]);
      tan_ = tan_2 * tan_2 / (1 + std::sqrt(1 + tan_2 * tan_2));
      cos_ = 1 / std::sqrt(tan_ * tan_ + 1);
      sin_ = cos_ * tan_;
    } else {
      cos_ = std::cos(PI / 4);
      sin_ = std::sin(PI / 4);
    }
    R = I;
    R[row_lmb][row_lmb] = R[col_lmb][col_lmb] = cos_;
    R[row_lmb][col_lmb] = sin_;
    R[col_lmb][row_lmb] = -sin_;
    U_mtx *= R;
    B = mymtx::MatrixTraspose(U_mtx) * A.prod_as_band(U_mtx, 1, 1);
  };
  U = mymtx::matrix::identity(n);
  while (iter++ < max_iter) {
    IndexOfMax(B, row, col);
    if (std::abs(B[row][col]) < JACOBI_UMBRAL) {
      break;
    }
    rotate(U, row, col);
    // printf("%d\n",iter);
  }
  for (size_t i = 0; i < n; ++i)
    e[i] = B[i][i];
}
void subspace_pow(const mymtx::matrix &A, mymtx::matrix &I_t,
                  mymtx::vector &eig) {
  size_t m = eig.size;
  mymtx::MatrixTraspose I(I_t);
  mymtx::matrix B(m, m);
  mymtx::matrix ro(m, m);
  mymtx::vector d(m);
  double val1;
  randomize(I_t);
  for (size_t i = 0; i < I_t.shape_y; i++) {
    I_t[i] /= I_t[i].distance();
  }
  size_t iter, row, col;
  double error = 1;
  auto eig_old = eig;
  auto IndexOfMax = [](const mymtx::matrix &inA) {
    double max = 0;
    for (size_t k = 0; k < inA.shape_y; ++k) {
      for (auto i = inA.begin(k); i < inA.end(k); ++i) {
        if (i.get_col() == i.get_row())
          continue;
        if (std::abs(*i) > std::abs(max)) {
          max = *i;
        }
      }
    }
    return max;
  };
  while (1) {
    iter = 0;
    for (size_t k = 0; k < m; ++k) {
      iter = 0;
      mymtx::vector V0 = I_t[k];
      while (iter < 40) {
        iter++;
        for (size_t i = 0; i < k; i++) {
          V0 -= I_t[i] * (I_t[i] * V0);
        }
        V0 = A.prod_as_band(V0, 2, 2);
        normalize(V0);
      }
      I_t[k] = V0;
    }
    B = I_t * A * I;
    eig_old = eig;
    jacobi_eigen(B, eig, ro, 10000);
    if ((eig - eig_old).distance() < ZERO_UMBRAL)
      break;
    I_t = mymtx::MatrixTraspose(ro) * I_t;
    for (size_t i = 0; i < I_t.shape_y; i++) {
      I_t[i] /= I_t[i].distance();
    }
  }
  for (size_t i = 0; i < I_t.shape_y; i++) {
    I_t[i] /= I_t[i].distance();
  }
}
void subspace_ipow(const mymtx::matrix &A, mymtx::matrix &I_t,
                   mymtx::vector &eig) {
  size_t m = eig.size;
  mymtx::MatrixTraspose I(I_t);
  mymtx::matrix B(m, m);
  mymtx::matrix Chol(A.shape_y, A.shape_x);
  mymtx::matrix ro(m, m);
  mymtx::vector d(m);
  double val1;
  randomize(I_t);
  for (size_t i = 0; i < I_t.shape_y; i++) {
    I_t[i] /= I_t[i].distance();
  }
  size_t iter, row, col;
  double error = 1;
  auto eig_old = eig;
  auto IndexOfMax = [](const mymtx::matrix &inA) {
    double max = 0;
    for (size_t k = 0; k < inA.shape_y; ++k) {
      for (auto i = inA.begin(k); i < inA.end(k); ++i) {
        if (i.get_col() == i.get_row())
          continue;
        if (std::abs(*i) > std::abs(max)) {
          max = *i;
        }
      }
    }
    return max;
  };
  factor_cholesky(A, Chol);
  while (1) {
    iter = 0;
    for (size_t k = 0; k < m; ++k) {
      iter = 0;
      mymtx::vector V0 = I_t[k];
      mymtx::vector V1(V0.size);
      while (iter < 40) {
        iter++;
        for (size_t i = 0; i < k; i++) {
          V0 -= I_t[i] * (I_t[i] * V0);
        }
        solve_cholesky(Chol, V1, V0);
        V0 = V1;
        normalize(V0);
      }
      I_t[k] = V0;
    }
    B = I_t * A * I;
    eig_old = eig;
    jacobi_eigen(B, eig, ro, 100000);
    if ((eig - eig_old).distance() < ZERO_UMBRAL)
      break;
    I_t = mymtx::MatrixTraspose(ro) * I_t;
    for (size_t i = 0; i < I_t.shape_y; i++) {
      I_t[i] /= I_t[i].distance();
    }
  }
  for (size_t i = 0; i < I_t.shape_y; i++) {
    I_t[i] /= I_t[i].distance();
  }
}
void rayleigh_method(const mymtx::matrix &A, mymtx::vector &V1, double &val) {
  const auto identity = mymtx::matrix::identity(V1.size);
  if (((A - identity * val) * V1).distance() < ZERO_UMBRAL)
    return;
  normalize(V1);
  mymtx::vector V0 = V1;
  val = V0 * (A * V0);
  mymtx::matrix B = (A - identity * val);
  mymtx::matrix Q(B.shape_y, B.shape_x), R(B.shape_y, B.shape_x);
  double old_val = val + 1;
  while ((B * V1).distance() > ZERO_UMBRAL &&
         std::abs(old_val - val) > ZERO_UMBRAL) {
    try {
      // crout(B,B,B);
      // solucion_crout(B,V1,V0);
      gauss(B, V1, V0);
      // qr_decomposition(B,Q,R);
      // solve_qr(Q,R,V1,V0);
    } catch (...) {
      return;
    }
    normalize(V1);
    V0 = V1;
    old_val = val;
    val = V0 * (A * V0);
    B = (A - identity * val);
  }
}
void solve_qr(const mymtx::matrix &Q, mymtx::matrix &R, mymtx::vector &var,
              const mymtx::vector &res) {
  mymtx::vector Qres = mymtx::MatrixTraspose(Q) * res;
  solucion_triangular_sup(R, var, Qres);
}
void conjugate_gradient(const mymtx::matrix &A, mymtx::vector &x,
                        mymtx::vector &b) {
  const int n = A.shape_y;
  mymtx::vector r = b;
  mymtx::vector p = r;
  int k = 0;

  mymtx::vector Rold(r);
  while (k < n) {
    Rold = r;
    mymtx::vector AP = A * p;
    double alpha = r * r / std::max(p * AP, ZERO_UMBRAL);
    x = x + alpha * p;
    r = r - alpha * AP;

    if (r.distance() < ZERO_UMBRAL)
      break; // Convergence test

    double beta = r * r / std::max(Rold * Rold, ZERO_UMBRAL);
    p = r + beta * p;
    k++;
  }
}
void conjugate_gradient_jacobi(const mymtx::matrix &A, mymtx::vector &x,
                               mymtx::vector &b) {
  double alpha, betha;
  auto r = b;
  mymtx::vector z(x.size);
  mymtx::vector w(x.size);
  solucion_diagonal(A, z, r);
  auto p = z;
  unsigned i{0};
  double rzold;
  while (i++ < x.size * 3) {
    rzold = r * z;
    w = A * p;
    alpha = (rzold) / (p * w);
    x = x + alpha * p;
    r = r - alpha * w;
    if (std::sqrt(r * r) < ZERO_UMBRAL)
      break;
    solucion_diagonal(A, z, r);
    betha = (r * z) / (rzold);
    p = r + (betha)*p;
  }
}
