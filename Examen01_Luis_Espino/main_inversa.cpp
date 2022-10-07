#include "matrix/funcion_matriz.hpp"
#include "matrix/real_matrix.hpp"
#include <cmath>
#include <cstdio>
#include <ostream>

#include "print_time.hpp"

using namespace mymtx;
#include <iostream>

void printm(const mymtx::RealMatrix &m, std::ostream &out) {
  for (size_t i = 0; i < m.shape_y; ++i) {
    for (auto j = m.begin(i); j < m.end(i); ++j) {
      out << *j << " ";
    }
    out << std::endl;
  }
}
double norm(const RealMatrix &m);
const int size = 5000;
int main(int argc, const char **argv) {
  // generar matriz pendiagonal
  double coef[] = {-4, -8, 40, -8, -4};
  auto A = RealMatrix::pendiag(size, coef);
  RealMatrix L(A.shape_y, A.shape_y);
  // generar vector solucion
  RealMatrix I(size, size);
  auto start = high_resolution_clock::now();

  factor_cholesky_as_band(A, L, 2);
#pragma omp parallel for
  for (size_t i = 0; i < A.shape_y; i++) {
    RealVector b(size);
    RealVector a(size);
    b = 0;
    b[i] = 1;
    solve_cholesky(L, a, b);
    I.column(i) = a;
  }
  auto end = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(end - start);
  (*(macros::out)) << "time: " << duration.count() << "micro s" << std::endl;
  RealMatrix::fwrite("A_inversa.mtx", I);
  // Error
  printf("Error de inversa:%lf\n", norm(A * I - RealMatrix::identity(size)));
  return 0;
}
double norm(const RealMatrix &m) {
  double sum = 0;
  for (int i = 0; i < m.shape_y; ++i) {
    for (auto j = m.begin(i); j < m.end(i); ++j) {
      sum += *j * *j;
    }
  }
  return std::sqrt(sum);
}