#include "../macros.hpp"
#include "matrix_like/funcion_matriz.hpp"
#include "matrix_like/matrix.hpp"
#include "matrix_like/print.cpp"
#include <algorithm>

#define SIZE 5

const double h = 2.0 / (SIZE + 1);

struct eigen_index {
  int index;
  double val;
};

int main(int argc, char const *argv[]) {
  // auto f = std::ofstream("jacobi_out.txt");
  // macros::out = &f;
  auto mtx = mymtx::matrix::tridiag(SIZE, 1, -2, 1);
  mtx *= -1;
  mymtx::matrix eigvec(mtx.shape_x, mtx.shape_x);
  mymtx::vector eigval(mtx.shape_x);
  ANNOUNCE_TEST("jacobi eigen")
  measure_time(jacobi_eigen(mtx, eigval, eigvec, 10000);)
      eigen_index indexs[mtx.shape_x];
  for (int i = 0; i < mtx.shape_x; i++) {
    indexs[i].index = i;
    indexs[i].val = eigval[i];
  }
  eigvec = mymtx::matrix::traspose(eigvec);
  std::sort(indexs, indexs + mtx.shape_x,
            [](const eigen_index &a, const eigen_index &b) {
              return std::abs(a.val) < std::abs(b.val);
            });
  strm_out("min value:" << eigval[indexs[0].index]);
  strm_out("max value:" << eigval[indexs[mtx.shape_x - 1].index]);
  mymtx::vector::fwrite("min.vec", eigvec[indexs[0].index]);
  mymtx::vector::fwrite("max.vec", eigvec[indexs[mtx.shape_x - 1].index]);
  if (mtx.shape_x < 10) {
    printf("values\n");
    printv(eigval);
    printf("vectors\n");
    printm(mymtx::MatrixTraspose(eigvec));
  }
  // strm_out("");
  // printm(mtx);
  return 0;
}
