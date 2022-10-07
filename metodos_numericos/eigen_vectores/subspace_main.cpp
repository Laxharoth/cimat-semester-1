#include "macros.hpp"
#include "matrix_like/funcion_matriz.hpp"
#include "matrix_like/matrix.hpp"
#include "matrix_like/print.cpp"
#include <algorithm>
#include <fstream>
#include <string>

#define SIZE 5

const double h = 2.0 / (SIZE + 1);

struct eigen_index {
  int index;
  double val;
};
mymtx::matrix readMtxtxt(const std::string filename);
mymtx::vector readVectxt(const std::string filename);

int main(int argc, char const *argv[]) {
  // auto f = std::ofstream("subspace_out.txt");
  // macros::out = &f;
  /*
  ANNOUNCE_TEST("Subespace Power Iter 3x3"){
  mymtx::matrix mtx = readMtxtxt("vector_in/Eigen_3x3.txt");
  const unsigned amt = 2;
  mymtx::matrix I_t(amt, mtx.shape_y);
  mymtx::vector eig(amt);
  measure_time(
      subspace_pow(mtx, I_t ,eig);
  )
  strm_out("vectors");
  printm( I_t, *macros::out );
  strm_out("values");
  printv(eig, *macros::out);
  mymtx::matrix::fwrite("vector.out/actual.pow.Eigen_3x3.mtx.txt",I_t);
  mymtx::vector::fwrite("vector.out/actual.pow.Eigen_3x3.vec.txt",eig);
  }
  ANNOUNCE_TEST("Subespace Inverse Power 3x3"){
  mymtx::matrix mtx = readMtxtxt("vector_in/Eigen_3x3.txt");
  const unsigned amt = 2;
  mymtx::matrix I_t(amt, mtx.shape_y);
  mymtx::vector eig(amt);
  measure_time(
      subspace_ipow(mtx, I_t ,eig);
  )
  strm_out("vectors");
  printm( I_t,*macros::out);
  strm_out("values");
  printv(eig,*macros::out);
  mymtx::matrix::fwrite("vector.out/actual.ipow.Eigen_3x3.mtx.txt",I_t);
  mymtx::vector::fwrite("vector.out/actual.ipow.Eigen_3x3.vec.txt",eig);
  }
  ANNOUNCE_TEST("Rayleigh 3x3"){
  mymtx::matrix mtx = readMtxtxt("vector_in/Eigen_3x3.txt");
  mymtx::vector V0{-0.2,-0.9,0.9};
  double val = 9.9;
  measure_time(
      rayleigh_method(mtx, V0, val);
  )
  strm_out("val" << std::endl << val);
  strm_out("vec");
  mymtx::vector valv(1);
  valv[0] = val;
  printv(V0,*macros::out);
  mymtx::vector::fwrite("vector.out/actual.ray.Eigen_vec_3x3.vec.txt",V0);
  mymtx::vector::fwrite("vector.out/actual.ray.Eigen_val_3x3.vec.txt",valv);
  }
  ANNOUNCE_TEST("QR 3x3"){
      mymtx::matrix mtx = readMtxtxt("vector_in/M_sys_3x3.txt");
      mymtx::vector V0 = readVectxt("vector_in/V_sys_3x1.txt");
      mymtx::matrix Q(mtx.shape_y, mtx.shape_x);
      mymtx::matrix R(mtx.shape_y, mtx.shape_x);
      mymtx::vector V1(V0.size);
      measure_time(
          qr_decomposition(mtx,Q,R);
          solve_qr(Q,R,V1,V0);
      )
      strm_out("solution");
      printv(V1,*macros::out);
      mymtx::vector::fwrite("vector.out/actual.QR.V_sol_3x1.vec.txt",V1);
  }
  ANNOUNCE_TEST("Gradient 3x3"){
      mymtx::matrix mtx = readMtxtxt("vector_in/M_sys_3x3.txt");
      mymtx::vector V0 = readVectxt("vector_in/V_sys_3x1.txt");
      mymtx::vector V1(V0.size);
      measure_time(
          conjugate_gradient(mtx,V1,V0);
      )
      strm_out("solution");
      printv(V1,*macros::out);
      mymtx::vector::fwrite("vector.out/actual.Grad.V_sol_3x1.vec.txt",V1);
  }
  ANNOUNCE_TEST("Gradient Jacobi 3x3"){
      mymtx::matrix mtx = readMtxtxt("vector_in/M_sys_3x3.txt");
      mymtx::vector V0 = readVectxt("vector_in/V_sys_3x1.txt");
      mymtx::vector V1(V0.size);
      measure_time(
          conjugate_gradient_jacobi(mtx,V1,V0);
      )
      strm_out("solution");
      printv(V1,*macros::out);
      mymtx::vector::fwrite("vector.out/actual.GradJ.V_sol_3x1.vec.txt",V1);
  }
  ANNOUNCE_TEST("Subespace Power Iter 50x50"){
  mymtx::matrix mtx = readMtxtxt("vector_in/Eigen_50x50.txt");
  const unsigned amt = 2;
  mymtx::matrix I_t(amt, mtx.shape_y);
  mymtx::vector eig(amt);
  measure_time(
      subspace_pow(mtx, I_t ,eig);
  )
  strm_out("vectors");
  printm( I_t, *macros::out );
  strm_out("values");
  printv(eig, *macros::out);
  mymtx::matrix::fwrite("vector.out/actual.pow.Eigen_50x50.mtx.txt",I_t);
  mymtx::vector::fwrite("vector.out/actual.pow.Eigen_50x50.vec.txt",eig);
  }
  ANNOUNCE_TEST("Subespace Inverse Power 50x50"){
  mymtx::matrix mtx = readMtxtxt("vector_in/Eigen_50x50.txt");
  const unsigned amt = 2;
  mymtx::matrix I_t(amt, mtx.shape_y);
  mymtx::vector eig(amt);
  measure_time(
      subspace_ipow(mtx, I_t ,eig);
  )
  strm_out("vectors");
  printm( I_t,*macros::out);
  strm_out("values");
  printv(eig,*macros::out);
  mymtx::matrix::fwrite("vector.out/actual.ipow.Eigen_50x50.mtx.txt",I_t);
  mymtx::vector::fwrite("vector.out/actual.ipow.Eigen_50x50.vec.txt",eig);
  }
  ANNOUNCE_TEST("Rayleigh 50x50"){
  mymtx::matrix mtx = readMtxtxt("vector_in/Eigen_50x50.txt");
  mymtx::vector
  V0{0.000351,1.67e-05,7.6446e-05,0.00046,4.0194e-05,9.6382e-05,0.00040,0.00030,0.00042,0.00013,0.00032,0.00069,0.00018,0.00094,0.00052,0.00070,9.0017e-05,0.00048,0.00159,0.0010,0.0018,0.00103,0.00067,0.0038,0.99,-0.0031,-0.00050,-0.0029,-0.00023,-0.00099,-0.0004,-0.0002,-0.0003,-0.00071,-0.00099,-0.00022,-0.00042,-0.00059,-0.00061,-9.4417e-05,-0.0001519,-0.00039747,-0.00037297,-0.00048534,-0.00015605,-0.00031083,-0.00026523,-0.00032853,-0.00039186,-0.00030112};
  double val = 230;
  measure_time(
      rayleigh_method(mtx, V0, val);
  )
  strm_out("val" << std::endl << val);
  strm_out("vec");
  mymtx::vector valv(1);
  valv[0] = val;
  printv(V0,*macros::out);
  mymtx::vector::fwrite("vector.out/actual.ray.Eigen_vec_50x50.vec.txt",V0);
  mymtx::vector::fwrite("vector.out/actual.ray.Eigen_val_50x50.vec.txt",valv);
  }
  ANNOUNCE_TEST("QR 125x125"){
      mymtx::matrix mtx = readMtxtxt("vector_in/M_sys_125x125.txt");
      mymtx::vector V0 = readVectxt("vector_in/V_sys_125x1.txt");
      mymtx::matrix Q(mtx.shape_y, mtx.shape_x);
      mymtx::matrix R(mtx.shape_y, mtx.shape_x);
      mymtx::vector V1(V0.size);
      measure_time(
          qr_decomposition(mtx,Q,R);
          solve_qr(Q,R,V1,V0);
      )
      strm_out("solution");
      printv(V1,*macros::out);
      mymtx::vector::fwrite("vector.out/actual.QR.V_sol_125x1.vec.txt",V1);
  }
  ANNOUNCE_TEST("Gradient 125x125"){
      mymtx::matrix mtx = readMtxtxt("vector_in/M_sys_125x125.txt");
      mymtx::vector V0 = readVectxt("vector_in/V_sys_125x1.txt");
      mymtx::vector V1(V0.size);
      measure_time(
          conjugate_gradient(mtx,V1,V0);
      )
      strm_out("solution");
      printv(V1,*macros::out);
      mymtx::vector::fwrite("vector.out/actual.Grad.V_sol_125x1.vec.txt",V1);
  }
  ANNOUNCE_TEST("Gradient Jacobi 125x125"){
      mymtx::matrix mtx = readMtxtxt("vector_in/M_sys_125x125.txt");
      mymtx::vector V0 = readVectxt("vector_in/V_sys_125x1.txt");
      mymtx::vector V1(V0.size);
      measure_time(
          conjugate_gradient_jacobi(mtx,V1,V0);
      )
      strm_out("solution");
      printv(V1,*macros::out);
      mymtx::vector::fwrite("vector.out/actual.GradJ.V_sol_125x1.vec.txt",V1);
  }
  // f.close();
  //*/
  mymtx::matrix A = mymtx::matrix::tridiag(10, -1, 2, -1);
  mymtx::matrix Q(A.shape_y, A.shape_x);
  mymtx::matrix R(A.shape_y, A.shape_x);
  qr_decomposition(A, Q, R);
  printm(Q);
  printm(R);
  printm(Q * R);
  return 0;
}
mymtx::matrix readMtxtxt(const std::string filename) {
  auto f_matriz = std::ifstream(filename.c_str());
  size_t rows;
  size_t cols;
  f_matriz >> rows;
  f_matriz >> cols;
  mymtx::matrix matriz(rows, cols);
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      f_matriz >> matriz[i][j];
    }
  }
  return matriz;
}
mymtx::vector readVectxt(const std::string filename) {
  auto f_matriz = std::ifstream(filename.c_str());
  size_t rows;
  size_t cols;
  f_matriz >> rows;
  f_matriz >> cols;
  mymtx::vector matriz(rows);
  for (size_t i = 0; i < rows; ++i) {
    f_matriz >> matriz[i];
  }
  return matriz;
}
