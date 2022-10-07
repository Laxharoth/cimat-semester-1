#include <fstream>
#include <iostream>

#include "functions.hpp"
#include "macros.hpp"
#include "matrix_like/funcion_matriz.hpp"
#include "matrix_like/matrix.hpp"

void fill_matrix(mymtx::matrix &matrix);
void fill_vector(mymtx::vector &vector, double Q, double q0, double qn);
void fill_x(mymtx::vector &vector, double min, double max);
void print_points(std::string filename, mymtx::vector &x, mymtx::vector &y);

int main(int argc, char const *argv[]) {
  double error;
  auto f = std::ofstream("out_gaus.txt");
  macros::out = &f;

  mymtx::matrix m10(10, 10);
  mymtx::vector v10(10);
  mymtx::vector s10(10);
  mymtx::vector x10(10);
  fill_matrix(m10);
  fill_vector(s10, 2, 0, 2);
  fill_x(x10, 0, 1);
  strm_out("Gauss Seidel 10:");
  measure_time(solve_gauss_seidel(m10, v10, s10, &error));
  strm_out("error:" << error);
  out_vector(v10, "solutionsi10.txt");
  print_points("10", x10, v10);

  mymtx::matrix m50(50, 50);
  mymtx::vector v50(50);
  mymtx::vector s50(50);
  mymtx::vector x50(50);
  fill_matrix(m50);
  fill_vector(s50, 2, 0, 2);
  fill_x(x50, 0, 1);
  strm_out("Gauss Seidel 50:");
  measure_time(solve_gauss_seidel(m50, v50, s50, &error));
  strm_out("error:" << error);

  out_vector(v50, "solutionsi50.txt");
  print_points("50", x50, v50);

  mymtx::matrix m100(200, 200);
  mymtx::vector v100(200);
  mymtx::vector s100(200);
  mymtx::vector x100(200);
  fill_matrix(m100);
  fill_vector(s100, 2, 0, 2);
  fill_x(x100, 0, 1);
  strm_out("Gauss Seidel 100:");
  measure_time(solve_gauss_seidel(m100, v100, s100, &error));
  strm_out("error:" << error);

  out_vector(v100, "solutionsi100.txt");
  print_points("100", x100, v100);

  mymtx::matrix m1000(1000, 1000);
  mymtx::vector v1000(1000);
  mymtx::vector s1000(1000);
  mymtx::vector x1000(1000);
  fill_matrix(m1000);
  fill_vector(s1000, 2, 0, 2);
  fill_x(x1000, 0, 1);
  strm_out("Gauss Seidel 1000:");
  measure_time(solve_gauss_seidel(m1000, v1000, s1000, &error));
  strm_out("error:" << error);

  out_vector(v1000, "solutionsi1000.txt");
  print_points("1000", x1000, v1000);

  mymtx::matrix m10000(10000, 10000);
  mymtx::vector v10000(10000);
  mymtx::vector s10000(10000);
  mymtx::vector x10000(10000);
  fill_matrix(m10000);
  fill_vector(s10000, 2, 0, 2);
  fill_x(x10000, 0, 1);
  strm_out("Gauss Seidel 10000:");
  measure_time(solve_gauss_seidel(m10000, v10000, s10000, &error));
  strm_out("error:" << error);

  out_vector(v10000, "solutionsi10000.txt");
  print_points("10000", x10000, v10000);

  return 0;
}

void fill_matrix(mymtx::matrix &matrix) {
  double h = 1.0 / (matrix.shape_y + 1);
  for (size_t i = 0; i < matrix.shape_y; i++) {
    matrix[i][i] = 2;
    if (i != 0)
      matrix[i][i - 1] = -1;
    if (i + 1 < matrix.shape_x)
      matrix[i][i + 1] = -1;
  }
  matrix *= (1.0 / (h));
}
void fill_vector(mymtx::vector &vector, double Q, double q0, double qn) {
  double h = 1.0 / (vector.size + 1);
  vector[0] = -2 * h;
  vector[vector.size - 1] = -2 * h + 2 / h;
  for (size_t i = 0; i < vector.size - 1; i++) {
    vector[i] = -2 * h;
  }
}

void fill_x(mymtx::vector &vector, double min, double max) {
  double increment = (max - min) / (vector.size - 1);
  for (size_t i = 0; i < vector.size; i++) {
    vector[i] = min + i * increment;
  }
}

void print_points(std::string filename, mymtx::vector &x, mymtx::vector &y) {
  auto filey = std::ofstream(("p" + filename + "yg.py").c_str());
  auto filex = std::ofstream(("p" + filename + "x.py").c_str());

  filey << "import numpy as np" << std::endl;
  filex << "import numpy as np" << std::endl;

  filey << "points" << filename << "y=np.array([";
  filex << "points" << filename << "x=np.array([";
  for (auto i = 0; i < x.size; ++i) {
    filex << x[i] << ",";
  }
  filey << 0 << ",";
  for (auto i = 0; i < y.size; ++i) {
    filey << y[i] << ",";
  }
  filey << 2 << ",";
  filey << "])";
  filex << "])";
  filey.close();
  filex.close();
}