#include "macros.hpp"
#include "matrix_like/funcion_matriz.hpp"
#include "matrix_like/matrix.hpp"

#define PI 3.141592653589793

const double h = 2.0 / (1001);
auto phi_i_j = [](const int i, const int j) {
  return std::sin(i * j * h * PI);
};

int main(int argc, char const *argv[]) {
  auto mtx = mymtx::matrix::tridiag(1000, 1, -2, 1);
  auto V0rig = mymtx::vector::normal(1000);
  auto V0 = V0rig;
  mymtx::vector V1(1000);
  double val = 0.0;
  ANNOUNCE_TEST("1 valor");
  strm_out("power iteration") measure_time(
      power_iteration(mtx, V1, 1E-20, val, 2, nullptr, nullptr, 10000););

  auto biggest_l = 3.9546266935713628;
  mymtx::vector biggest_v(1000);
  for (auto i = biggest_v.begin(); i != biggest_v.end(); ++i)
    *i = phi_i_j(1, i.get_col());
  normalize(biggest_v);
  normalize(V1);

  strm_out("error biggest eigen");
  strm_out("value");
  strm_out(std::abs(std::abs(val) - std::abs(biggest_l))) strm_out("vector");
  int i = 0;
  double error = 0;
  double contrl = 0;
  for (auto i = V1.begin(), j = biggest_v.begin(); i != V1.end(); ++i, ++j) {
    error += std::abs((*i) - (*j));
    contrl += std::abs(*j);
  }
  strm_out(error / contrl);

  V0 = V0rig;
  strm_out("inverse power iteration")
      measure_time(inverse_power_iteration(mtx, V1, 1E-20, val, 1, nullptr,
                                           nullptr, 10000););

  auto smallest = 0.1742732695496015;
  mymtx::vector smallest_v(1);

  strm_out("error smallest eigen");
  strm_out("value");
  strm_out(std::abs(std::abs(val) - std::abs(smallest))) strm_out("vector");
  error = 0;
  contrl = 0;
  for (auto i = V1.begin(), j = smallest_v.begin(); i != V1.end(); ++i, ++j) {
    error += std::abs((*i) - (*j));
    contrl += std::abs(*j);
  }
  strm_out(error / contrl);
  return 0;
}
