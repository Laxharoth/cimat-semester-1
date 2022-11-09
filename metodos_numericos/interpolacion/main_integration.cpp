#include "function_wrapper/function_wrapper.hpp"
#include "interpolation.hpp"
#include "macros.hpp"
#include "matrix_like/funcion_matriz.hpp"
#include "matrix_like/matrix.hpp"
#include "matrix_like/print.cpp"
#include <cmath>
#include <cstdio>

class exp : public FunctionWrapper {
  double eval(const double &x) const { return x * std::exp(x); }
  double eval(const double &x) { return x * std::exp(x); }
};

class polyn : public FunctionWrapper {
  double eval(const double &x) const {
    return (x * x + 3 * x + 3) * (2 * x * x - 9 * x + 7);
  }
  double eval(const double &x) {
    return (x * x + 3 * x + 3) * (2 * x * x - 9 * x + 7);
  }
};

int main(int argc, const char **argv) {
  class exp e;
  polyn p;

  const double from = 0;
  const double to = 5;
  FunctionWrapper *f = &p;
  {
    strm_out("polynomial: (x^2+3x+3)(2x^2-9x+7) from 0 to 5");
    const double integral_expected = 2735 / 12.0;
    const double integral_actual = gaussian_cuadrature(*f, from, to, 4);
    strm_out("expected:" << integral_expected);
    strm_out("actual:" << integral_actual);
    strm_out("error:" << std::abs(integral_actual - integral_expected));
  }
  f = &e;
  {
    strm_out("exponential xe^x from 0 to 5");
    const double integral_expected = 1 + 4 * std::exp(5);
    const double integral_actual = gaussian_cuadrature(*f, from, to, 4);
    strm_out("expected:" << integral_expected);
    strm_out("actual:" << integral_actual);
    strm_out("error:" << std::abs(integral_actual - integral_expected));
  }
  return 0;
}