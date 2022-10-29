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
  double eval(const double &x) const { return x + x * x; }
  double eval(const double &x) { return x + x * x; }
};

int main(int argc, const char **argv) {
  class exp e;
  polyn p;

  FunctionWrapper *f = &p;
  strm_out("poly");
  measure_time(strm_out(
      "Newton-Cotes grade 2    : " << integral_newton_cotes(*f, 0, 1, 2)));
  measure_time(strm_out("Richardson Extrapolation: "
                        << richardson_extrapolation(*f, 0, 1, 10, 1e-5, 2)));
  measure_time(strm_out(
      "Romberg Method          : " << romberg_method(*f, 0, 1, 10, 1e-5)));

  f = &e;
  strm_out("exponencial");
  measure_time(strm_out(
      "Newton-Cotes grade 2    : " << integral_newton_cotes(*f, 0, 1, 2)));
  measure_time(strm_out("Richardson Extrapolation: "
                        << richardson_extrapolation(*f, 0, 1, 10, 1e-5, 6)));
  measure_time(strm_out(
      "Romberg Method          : " << romberg_method(*f, 0, 1, 10, 1e-5)));
  return 0;
}