#ifndef NEWTON_CPP
#define NEWTON_CPP
#include "function_wrapper/function_wrapper.hpp"
#include <cmath>

#ifndef TOLERANCIA
#define TOLERANCIA 10e-18
#endif

double newton(FunctionWrapper *funcion, double x_inicial, int max_iter,
              int *iteraciones, double *error) {
  auto calc_err = [](const double &anterior, const double &actual) {
    return std::abs((actual - anterior)) / std::abs(actual);
  };
  int inner_iter{0};
  int *real_iter;
  double y_actual{};
  Derivative derivada(funcion);
  real_iter = &inner_iter;
  if (iteraciones != nullptr) {
    (*iteraciones) = 0;
    real_iter = iteraciones;
  }
  while ((*real_iter) < max_iter) {
    ++(*real_iter);
    y_actual = funcion->eval(x_inicial);
    if (std::abs(y_actual) <= TOLERANCIA)
      return x_inicial;
    if (error != nullptr) {
      (*error) = calc_err(x_inicial, x_inicial - funcion->eval(x_inicial) /
                                                     derivada.eval(x_inicial));
    }
    x_inicial -= funcion->eval(x_inicial) / derivada.eval(x_inicial);
  }
  return x_inicial;
}

#endif /* NEWTON_CPP */
