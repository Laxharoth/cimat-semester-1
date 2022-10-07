#include "biseccion.cpp"
#include "function_wrapper/function_wrapper.hpp"
#include <cmath>
#include <memory>
#include <vector>

#ifndef TOLERANCIA
#define TOLERANCIA 10e-18
#endif

const double SEARCH_MIN_MAX_DELTA = 10E-2;

typedef enum {
  MIN = 1,
  MAX = 1 << 1,
  SILLA = 1 << 2,
} MIN_MAX_TYPE;

struct min_max_value {
  double x_value;
  MIN_MAX_TYPE type;
};

std::unique_ptr<std::vector<min_max_value>>
search_min_max(FunctionWrapper *fw, const double &x_inferior,
               const double &x_superior) {
  std::unique_ptr<std::vector<min_max_value>> ptr_search_min_max =
      std::make_unique<std::vector<min_max_value>>();
  std::vector<min_max_value> &puntos = *ptr_search_min_max;

  Derivative derivada(fw);
  Derivative derivada_2(&derivada);

  for (double x_izquierda = x_inferior; x_izquierda < x_superior;
       x_izquierda += SEARCH_MIN_MAX_DELTA) {
    double y_izquierda = derivada.eval(x_izquierda);
    double y_derecha = derivada.eval(x_izquierda + SEARCH_MIN_MAX_DELTA);
    if (y_izquierda * y_derecha > 0.0)
      continue;
    double raiz =
        biseccion(&derivada, x_izquierda, x_izquierda + SEARCH_MIN_MAX_DELTA,
                  1000, nullptr, nullptr);
    double aux_seleccionar_tipo = derivada_2.eval(raiz);
    if (aux_seleccionar_tipo > TOLERANCIA) {
      puntos.push_back(min_max_value{raiz, MIN});
      continue;
    }
    if (aux_seleccionar_tipo < TOLERANCIA) {
      puntos.push_back(min_max_value{raiz, MAX});
      continue;
    }
    puntos.push_back(min_max_value{raiz, SILLA});
  }
  return ptr_search_min_max;
}
