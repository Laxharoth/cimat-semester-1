#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP
#include "function_wrapper/function_wrapper.hpp"
#include "matrix_like/funcion_matriz.hpp"
#include "matrix_like/matrix.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

class PolyFunction : public FunctionWrapper {
  const mymtx::vector A;

public:
  double eval(const double &x);
  double eval(const double &x) const;
  PolyFunction(const mymtx::vector &v);
};
class LagramFunction : public FunctionWrapper {
  mymtx::vector X;
  mymtx::vector Y;

public:
  double eval(const double &x);
  double eval(const double &x) const;
  LagramFunction(const mymtx::vector &X, const mymtx::vector &Y);
};
class NewtonPolyFunction : public FunctionWrapper {
  mymtx::vector X;
  mymtx::vector Y;
  mymtx::vector as;

public:
  double eval(const double &x);
  double eval(const double &x) const;
  NewtonPolyFunction(const mymtx::vector &X, const mymtx::vector &Y);
};
class MultiFunctionWrapper : public FunctionWrapper {
  std::vector<FunctionWrapper *> fns;
  mymtx::vector coef;

public:
  double eval(const double &x);
  double eval(const double &x) const;
  MultiFunctionWrapper(std::vector<FunctionWrapper *> fns, mymtx::vector coef);
};
class LineSpline : public FunctionWrapper {
  std::vector<point> points;

public:
  double eval(const double &x);
  double eval(const double &x) const;
  LineSpline(const std::vector<point> &points);
};
class CuadraticSpline : public FunctionWrapper {
  std::vector<point> points;
  mymtx::vector Si;

public:
  double eval(const double &x);
  double eval(const double &x) const;
  CuadraticSpline(const std::vector<point> &points);
};

class CubicSpline : public FunctionWrapper {
  std::vector<point> points;
  mymtx::vector Si;

public:
  double eval(const double &x);
  double eval(const double &x) const;
  CubicSpline(const std::vector<point> &points);
};

PolyFunction interpolate_line(const mymtx::vector &X, const mymtx::vector &Y);
PolyFunction interpolate_poly(const mymtx::vector &X, const mymtx::vector &Y,
                              unsigned int grade);
MultiFunctionWrapper interpolate_funcs(const mymtx::vector &X,
                                       const mymtx::vector &Y,
                                       std::vector<FunctionWrapper *> fns);
PolyFunction interpolate_poly_2(const mymtx::vector &X, const mymtx::vector &Y);
LagramFunction interpolate_lagram(const mymtx::vector &X,
                                  const mymtx::vector &Y);
NewtonPolyFunction interpolate_newton(const mymtx::vector &X,
                                      const mymtx::vector &Y);
// splines

namespace mymtx {
vector map(const vector &v, FunctionWrapper *fn);
vector map(const vector &v, const FunctionWrapper &fn);
} // namespace mymtx

#endif /* INTERPOLATION_HPP */
