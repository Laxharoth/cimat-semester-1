#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP
#include "function_wrapper/function_wrapper.hpp"
#include "matrix_like/funcion_matriz.hpp"
#include "matrix_like/matrix.hpp"
#include <algorithm>
#include <cmath>
#include <vector>
#define EP 1E-8

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
class FiniteElement : public FunctionWrapper {
  std::vector<point> points;
  mymtx::vector phi;
  double increment;

public:
  FiniteElement(const std::vector<point> &points, const unsigned int nodes,
                const double lambda);
  double eval(const double &x) const;
  double eval(const double &x);
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

enum sign { POSITIVE = 1, NEGATIVE = -1 };

struct fn_interval {
  const double x0;
  const double x1;
  sign area_sign;
};

double bisection(FunctionWrapper &fn, double, double);
double area_montecarlo(FunctionWrapper &fn, const double x0, const double x1);
std::vector<double> NewtonMultivar(MultiVarFunctionWrapper &fn,
                                   const std::vector<double> &start_guess);
double volum_montecarlo(MultiVarFunctionWrapper &fn, const point p0,
                        const point p1);
#endif /* INTERPOLATION_HPP */
