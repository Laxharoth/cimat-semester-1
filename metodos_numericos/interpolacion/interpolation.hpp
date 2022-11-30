#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP
#include "function_wrapper/function_wrapper.hpp"
#include "matrix_like/funcion_matriz.hpp"
#include "matrix_like/matrix.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <functional>
#include <tuple>
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
class FiniteElement2 : public MultiVarFunctionWrapper {
  std::vector<point3d> points;
  mymtx::vector phi;
  double increment;
  point start;
  double increment_x, increment_y;
  unsigned int nodes_x;
  unsigned int nodes_y;
  double N1(double ro, double nu) const;
  double N2(double ro, double nu) const;
  double N3(double ro, double nu) const;
  double N4(double ro, double nu) const;

public:
  FiniteElement2(const std::vector<point3d> &points, const point start,
                 const point end, const unsigned int nodes_x,
                 const unsigned int nodes_y, const double lambda_x,
                 const double lambda_y);
  std::vector<double> eval(const std::vector<double> &x) override;
  std::vector<double> eval(const std::vector<double> &x) const override;
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

double bisection(FunctionWrapper &fn, double, double);
double area_montecarlo(FunctionWrapper &fn, const double x0, const double x1);
std::vector<double> NewtonMultivar(MultiVarFunctionWrapper &fn,
                                   const std::vector<double> &start_guess);
double volum_montecarlo(MultiVarFunctionWrapper &fn, const point p0,
                        const point p1);
double integral_newton_cotes(FunctionWrapper &fn, const double from,
                             const double to, const unsigned grade);
typedef std::function<double(FunctionWrapper &, const double, double,
                             const unsigned)>
    aproximation_function;
double richardson_extrapolation(FunctionWrapper &f, aproximation_function aprox,
                                double from, double to, const int maxRows,
                                const double tolerance,
                                const unsigned int grade);
double romberg_method(FunctionWrapper &fn, const double a, const double b,
                      const int max_iter, const double toler);
double gaussian_cuadrature(FunctionWrapper &fn, const double from,
                           const double to, const unsigned grade);
#endif /* INTERPOLATION_HPP */
