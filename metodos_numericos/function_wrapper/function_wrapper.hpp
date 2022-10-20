#ifndef FUNCTION_WRAPPER_HPP
#define FUNCTION_WRAPPER_HPP

#include <functional>
#include <vector>

#ifndef DELTA_X
#define DELTA_X 0.00001
#endif

class FunctionWrapper {
public:
  virtual double eval(const double &x) = 0;
  virtual double eval(const double &x) const = 0;
  double operator()(const double &x);
  double operator()(const double &x) const;
};
class LambdaWrapper : public FunctionWrapper {
  std::function<double(const double &)> fn;

public:
  LambdaWrapper(std::function<double(const double &)> fn);
  double eval(const double &x);
  double eval(const double &x) const;
};
class Derivative : public FunctionWrapper {
  FunctionWrapper *original_function;

public:
  Derivative(FunctionWrapper *original);
  double eval(const double &x);
  double eval(const double &x) const;
};
class Count : public FunctionWrapper {
  double current, step;

public:
  Count();
  Count(double start);
  Count(double start, double step);
  Count(double start, double end, unsigned int steps);
  double eval(const double &x);
  double eval(const double &x) const;
};
struct point {
  double x;
  double y;
};
class MultiVarFunctionWrapper {
  virtual std::vector<double> eval(const std::vector<double> &x) = 0;
  virtual std::vector<double> eval(const std::vector<double> &x) const = 0;

public:
  std::vector<double> operator()(const std::vector<double> &x);
  std::vector<double> operator()(const std::vector<double> &x) const;
};
class Gradient : public MultiVarFunctionWrapper {
  MultiVarFunctionWrapper *const original_function;

public:
  Gradient(MultiVarFunctionWrapper *original);
  virtual std::vector<double> eval(const std::vector<double> &x);
  virtual std::vector<double> eval(const std::vector<double> &x) const;
};
double norm(std::vector<double> &);
std::vector<double> operator-(const std::vector<double> &x0,
                              const std::vector<double> &x1);
std::vector<double> operator*(const std::vector<double> &x0, const double x1);
std::vector<double> operator*(const double x1, const std::vector<double> &x0);
std::vector<double> operator/(const std::vector<double> &x0, const double x1);
std::vector<double> operator/(const double x1, const std::vector<double> &x0);
#endif /* FUNCTION_WRAPPER_HPP */
