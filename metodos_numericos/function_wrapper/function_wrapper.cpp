#ifndef FUNCTION_WRAPPER_CPP
#define FUNCTION_WRAPPER_CPP
#include "function_wrapper.hpp"
#include <exception>
#include <vector>

double FunctionWrapper::operator()(const double &x) { return this->eval(x); }
double FunctionWrapper::operator()(const double &x) const {
  return this->eval(x);
}
Derivative::Derivative(FunctionWrapper *original)
    : original_function(original) {}

double Derivative::eval(const double &x) {
  return const_cast<Derivative *>(this)->eval(x);
}
double Derivative::eval(const double &x) const {
  return (original_function->eval(x + DELTA_X) - original_function->eval(x)) /
         DELTA_X;
}
class cant_be_const : std::exception {
  const char *what() const noexcept { return "Instance can't be const"; }
};
Count::Count() : current(-1), step(1) {}
Count::Count(double start) : current(start - 1), step(1) {}
Count::Count(double start, double step) : current(start - step), step(step) {}
Count::Count(double start, double end, unsigned int steps)
    : step((end - start) / steps), current(0) {
  current = start - step;
}
double Count::eval(const double &x) { return current += step; }
double Count::eval(const double &x) const { throw cant_be_const(); }

LambdaWrapper::LambdaWrapper(std::function<double(const double &)> fn)
    : fn(fn) {}
double LambdaWrapper::eval(const double &x) { return fn(x); }
double LambdaWrapper::eval(const double &x) const { return fn(x); }
std::vector<double>
MultiVarFunctionWrapper::operator()(const std::vector<double> &x) {
  return eval(x);
}
std::vector<double>
MultiVarFunctionWrapper::operator()(const std::vector<double> &x) const {
  return eval(x);
}
Gradient::Gradient(MultiVarFunctionWrapper *const original)
    : original_function(original) {}
std::vector<double> Gradient::eval(const std::vector<double> &x) {
  std::vector<double> result;
  result.reserve(x.size());
  for (unsigned int i = 0; i < x.size(); ++i) {
    auto cpy = x;
    cpy[i] += DELTA_X;
    result[i] =
        ((*original_function)(cpy)[0] - (*original_function)(x)[0]) / DELTA_X;
  }
  return result;
}
std::vector<double> Gradient::eval(const std::vector<double> &x) const {
  std::vector<double> result;
  result.reserve(x.size());
  for (unsigned int i = 0; i < x.size(); ++i) {
    auto cpy = x;
    cpy[i] += DELTA_X;
    result[i] =
        ((*original_function)(cpy)[0] - (*original_function)(x)[0]) / DELTA_X;
  }
  return result;
}

#endif /* FUNCTION_WRAPPER_CPP */
