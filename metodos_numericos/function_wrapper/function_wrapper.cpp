#ifndef FUNCTION_WRAPPER_CPP
#define FUNCTION_WRAPPER_CPP
#include "function_wrapper.hpp"

double FunctionWrapper::operator()(const double &x) { return this->eval(x); }
double FunctionWrapper::operator()(const double &x) const {
  return this->eval(x);
}
Derivative::Derivative(FunctionWrapper *original)
    : original_function(original) {
  this->dy = [this](const double &x) {
    return (this->original_function->eval(x + DELTA_X) -
            this->original_function->eval(x)) /
           DELTA_X;
  };
}
class NotDerivativeStrategy : public std::exception {
  char *err;

public:
  NotDerivativeStrategy(Derivative::DerivativeStrategy s) {
    std::stringstream ss;
    ss << s << "is not a valid Derivative Strategy";
    err = new char[100];
    std::strcpy(err, ss.str().c_str());
  };
  ~NotDerivativeStrategy() { delete[] err; }
  const char *what() const noexcept { return err; }
};
Derivative::Derivative(FunctionWrapper *original, DerivativeStrategy strategy)
    : original_function(original) {
  switch (strategy) {
  case DerivativeStrategy::FORWARD:
    this->dy = [this](const double &x) {
      return (this->original_function->eval(x + DELTA_X) -
              this->original_function->eval(x)) /
             DELTA_X;
    };
    break;
  case DerivativeStrategy::BACKWARD:
    this->dy = [this](const double &x) {
      return (this->original_function->eval(x) -
              this->original_function->eval(x - DELTA_X)) /
             DELTA_X;
    };
    break;
  case DerivativeStrategy::CENTRAL:
    this->dy = [this](const double &x) {
      return (this->original_function->eval(x + DELTA_X) -
              this->original_function->eval(x - DELTA_X)) /
             (DELTA_X * 2);
    };
    break;
  case DerivativeStrategy::FIVE_POINT_1ST:
    this->dy = [this](const double &x) {
      auto &f = *(this->original_function);
      return (-f(x + 2 * DELTA_X) + 8 * f(x + DELTA_X) - 8 * f(x - DELTA_X) +
              f(x - 2 * DELTA_X)) /
             (12 * DELTA_X);
    };
    break;
  case DerivativeStrategy::FIVE_POINT_2ND:
    this->dy = [this](const double &x) {
      auto &f = *(this->original_function);
      return (f(x - DELTA_X) - 2 * f(x) + f(x + DELTA_X)) / (DELTA_X * DELTA_X);
    };
    break;
  case DerivativeStrategy::FIVE_POINT_3RD:
    this->dy = [this](const double &x) {
      auto &f = *(this->original_function);
      return (f(x + 2 * DELTA_X) - 2 * f(x + DELTA_X) + 2 * f(x - DELTA_X) -
              f(x - 2 * DELTA_X)) /
             (2 * DELTA_X * DELTA_X * DELTA_X);
    };
    break;
  case DerivativeStrategy::FIVE_POINT_4TH:
    this->dy = this->dy = [this](const double &x) {
      auto &f = *(this->original_function);
      return (f(x + 2 * DELTA_X) - 4 * f(x + DELTA_X) + 6 * f(x) -
              4 * f(x - DELTA_X) - f(x - 2 * DELTA_X)) /
             (DELTA_X * DELTA_X);
    };
    break;
  default:
    throw NotDerivativeStrategy(strategy);
  }
}

double Derivative::eval(const double &x) {
  return const_cast<const Derivative *>(this)->eval(x);
}
double Derivative::eval(const double &x) const { return dy(x); }
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
  std::vector<double> result(x.size());
  for (unsigned int i = 0; i < x.size(); ++i) {
    auto cpy = x;
    cpy[i] += DELTA_X;
    result[i] =
        ((*original_function)(cpy)[0] - (*original_function)(x)[0]) / DELTA_X;
  }
  return result;
}
std::vector<double> Gradient::eval(const std::vector<double> &x) const {
  std::vector<double> result(x.size());
  for (unsigned int i = 0; i < x.size(); ++i) {
    auto cpy = x;
    cpy[i] += DELTA_X;
    result[i] =
        ((*original_function)(cpy)[0] - (*original_function)(x)[0]) / DELTA_X;
  }
  return result;
}
std::vector<double> operator*(const std::vector<double> &x0, const double x1) {
  auto cpy = x0;
  for (auto &&v : cpy) {
    v *= x1;
  }
  return cpy;
}
std::vector<double> operator*(const double x1, const std::vector<double> &x0) {
  return x0 * x1;
}
std::vector<double> operator/(const std::vector<double> &x0, const double x1) {
  const double recip = 1 / x1;
  return x0 * recip;
}
std::vector<double> operator/(const double x1, const std::vector<double> &x0) {
  return x0 / x1;
}
std::vector<double> operator-(const std::vector<double> &x0,
                              const std::vector<double> &x1) {
  std::vector<double> cpy = x0;
  for (unsigned int i = 0; i < cpy.size(); ++i) {
    cpy[i] -= x1[i];
  }
  return cpy;
}
#endif /* FUNCTION_WRAPPER_CPP */
