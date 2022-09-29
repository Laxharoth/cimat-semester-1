#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP
#include "matrix_like/real_matrix.hpp"
#include "function_wrapper/function_wrapper.hpp"
#include "funcion_matriz.hpp"
#include <vector>
#include <cmath>

class PolyFunction : public FunctionWrapper{
    const mymtx::RealVector A;
    public:
    double eval(const double &x);
    double eval(const double &x) const;
    PolyFunction(const mymtx::RealVector &v);
};
class LagramFunction : public FunctionWrapper{
    mymtx::RealVector X;
    mymtx::RealVector Y;
    public: double eval(const double &x);
    double eval(const double &x) const;
    LagramFunction(const mymtx::RealVector &X, const mymtx::RealVector &Y);
};
class NewtonPolyFunction : public FunctionWrapper{
    mymtx::RealVector X;
    mymtx::RealVector Y;
    mymtx::RealVector as;
    public: double eval(const double &x);
    double eval(const double &x) const;
    NewtonPolyFunction(const mymtx::RealVector &X, const mymtx::RealVector &Y);
};
class MultiFunctionWrapper : public FunctionWrapper{
    std::vector<FunctionWrapper *> fns;
    mymtx::RealVector coef;
    public: double eval(const double &x);
    double eval(const double &x) const;
    MultiFunctionWrapper(std::vector<FunctionWrapper *> fns, mymtx::RealVector coef);
};

PolyFunction interpolate_line(const mymtx::RealVector &X, const mymtx::RealVector &Y);
PolyFunction interpolate_poly(const mymtx::RealVector &X, const mymtx::RealVector &Y, unsigned int grade);
MultiFunctionWrapper interpolate_funcs(const mymtx::RealVector &X, const mymtx::RealVector &Y, std::vector<FunctionWrapper *>fns);
PolyFunction interpolate_poly_2(const mymtx::RealVector &X, const mymtx::RealVector &Y);
LagramFunction interpolate_lagram(const mymtx::RealVector &X, const mymtx::RealVector &Y);
NewtonPolyFunction interpolate_newton(const mymtx::RealVector &X, const mymtx::RealVector &Y);

namespace mymtx{
    RealVector map(const RealVector &v, FunctionWrapper *fn);
    RealVector map(const RealVector &v, const FunctionWrapper &fn);
}

#endif /* INTERPOLATION_HPP */
