#ifndef FUNCTION_WRAPPER_HPP
#define FUNCTION_WRAPPER_HPP

#ifndef DELTA_X
#define DELTA_X 0.00001
#endif

class FunctionWrapper{
    public:
    virtual double eval(const double &x) = 0;
};

class Derivative : public FunctionWrapper{
    FunctionWrapper *original_function;
    public:
    Derivative(FunctionWrapper *original);
    double eval(const double &x);
};

#endif /* FUNCTION_WRAPPER_HPP */
