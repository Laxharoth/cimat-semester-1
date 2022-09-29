#ifndef FUNCTION_WRAPPER_HPP
#define FUNCTION_WRAPPER_HPP

#ifndef DELTA_X
#define DELTA_X 0.00001
#endif

class FunctionWrapper{
    public:
    virtual double eval(const double &x) = 0;
    virtual double eval(const double &x) const = 0;
};

class Derivative : public FunctionWrapper{
    FunctionWrapper *original_function;
    public:
    Derivative(FunctionWrapper *original);
    double eval(const double &x);
    double eval(const double &x) const;
};

class Count : public FunctionWrapper{
    double current,step;
    public:
    Count();
    Count(double start);
    Count(double start,double step);
    Count(double start, double end,unsigned int steps);
    double eval(const double &x);
    double eval(const double &x) const;
};

#endif /* FUNCTION_WRAPPER_HPP */
