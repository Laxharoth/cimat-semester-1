#ifndef FUNCTION_WRAPPER_CPP
#define FUNCTION_WRAPPER_CPP
#include "function_wrapper.hpp"

Derivative::Derivative(FunctionWrapper* original):original_function(original){}

double Derivative::eval(const double &x){
    return ( original_function->eval( x + DELTA_X ) - original_function->eval(x)  ) / DELTA_X;
}

Count::Count():current(-1),step(1){}
Count::Count(double start):current(start-1),step(1){}
Count::Count(double start,double step):current(start-step),step(step){}
Count::Count(double start, double end,unsigned int steps):step((end-start)/steps),current(start-step){}
double Count::eval(const double &x){ return current+=step; }
double Count::eval(const double &x) const { return current; }
#endif /* FUNCTION_WRAPPER_CPP */
