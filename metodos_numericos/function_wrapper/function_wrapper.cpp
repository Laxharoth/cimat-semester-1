#ifndef FUNCTION_WRAPPER_CPP
#define FUNCTION_WRAPPER_CPP
#include "function_wrapper.hpp"

Derivative::Derivative(FunctionWrapper* original):original_function(original){}

double Derivative::eval(const double &x){
    return ( original_function->eval( x + DELTA_X ) - original_function->eval(x)  ) / DELTA_X;
}

#endif /* FUNCTION_WRAPPER_CPP */
