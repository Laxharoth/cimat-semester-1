#ifndef NEWTON_CPP
#define NEWTON_CPP
#include "function_wrapper/function_wrapper.hpp"
#include <cmath>

#ifndef TOLERANCIA
#define TOLERANCIA 10e-18
#endif

double newton( FunctionWrapper *funcion, double x_inicial ){
    Derivative derivada(funcion);
    double y_actual{};
    while(1){
        y_actual = funcion->eval(x_inicial);
        if(std::abs(y_actual) <= TOLERANCIA) return y_actual;
        x_inicial -= funcion->eval(x_inicial) / derivada.eval(x_inicial);
    }
}

#endif /* NEWTON_CPP */
