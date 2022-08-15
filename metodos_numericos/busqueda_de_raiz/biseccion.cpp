#ifndef BISECCION_CPP
#define BISECCION_CPP
#include "function_wrapper/function_wrapper.hpp"
#include <algorithm>
#include <cmath>

#ifndef TOLERANCIA
#define TOLERANCIA 10e-18
#endif

const char* ERROR_BAD_LIMITS = "5672ba31-9e2a-4796-bf7e-dc61acd048bb";

double biseccion( FunctionWrapper *funcion, double x_inferior, double x_superior ){
    if( x_inferior == x_superior ) throw ERROR_BAD_LIMITS;
    if( x_inferior > x_superior ){ std::swap( x_inferior, x_superior); }
    double y_inferior = funcion->eval( x_inferior );
    double y_superior = funcion->eval( x_superior );	
    if( std::abs(y_inferior) <= TOLERANCIA) return x_inferior;
    if( std::abs(y_superior) <= TOLERANCIA) return x_superior;
    double x_medio{}, y_medio{};
    while(1){
        if( y_inferior * y_superior > 0 ) throw ERROR_BAD_LIMITS;
        x_medio = ( x_inferior + x_superior ) / 2;
        y_medio = funcion->eval(x_medio);
        if(std::abs(y_medio) <= TOLERANCIA) return x_medio;
        if(y_medio * y_inferior > 0) x_inferior = x_medio;
        else x_superior = x_medio;
    }
}

#endif /* BISECCION_CPP */
