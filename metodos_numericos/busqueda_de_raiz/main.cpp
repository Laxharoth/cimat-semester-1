#include "biseccion.cpp"
#include "newton.cpp"
#include "search_min_max.cpp"
#include "fparser/fparser.hh"
#include "function_wrapper/parser_wrapper.hpp"

#include <string>
#include <iostream>

void run_test( FunctionWrapper *func,const char*, double lim_inf, double lim_sup, double x_init );
void run_test_Newton(FunctionWrapper *func, double x_init,int iteraciones, double &error, double &res);
int main(int argc, char** argv){
    const char *f1 = "x^3-21*x^2+120*x-100";
    const char *f1_2 = "0.99*x^3-21*x^2+120*x-100";
    const char *f1_3 = "1.01*x^3-21*x^2+120*x-100";
    const char *f2 = "2 - log(x) / x";
    const char *f3 = "log( x^2 +1) - 2.71828^(0.4*x)*cos(3.1416*x)";
    FunctionParser parser;
    FunctionParserAdapter_to_FWrapper fp( &parser );

    parser.Parse(f1, "x");
    run_test(&fp,f1,0.8,1.2,1.2);
    parser.Parse(f1_2, "x");
    run_test(&fp,f1_2,0.8,1.2,1.2);
    parser.Parse(f1_3, "x");
    run_test(&fp,f1_3,0.8,1.2,1.2);
    parser.Parse(f2, "x");
    std::cout << "Function:" << f2 << std::endl;
    double error, res;
    run_test_Newton(&fp,1,1000,error,res);
    parser.Parse(f3, "x");
    run_test(&fp,f3,3.6,3.8,3.8);
    return 0;
}
void run_test_Newton(FunctionWrapper *func, double x_init,int iteraciones, double &error, double &res){
    std::cout << "Newton:" << std::endl;
    std::cout << "x inicial:" <<x_init << std::endl;
    std::cout << "valore iniciales:" <<func->eval(x_init)<< std::endl;
    res = newton(func, x_init,10000, &iteraciones,&error);
    std::cout << "raiz encontrada:" << res << std::endl;
    std::cout << "iteraciones: " << iteraciones << std::endl;
    std::cout << "error relativo: " << error << std::endl;
}
void run_test( FunctionWrapper *func,const char* fn_str, double lim_inf, double lim_sup, double x_init ){
    double res{};
    double error{};
    int iteraciones{};
    std::cout << "Function:" << fn_str << std::endl;
    std::cout << "biseccion:" << std::endl;
    std::cout << "intervalo inicial:" <<"["<<lim_inf<<","<<lim_sup<<"]" << std::endl;
    std::cout << "valores iniciales:" <<"("<<func->eval(lim_inf)<<","<<func->eval(lim_sup)<<")" << std::endl;
    res = biseccion(func, lim_inf, lim_sup,10000,&iteraciones,&error);
    std::cout << "raiz encontrada:" << res << std::endl;
    std::cout << "iteraciones" << iteraciones << std::endl;
    std::cout << "error relativo: " << error << std::endl;
    run_test_Newton(func,x_init,iteraciones,error,res);
}