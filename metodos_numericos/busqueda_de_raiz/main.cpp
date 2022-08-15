#include "biseccion.cpp"
#include "newton.cpp"
#include "search_min_max.cpp"
#include "fparser/fparser.hh"
#include "function_wrapper/parser_wrapper.hpp"

#include <string>
#include <iostream>

int main(int argc, char** argv){
    const char* funcion_str = argv[1];
    const double x_inferior = std::stod(argv[2]);
    const double x_superior = std::stod(argv[3]);

    std::cout << funcion_str << std::endl;

    FunctionParser parser;
    FunctionParserAdapter_to_FWrapper fp( &parser );

    parser.Parse(funcion_str, "x");
    try{
        double raiz_biseccion = biseccion( &fp, x_inferior, x_superior );
        std::cout << "Raiz encontrada en (" << raiz_biseccion <<"," << fp.eval(raiz_biseccion) <<  ") por biseccion" << std::endl;
        double raiz_newton = newton( &fp, x_superior);
        std::cout << "Raiz encontrada en (" << raiz_newton <<"," << fp.eval(raiz_newton) <<  ") por Newton" << std::endl;
        const char* min_max_funcion = "sin(x)";
        parser.Parse(min_max_funcion,"x");
        auto min_max_values = search_min_max(&fp, x_inferior, 10);
        std::cout << "minimos y maximos encontrados en "<< min_max_funcion <<" desde " << x_inferior << "hasta" << x_superior << std::endl;
        for(auto &&value : *min_max_values) {
            std::cout << "raiz : " << value.x_value<< ", ";
            switch(value.type) {
                case MIN: std::cout << "type : "<< "min" << std::endl; break;
                case MAX: std::cout << "type : "<< "max" << std::endl; break;
                case SILLA: std::cout << "type : "<< "silla" << std::endl; break;
            }
        }
    }
    catch(const char* e){
        if(e == ERROR_BAD_LIMITS){
            std::cout << "No se puede encontrar una raiz_biseccion con los limites proporcionados";
            return 1;
        }
        throw e;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
    return 0;
}
