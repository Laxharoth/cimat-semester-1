#include "funcion_matriz.hpp"
#include "gauss_jordan.hpp"
#include "factorizacion.hpp"
#include "matriz_banda.hpp"

int main(int argc, char **argv){
    MatrizBanda m(5, 5, 10);
    m[0][10] = 100;
    
    return 0;
}