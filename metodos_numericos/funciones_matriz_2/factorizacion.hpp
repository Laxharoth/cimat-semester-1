#ifndef FACTORIZACION_HPP
#define FACTORIZACION_HPP

#include <exception>

void metodo_de_crout(matrix_like<double> &matriz, matrix_like<double> &matriz_inferior, matrix_like<double> &matriz_superior, const int &size);
void metodo_de_doolittle(matrix_like<double> &matriz, matrix_like<double> &matriz_inferior, matrix_like<double> &matriz_superior, const int &size);
void factorizacion_LDU(matrix_like<double> &matriz, matrix_like<double> &matriz_inferior, matrix_like<double> &matriz_diagonal, matrix_like<double> &matriz_superior, const int &size);

class cant_factor_exception : public std::exception{
    using std::exception::exception;
    public:
    virtual const char* 
        what() 
        const throw();
};

#endif /* FACTORIZACION_HPP */
