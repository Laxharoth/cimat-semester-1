#ifndef FACTORIZACION_HPP
#define FACTORIZACION_HPP

#include <exception>

void metodo_de_crout(double **matriz, double **matriz_inferior, double **matriz_superior, const int &size);
void metodo_de_doolittle(double **matriz, double **matriz_inferior, double **matriz_superior, const int &size);
void factorizacion_LDU(double **matriz, double **matriz_inferior,double **matriz_diagonal, double **matriz_superior, const int &size);

class cant_factor_exception : public std::exception{
    using std::exception::exception;
    public:
    virtual const char* 
        what() 
        const throw();
};

#endif /* FACTORIZACION_HPP */
