#ifndef FUNCION_MATRIZ_HPP
#define FUNCION_MATRIZ_HPP

void solucion_diagonal(double *matriz, double *vector, double *result, int size);
double determinante_diagonal(double *matriz_diagonal, const int &size);
void inversa_diagonal(double *matriz_diagonal,double *inversa, const int &size);
double determinante_triangular(double **matriz_triangular, const int &size);
void solucion_triangular_inf( double **matriz, double *vector, double *result, const int &size);
void solucion_triangular_sup( double **matriz, double *incognitas, double *result, const int &row);

#endif /* FUNCION_MATRIZ_HPP */
