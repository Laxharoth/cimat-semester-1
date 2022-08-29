#include "funcion_matriz.hpp"

void solucion_diagonal(matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result){
    const size_t size = matriz.get_shape_x();
    for(int i=0; i<size; ++i){
        incognitas[i] = result[i] / matriz[i][i];
    }
}
double determinante_diagonal(matrix_like<double> &matriz_diagonal){
    const size_t size = matriz_diagonal.get_shape_x();
    double result{1};
    for(int i=0; i<size; ++i) result*=matriz_diagonal[i][i];
    return result;
}
void inversa_diagonal(matrix_like<double> &matriz_diagonal,array_like<double> &inversa){
    const size_t size = matriz_diagonal.get_shape_x();
    for(int i=0; i<size; ++i){
        inversa[i] = 1 / matriz_diagonal[i][i];
    }
}
double determinante_triangular(matrix_like<double> &matriz_triangular){
    const size_t size = matriz_triangular.get_shape_x();
    double result{1};
    for(int i=0; i<size; ++i) result*=matriz_triangular[i][i];
    return result;
}
void solucion_triangular_inf( matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result){
    const size_t size = matriz.get_shape_x();
    for(int i=0; i<size; ++i){
        incognitas[i] = result[i];
        auto iter = matriz[i].begin();
        for(int j=0; j<=i - 1; ++j){
            incognitas[i] -= matriz[i][j] * incognitas[j];
        }
        incognitas[i] /= matriz[i][i];
    }
}
void solucion_triangular_sup( matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result){
    const size_t size = matriz.get_shape_x();
    for(int i = size - 1; i>=0; --i){
        incognitas[i] = result[i];
        for(int j=i+1; j< size; ++j){
            incognitas[i] -= matriz[i][j] * incognitas[j];
        }
        incognitas[i] /= matriz[i][i];
    }
}
void gauss( matrix_like<double> &matriz, array_like<double> &variables, array_like<double> &resultados){
    const size_t size = matriz.get_shape_x();
    for(int i=0; i<size; ++i){
        const double divide_privote = matriz[i][i];
        for(auto j=matriz[i].begin()+i; j<matriz[i].rend(); ++j){
            *j /= divide_privote;
        }
        resultados[i] /= divide_privote;
        for(int j=i+1; j < size; ++j){
            auto iter_pivote = matriz[i].begin()+(i+1);
            auto iter = matriz[j].begin()+i;
            auto end  = matriz[j].rend();
            const double coeficiente_eliminar = *iter;
            *iter = 0;
            for(++iter; iter < end; ++iter,++iter_pivote){
                *iter -= coeficiente_eliminar * (*iter_pivote);
            }
            resultados[j] -= coeficiente_eliminar * resultados[i] / matriz[i][i];
        }
    }
    solucion_triangular_sup(matriz, variables, resultados);
}
void solucion_LDU( matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result){
    auto matriz_doolittle = LDU_wrapper::from(&matriz, LDU_wrapper::DOOLITTLE);
    auto matriz_crout = LDU_wrapper::from(&matriz, LDU_wrapper::CROUT);
    mymtx::vector<double> aux1(incognitas.get_size());
    mymtx::vector<double> aux2(incognitas.get_size());
    solucion_triangular_inf( matriz_doolittle, aux1, result);
    solucion_diagonal( matriz, aux2, aux1);
    solucion_triangular_sup( matriz_crout, incognitas, aux2);
}
void solucion_crout( matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result){
    auto matriz_crout = LDU_wrapper::from(&matriz, LDU_wrapper::CROUT);
    mymtx::vector<double> aux1(incognitas.get_size());
    solucion_triangular_inf( matriz, aux1, result);
    solucion_triangular_sup( matriz_crout, incognitas, aux1);
}
void solucion_doolittle( matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result){
    auto matriz_doolittle = LDU_wrapper::from(&matriz, LDU_wrapper::DOOLITTLE);
    mymtx::vector<double> aux1(incognitas.get_size());
    solucion_triangular_inf( matriz_doolittle, aux1, result );
    solucion_triangular_sup( matriz, incognitas, aux1);
}
