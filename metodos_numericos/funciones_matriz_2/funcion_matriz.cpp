#include "funcion_matriz.hpp"
#include "factorizacion.hpp"

void solucion_diagonal(mymtx::RealMatrix &matriz, mymtx::RealVector &incognitas, mymtx::RealVector &result){
    const size_t size = matriz.shape_x;
    for(int i=0; i<size; ++i){
        incognitas[i] = result[i] / matriz[i][i];
    }
}
double determinante_diagonal(mymtx::RealMatrix &matriz_diagonal){
    const size_t size = matriz_diagonal.shape_x;
    double result{1};
    for(int i=0; i<size; ++i) result*=matriz_diagonal[i][i];
    return result;
}
void inversa_diagonal(mymtx::RealMatrix &matriz_diagonal,mymtx::RealVector &inversa){
    const size_t size = matriz_diagonal.shape_x;
    for(int i=0; i<size; ++i){
        inversa[i] = 1 / matriz_diagonal[i][i];
    }
}
double determinante_triangular(mymtx::RealMatrix &matriz_triangular){
    const size_t size = matriz_triangular.shape_x;
    double result{1};
    for(int i=0; i<size; ++i) result*=matriz_triangular[i][i];
    return result;
}
void solucion_triangular_inf( mymtx::RealMatrix &matriz, mymtx::RealVector &incognitas, mymtx::RealVector &result){
    solucion_triangular_inf(matriz,incognitas,result,true);
}
void solucion_triangular_inf( mymtx::RealMatrix &matriz, mymtx::RealVector &incognitas, mymtx::RealVector &result,bool compute_diagonale){
    const size_t size = matriz.shape_x;
    for(int i=0; i<size; ++i){
        incognitas[i] = result[i];
        auto iter = matriz[i].begin();
        for(int j=0; j<=i - 1; ++j){
            incognitas[i] -= matriz[i][j] * incognitas[j];
        }
        if(compute_diagonale)
            incognitas[i] /= matriz[i][i];
    }
}
void solucion_triangular_sup( mymtx::RealMatrix &matriz, mymtx::RealVector &incognitas, mymtx::RealVector &result){
    solucion_triangular_sup(matriz,incognitas,result,true);
}
void solucion_triangular_sup( mymtx::RealMatrix &matriz, mymtx::RealVector &incognitas, mymtx::RealVector &result, bool compute_diagonale){
    const size_t size = matriz.shape_x;
    for(int i = size - 1; i>=0; --i){
        incognitas[i] = result[i];
        for(int j=i+1; j< size; ++j){
            incognitas[i] -= matriz[i][j] * incognitas[j];
        }
        if(compute_diagonale)
            incognitas[i] /= matriz[i][i];
    }
}
void gauss( mymtx::RealMatrix &matriz, mymtx::RealVector &variables, mymtx::RealVector &resultados){
    const size_t size = matriz.shape_x;
    for(int i=0; i<size; ++i){
        const double divide_privote = matriz[i][i];
        for(auto j=matriz[i].begin()+i; j<matriz[i].end(); ++j){
            *j /= divide_privote;
        }
        resultados[i] /= divide_privote;
        for(int j=i+1; j < size; ++j){
            auto iter_pivote = matriz[i].begin()+(i+1);
            auto iter = matriz[j].begin()+i;
            auto end  = matriz[j].end();
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
void solucion_LDU( mymtx::RealMatrix &matriz, mymtx::RealVector &incognitas, mymtx::RealVector &result){
    mymtx::RealVector aux1(incognitas.size);
    mymtx::RealVector aux2(incognitas.size);
    solucion_triangular_inf( matriz, aux1, result,false);
    solucion_diagonal( matriz, aux2, aux1);
    solucion_triangular_sup( matriz, incognitas,aux2, false);
}
void solucion_crout( mymtx::RealMatrix &matriz, mymtx::RealVector &incognitas, mymtx::RealVector &result){
    mymtx::RealVector aux1(incognitas.size);
    solucion_triangular_inf( matriz, aux1, result);
    solucion_triangular_sup( matriz, incognitas, aux1, false);
}
void solucion_doolittle( mymtx::RealMatrix &matriz, mymtx::RealVector &incognitas, mymtx::RealVector &result){
    mymtx::RealVector aux1(incognitas.size);
    solucion_triangular_inf( matriz, aux1, result, false );
    solucion_triangular_sup( matriz, incognitas, aux1);
}
double normalize(RealVector &vec){
    double sum{0};
    for(auto i = vec.begin(); i != vec.end(); ++i){
        sum += *i**i;
    }
    sum = std::sqrt(sum);
    for(auto i = vec.begin(); i != vec.end(); ++i){
        *i/=sum;
    }
    return sum;
}
void power_iteration(const RealMatrix &A, RealVector &V0, RealVector &V1, const double tolerance, double &value, size_t n_values, RealVector *vec_holder, double *val_holder){
    double error = 1E20;
    double old_val{0};
    size_t found{0};
    RealVector original = V0;
    double norm;
    double *component = new double[n_values];
    RealMatrix working_m = A;
    for(size_t k=0; k<n_values; ++k){
        while(error > tolerance){
            V1 = working_m*V0;
            value = V0*V1;
            error = std::abs(old_val-value);
            old_val = value;
            norm = normalize(V1);
            V0 = V1;
        }
        if(vec_holder == nullptr || val_holder == nullptr || k >= n_values){ return; }
        V0 = original;
        component[found] = V1*V0;
        vec_holder[found] = V1 * norm;
        val_holder[found++] = value;
        working_m -= V1.cross_product(V1) * (value * norm);
        error = 1;
    }
    delete[] component;
}

void inverse_power_iteration(const RealMatrix &A, RealVector &V0, RealVector &V1, const double tolerance, double &value, size_t n_values, RealVector *vec_holder, double *val_holder){
    double error = 1E20;
    double old_val{0};
    size_t found{0};
    RealVector original = V0;
    double *component = new double[n_values];
    RealMatrix working_m=A;
    metodo_de_crout( working_m,working_m,working_m );
    for(size_t k=0; k<n_values; ++k){
        while(error > tolerance){
            solucion_crout(working_m,V1,V0);
            value = 1 / (V0*V1);
            error = std::abs(old_val-value);
            old_val = value;
            normalize(V1);
            V0 = V1;
        }
        if(vec_holder == nullptr || val_holder == nullptr || k >= n_values){ return; }
        V0 = original;
        component[found] = V1*V0;
        vec_holder[found] = V1;
        val_holder[found++] = value;
        working_m -= V1.cross_product(V1) * value;
        error = 1;
    }
    delete[] component;
}

