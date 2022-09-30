#include "functions.hpp"

#include <cmath>
#include <cstdio>

void factor_cholesky(mymtx::RealMatrix &matix, mymtx::RealMatrix &triangular){
    auto reduce_triangular = [&](const size_t i, const size_t j){
        double &current = triangular[i][j] = matix[i][j];
        auto l_ik =triangular[i].begin();
        auto l_jk =triangular[j].begin();
        for (size_t k = 0; j != 0 && k < j; ++k, ++l_ik, ++l_jk){
            current -= (*l_ik) * (*l_jk);
        }
        current /= triangular[j][j];
        triangular[j][i] = current;
    };
    auto reduce_diagonal = [&](const size_t i){
        double &current = triangular[i][i] = matix[i][i];
        auto l_ik =triangular[i].begin();
        for (size_t k = 0; i != 0 && k < i; ++k, ++l_ik){
            current -= (*l_ik) * (*l_ik);
        }
        current = sqrt(current);
    };
    for(size_t i=0; i<matix.shape_y; ++i){
        for(size_t j=0;j<i; ++j){
            reduce_triangular(i,j);
        }
        reduce_diagonal(i);
    }
}
void factor_cholesky_tridiag(mymtx::RealMatrix &matix, mymtx::RealMatrix &triangular){
    triangular[0][0] = std::sqrt(matix[0][0]);
    for(int i=1; i<matix.shape_y; i++){
        triangular[i][i-1] = matix[i][i-1]/triangular[i-1][i-1];
        triangular[i-1][i] = triangular[i][i-1];
        double substract = triangular[i][i-1];
        double &red =triangular[i][i] = matix[i][i]-substract*substract;
        red = std::sqrt(red);
    }
}
#define ITER 1000
#define TOLER 1E-10
void solve_jacobi(mymtx::RealMatrix &matix, mymtx::RealVector &variables, mymtx::RealVector &solutions, double *error){
    double sum, numerator, denominator,toler;
    const size_t n = solutions.size;
    size_t i,j;
    mymtx::RealVector vv(n);
    mymtx::RealVector &vn=variables;
    for( size_t i = 0; i < n; ++i ){
        vv[i] = solutions[i] / matix[i][i];
    }
for(  unsigned int iter = 0; iter < ITER; ++iter){
    numerator = denominator = 0;
    for( i = 0; i < n; ++i){
        sum =  0;
        for(j = 0; j < n; ++j){
            if(j != i){ sum += matix[i][j] * vv[j]; }
        }
        vn[i] = (solutions[i] - sum) / matix[i][i];
        numerator += std::abs(vn[i]-vv[i]);
        denominator += std::abs(vn[i]);
    }
    toler = numerator / denominator;
    if(toler < TOLER){ if(error!=nullptr)*error = toler; return;}
    for( size_t i = 0; i < n; ++i ){
        vv[i] = vn[i];
    }
}
}
void solve_gauss_seidel(mymtx::RealMatrix &matix, mymtx::RealVector &variables, mymtx::RealVector &solutions, double *error){
    double sum, numerator, denominator,toler;
    const size_t n = solutions.size;
    size_t i,j;
    mymtx::RealVector vv(n);
    mymtx::RealVector &vn=variables;
    for( size_t i = 0; i < n; ++i ){
        vv[i] = vn[i] = solutions[i] / matix[i][i];
    }
for(  unsigned int iter = 0; iter < ITER; ++iter){
    numerator = denominator = 0;
    for( i = 0; i < n; ++i){
        sum =  0;
        for(j = 0; j < n; ++j){
            if(j != i){ sum += matix[i][j] * vn[j]; }
        }
        vn[i] = (solutions[i] - sum) / matix[i][i];
        numerator += std::abs(vn[i]-vv[i]);
        denominator += std::abs(vn[i]);
    }
    toler = numerator / denominator;
    if(toler < TOLER){ if(error!=nullptr)*error = toler; return;}
    for( size_t i = 0; i < n; ++i ){
        vv[i] = vn[i];
    }
}
}

void solve_cholesky(mymtx::RealMatrix &cholesky_factored,mymtx::RealVector &variables, mymtx::RealVector &solutions){
    mymtx::RealVector tmp(solutions.size);
    solucion_triangular_inf(cholesky_factored,tmp,solutions);
    solucion_triangular_sup(cholesky_factored,variables,tmp);
}
