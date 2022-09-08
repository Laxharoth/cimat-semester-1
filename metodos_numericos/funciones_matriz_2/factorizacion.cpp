#include "matrix_like/matrix_like.tcc"
#include "factorizacion.hpp"

void crout(RealMatrix &A_mtx, RealMatrix &L_mtx, RealMatrix &U_mtx){
    const size_t size = A_mtx.shape_y;
    auto calcular_factor_inferior = [&](const int &i, const int &j){
        L_mtx[i][j] = A_mtx[i][j];
        for(int k = 0; k <= j-1; ++k)
            L_mtx[i][j] -= L_mtx[i][k] * U_mtx[k][j]; 
    };
    auto calcular_factor_superior = [&](const int &i, const int &j){
        U_mtx[i][j] = A_mtx[i][j];
        for(int k = 0; k <= i-1; ++k)
            U_mtx[i][j] -= L_mtx[i][k] * U_mtx[k][j];         
        U_mtx[i][j] /= L_mtx[i][i];
    };
    for(int i=0; i<size; ++i){
        for(int j=0; j<i; ++j){
            calcular_factor_inferior(i,j);
            calcular_factor_superior(j,i);
        }
        calcular_factor_inferior(i,i);
        if(L_mtx[i][i] == 0) throw cant_factor_exception();
    }
}
void doolittle(RealMatrix &A_mtx, RealMatrix &L_mtx, RealMatrix &U_mtx){
    const size_t size = A_mtx.shape_y;
    auto calcular_factor_superior = [&](const int &i, const int &j){
        U_mtx[i][j] = A_mtx[i][j];
        for(int k = 0; k <=i-1; ++k)
            U_mtx[i][j] -= L_mtx[i][k] * U_mtx[k][j]; 
    };
    auto calcular_factor_inferior = [&](const int &i, const int &j){
        L_mtx[i][j] = A_mtx[i][j];
        for(int k = 0; k <=j-1; ++k)
            L_mtx[i][j] -= L_mtx[i][k] * U_mtx[k][j]; 
        L_mtx[i][j] /= U_mtx[j][j];
    };
    for(int j=0; j<size; ++j){
        for(int i=0; i<j; ++i){
            calcular_factor_superior(i,j);
            calcular_factor_inferior(j,i);
        }
        calcular_factor_superior(j,j);
        if(L_mtx[j][j] == 0) throw cant_factor_exception();
    }
}
void LDU_factor(RealMatrix &A_mtx, RealMatrix &L_mtx, RealMatrix &D_mtx, RealMatrix &U_mtx){
    const size_t size = A_mtx.shape_y;
    auto calcular_factor_inferior = [&](const int &i, const int &j){
        L_mtx[i][j] = A_mtx[i][j];
        for(int k = 0; k <= j-1; ++k)
            L_mtx[i][j] -= D_mtx[k][k] * L_mtx[i][k] * U_mtx[k][j]; 
        L_mtx[i][j]/= D_mtx[j][j];
    };
    auto calcular_factor_superior = [&](const int &i, const int &j){
        U_mtx[i][j] = A_mtx[i][j];
        for(int k = 0; k <= i-1; ++k)
            U_mtx[i][j] -= D_mtx[k][k] * L_mtx[i][k] * U_mtx[k][j];
        U_mtx[i][j] /= D_mtx[i][i];
    };
    auto calcular_factor_diagonal = [&](const int &i){
        D_mtx[i][i] = A_mtx[i][i];
        for(int k = 0; k <= i-1; ++k)
            D_mtx[i][i] -= D_mtx[k][k] * L_mtx[i][k] * U_mtx[k][i]; 
    };
    for(int i=0; i<size; ++i){
        for(int j=0; j<i; ++j){
            calcular_factor_inferior(i,j);
            calcular_factor_superior(j,i);
        }
        calcular_factor_diagonal(i);
        if(L_mtx[i][i] == 0) throw cant_factor_exception();
    }
}
const char* cant_factor_exception::what() const throw(){
    return "zero found in diagonal";
}

void crout_tridiagonal(RealMatrix &A_mtx, RealMatrix &L_mtx, RealMatrix &U_mtx){
    L_mtx[0][0] = A_mtx[0][0];
    if(L_mtx[0][0] == 0) throw cant_factor_exception();
    for (size_t i = 1; i < A_mtx.shape_y; i++){
        L_mtx[i][i-1] = A_mtx[i][i-1];
        U_mtx[i-1][i] = A_mtx[i-1][i] / L_mtx[i-1][i-1];
        L_mtx[i][i] = A_mtx[i][i] - L_mtx[i][i-1]*U_mtx[i-1][i];
        if(L_mtx[i-1][i-1] == 0) throw cant_factor_exception();
    }
}
void doolittle_tridiagonal(RealMatrix &A_mtx, RealMatrix &L_mtx, RealMatrix &U_mtx){
    U_mtx[0][0] = A_mtx[0][0];
    if(U_mtx[0][0] == 0) throw cant_factor_exception();
    for (size_t i = 1; i < A_mtx.shape_y; i++){
        U_mtx[i-1][i] = A_mtx[i-1][i];
        L_mtx[i][i-1] = A_mtx[i][i-1] / L_mtx[i-1][i-1];
        U_mtx[i][i] = A_mtx[i][i] - L_mtx[i][i-1]*U_mtx[i-1][i];
        if(U_mtx[i-1][i-1] == 0) throw cant_factor_exception();
    }
}
void LDU_factor_tridiagonal(RealMatrix &A_mtx, RealMatrix &L_mtx, RealMatrix &D_mtx, RealMatrix &U_mtx){
    D_mtx[0][0] = A_mtx[0][0];
    for (size_t i = 1; i < A_mtx.shape_y; i++){
        U_mtx[i-1][i] = A_mtx[i-1][i] / D_mtx[i][i];
        L_mtx[i][i-1] = A_mtx[i][i-1] / D_mtx[i][i];
        D_mtx[i][i] = A_mtx[i][i] - D_mtx[i][i]*L_mtx[i][i-1]*U_mtx[i-1][i];
        if(D_mtx[i][i] == 0) throw cant_factor_exception();
    }
}

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
