#include "factorizacion.hpp"

void crout(RealMatrix &A_mtx, RealMatrix &L_mtx, RealMatrix &U_mtx){
    const size_t size = A_mtx.shape_y;
    auto calcular_factor_inferior = [&](const int &i, const int &j){
        L_mtx(i,j) = A_mtx(i,j);
        for(int k = 0; k <= j-1; ++k)
            L_mtx(i,j) -= L_mtx(i,k) * U_mtx(k,j); 
    };
    auto calcular_factor_superior = [&](const int &i, const int &j){
        U_mtx(i,j) = A_mtx(i,j);
        for(int k = 0; k <= i-1; ++k)
            U_mtx(i,j) -= L_mtx(i,k) * U_mtx(k,j);         
        U_mtx(i,j) /= L_mtx(i,i);
    };
    for(int i=0; i<size; ++i){
        for(int j=0; j<i; ++j){
            calcular_factor_inferior(i,j);
            calcular_factor_superior(j,i);
        }
        calcular_factor_inferior(i,i);
        if(L_mtx(i,i) == 0) throw cant_factor_exception();
    }
}
void doolittle(RealMatrix &A_mtx, RealMatrix &L_mtx, RealMatrix &U_mtx){
    const size_t size = A_mtx.shape_y;
    auto calcular_factor_superior = [&](const int &i, const int &j){
        U_mtx(i,j) = A_mtx(i,j);
        for(int k = 0; k <=i-1; ++k)
            U_mtx(i,j) -= L_mtx(i,k) * U_mtx(k,j); 
    };
    auto calcular_factor_inferior = [&](const int &i, const int &j){
        L_mtx(i,j) = A_mtx(i,j);
        for(int k = 0; k <=j-1; ++k)
            L_mtx(i,j) -= L_mtx(i,k) * U_mtx(k,j); 
        L_mtx(i,j) /= U_mtx(j,j);
    };
    for(int j=0; j<size; ++j){
        for(int i=0; i<j; ++i){
            calcular_factor_superior(i,j);
            calcular_factor_inferior(j,i);
        }
        calcular_factor_superior(j,j);
        if(L_mtx(j,j) == 0) throw cant_factor_exception();
    }
}
void LDU_factor(RealMatrix &A_mtx, RealMatrix &L_mtx, RealMatrix &D_mtx, RealMatrix &U_mtx){
    const size_t size = A_mtx.shape_y;
    auto calcular_factor_inferior = [&](const int &i, const int &j){
        L_mtx(i,j) = A_mtx(i,j);
        for(int k = 0; k <= j-1; ++k)
            L_mtx(i,j) -= D_mtx(k,k) * L_mtx(i,k) * U_mtx(k,j); 
        L_mtx(i,j)/= D_mtx(j,j);
    };
    auto calcular_factor_superior = [&](const int &i, const int &j){
        U_mtx(i,j) = A_mtx(i,j);
        for(int k = 0; k <= i-1; ++k)
            U_mtx(i,j) -= D_mtx(k,k) * L_mtx(i,k) * U_mtx(k,j);
        U_mtx(i,j) /= D_mtx(i,i);
    };
    auto calcular_factor_diagonal = [&](const int &i){
        D_mtx(i,i) = A_mtx(i,i);
        for(int k = 0; k <= i-1; ++k)
            D_mtx(i,i) -= D_mtx(k,k) * L_mtx(i,k) * U_mtx(k,i); 
    };
    for(int i=0; i<size; ++i){
        for(int j=0; j<i; ++j){
            calcular_factor_inferior(i,j);
            calcular_factor_superior(j,i);
        }
        calcular_factor_diagonal(i);
        if(L_mtx(i,i) == 0) throw cant_factor_exception();
    }
}
const char* cant_factor_exception::what() const throw(){
    return "zero found in diagonal";
}
void crout_as_band(RealMatrix &A_mtx, RealMatrix &L_mtx, RealMatrix &U_mtx,
        int height, int width){
    const size_t size = A_mtx.shape_y;
    auto calcular_factor_inferior = [&](const int &i, const int &j, unsigned int begin_row, unsigned int begin_col){
        L_mtx(i,j) = A_mtx(i,j);
        for(int k = std::min(begin_row,begin_col); k <= j-1; ++k)
            L_mtx(i,j) -= L_mtx(i,k) * U_mtx(k,j); 
    };
    auto calcular_factor_superior = [&](const int &i, const int &j, unsigned int begin_row, unsigned int begin_col){
        U_mtx(i,j) = A_mtx(i,j);
        for(int k = std::min(begin_row,begin_col); k <= i-1; ++k)
            U_mtx(i,j) -= L_mtx(i,k) * U_mtx(k,j);         
        U_mtx(i,j) /= L_mtx(i,i);
    };
    size_t begin_row,begin_col;
    for(int i=0; i<size; ++i){
        begin_row = std::max(i - height , 0);
        begin_col = std::max(i - width , 0);
        for(int j=std::min(begin_col,begin_row); j<i; ++j){
            if(j >= begin_row)
                calcular_factor_inferior(i,j,begin_row,begin_col);
            if(j >= begin_col)
                calcular_factor_superior(j,i,begin_row,begin_col);
        }
        calcular_factor_inferior(i,i,begin_row,begin_col);
    }
}
void doolittle_as_band(RealMatrix &A_mtx, RealMatrix &L_mtx, RealMatrix &U_mtx,
        int height, int width){
    const size_t size = A_mtx.shape_y;
    auto calcular_factor_superior = [&](const int &i, const int &j,int begin_row, int begin_col){
        U_mtx(i,j) = A_mtx(i,j);
        for(int k = std::min(begin_row,begin_col); k <=i-1; ++k)
            U_mtx(i,j) -= L_mtx(i,k) * U_mtx(k,j); 
    };
    auto calcular_factor_inferior = [&](const int &i, const int &j,int begin_row, int begin_col){
        L_mtx(i,j) = A_mtx(i,j);
        for(int k = std::min(begin_row,begin_col); k <=j-1; ++k)
            L_mtx(i,j) -= L_mtx(i,k) * U_mtx(k,j); 
        L_mtx(i,j) /= U_mtx(j,j);
    };
    size_t begin_row,begin_col;
    for(int j=0; j<size; ++j){
        begin_col = std::max(j - width , 0);
        begin_row = std::max(j - height , 0);
        for(int i=std::min(begin_row,begin_col); i<j; ++i){
            if(i >= begin_col)
                calcular_factor_superior(i,j,begin_row,begin_col);
            if(i >= begin_row)
                calcular_factor_inferior(j,i,begin_row,begin_col);
        }
        calcular_factor_superior(j,j,begin_row,begin_col);
    }
}
void crout_tridiagonal(RealMatrix &A_mtx, RealMatrix &L_mtx, RealMatrix &U_mtx){
    L_mtx(0,0) = A_mtx(0,0);
    if(L_mtx(0,0) == 0) throw cant_factor_exception();
    for (size_t i = 1; i < A_mtx.shape_y; i++){
        L_mtx[i][i-1] = A_mtx[i][i-1];
        U_mtx[i-1][i] = A_mtx[i-1][i] / L_mtx(i-1,i-1);
        L_mtx(i,i) = A_mtx(i,i) - L_mtx[i][i-1]*U_mtx[i-1][i];
        if(L_mtx(i-1,i-1) == 0) throw cant_factor_exception();
    }
}
void doolittle_tridiagonal(RealMatrix &A_mtx, RealMatrix &L_mtx, RealMatrix &U_mtx){
    U_mtx(0,0) = A_mtx(0,0);
    if(U_mtx(0,0) == 0) throw cant_factor_exception();
    for (size_t i = 1; i < A_mtx.shape_y; i++){
        U_mtx[i-1][i] = A_mtx[i-1][i];
        L_mtx[i][i-1] = A_mtx[i][i-1] / L_mtx(i-1,i-1);
        U_mtx(i,i) = A_mtx(i,i) - L_mtx[i][i-1]*U_mtx[i-1][i];
        if(U_mtx(i-1,i-1) == 0) throw cant_factor_exception();
    }
}
void LDU_factor_tridiagonal(RealMatrix &A_mtx, RealMatrix &L_mtx, RealMatrix &D_mtx, RealMatrix &U_mtx){
    D_mtx(0,0) = A_mtx(0,0);
    for (size_t i = 1; i < A_mtx.shape_y; i++){
        U_mtx[i-1][i] = A_mtx[i-1][i] / D_mtx(i,i);
        L_mtx[i][i-1] = A_mtx[i][i-1] / D_mtx(i,i);
        D_mtx(i,i) = A_mtx(i,i) - D_mtx(i,i)*L_mtx[i][i-1]*U_mtx[i-1][i];
        if(D_mtx(i,i) == 0) throw cant_factor_exception();
    }
}

void factor_cholesky(const mymtx::RealMatrix &matix, mymtx::RealMatrix &triangular){
    auto reduce_triangular = [&](const size_t i, const size_t j){
        double &current = triangular(i,j) = matix(i,j);
        auto l_ik =triangular[i].begin();
        auto l_jk =triangular[j].begin();
        for (size_t k = 0; j != 0 && k < j; ++k, ++l_ik, ++l_jk){
            current -= (*l_ik) * (*l_jk);
        }
        current /= triangular(j,j);
        triangular(j,i) = current;
    };
    auto reduce_diagonal = [&](const size_t i){
        double &current = triangular(i,i) = matix(i,i);
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
void factor_cholesky_as_band(const mymtx::RealMatrix &matix, mymtx::RealMatrix &triangular, 
        int height){
    auto reduce_triangular = [&](const size_t i, const size_t j,int begin_col){
        double &current = triangular(i,j) = matix(i,j);
        auto l_ik =triangular[i].begin()+begin_col;
        auto l_jk =triangular[j].begin()+begin_col;
        for (size_t k = begin_col; j != 0 && k < j; ++k, ++l_ik, ++l_jk){
            current -= (*l_ik) * (*l_jk);
        }
        current /= triangular(j,j);
        triangular(j,i) = current;
    };
    auto reduce_diagonal = [&](const size_t i,int begin_col){
        double &current = triangular(i,i) = matix(i,i);
        auto l_ik =triangular[i].begin()+begin_col;
        for (size_t k = begin_col; i != 0 && k < i; ++k, ++l_ik){
            current -= (*l_ik) * (*l_ik);
        }
        current = sqrt(current);
    };
    size_t begin_col;
    for(size_t i=0; i<matix.shape_y; ++i){
        begin_col = std::max((int)i - height , 0);
        for(size_t j=begin_col;j<i; ++j){
            reduce_triangular(i,j,begin_col);
        }
        reduce_diagonal(i,begin_col);
    }
}
void factor_cholesky_tridiag(mymtx::RealMatrix &matix, mymtx::RealMatrix &triangular){
    triangular(0,0) = std::sqrt(matix(0,0));
    for(int i=1; i<matix.shape_y; i++){
        triangular[i][i-1] = matix[i][i-1]/triangular(i-1,i-1);
        triangular[i-1][i] = triangular[i][i-1];
        double substract = triangular[i][i-1];
        double &red =triangular(i,i) = matix(i,i)-substract*substract;
        red = std::sqrt(red);
    }
}
