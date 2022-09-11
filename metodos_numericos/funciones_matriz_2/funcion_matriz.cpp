#include "funcion_matriz.hpp"
#define PI 3.141592653589793
#define ZERO_UMBRAL 10e-15

class randgen{
    std::random_device *rd;
    std::mt19937 *gen;
    std::uniform_real_distribution<double> *dis;
    randgen(){
        rd = new std::random_device();
        gen= new std::mt19937((*rd)());
        dis = new std::uniform_real_distribution<double>();
    }
    ~randgen(){
        delete rd;
        delete gen;
        delete dis;
    }
    public:
    double generate(){ return (*dis)(*gen); }
    static randgen &get_randgen(){
        static randgen rng;
        return rng;
    }
};

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
void solucion_triangular_inf_as_band( mymtx::RealMatrix &matriz, mymtx::RealVector &incognitas, mymtx::RealVector &result,bool compute_diagonale, size_t heigh){
    const size_t size = matriz.shape_x;
    for(int i=0; i<size; ++i){
        incognitas[i] = result[i];
        auto iter = matriz[i].begin();
        size_t start = 0;
        if( i > heigh ) start = i - heigh;
        for(int j=start; j<=i - 1; ++j){
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
void solucion_triangular_sup_as_band( mymtx::RealMatrix &matriz, mymtx::RealVector &incognitas, mymtx::RealVector &result, bool compute_diagonale, size_t width){
    const size_t size = matriz.shape_x;
    for(int i = size - 1; i>=0; --i){
        incognitas[i] = result[i];
        size_t end = size;
        if( i < width) end = width - i;
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
void solve_crout_as_band( mymtx::RealMatrix &matriz, mymtx::RealVector &incognitas, mymtx::RealVector &result, const size_t heigh, const size_t width){
    mymtx::RealVector aux1(incognitas.size);
    solucion_triangular_inf_as_band( matriz, aux1, result,true,heigh);
    solucion_triangular_sup_as_band( matriz, incognitas, aux1, false,width);
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
void randomize(RealVector& vec){
    randgen &rng = randgen::get_randgen();
    for(auto i = vec.begin(); i != vec.end(); ++i){
        *i= rng.generate();
    }
}
void power_iteration(const RealMatrix &A, RealVector &V0, RealVector &V1, const double tolerance, double &value, size_t n_values, RealMatrix *_vec_holder, RealVector *_val_holder){
    auto &vec_holder = *_vec_holder;
    auto &val_holder = *_val_holder;
    double error = 1E20;
    double old_val{0};
    size_t found{0};
    for(size_t k=0; k<n_values; ++k){
        randomize(V0);
        normalize(V0);
        while(error > tolerance){
            for (size_t i = 0; i < found; i++){
                V0 -= vec_holder[i] * (V0*vec_holder[i]);
            }
            V1 = A.prod_as_band(V0,1,1);
            value = V0*V1;
            error = std::abs(old_val-value);
            old_val = value;
            normalize(V1);
            V0 = V1;
        }
        if(_vec_holder == nullptr || _val_holder == nullptr || k >= n_values){ return; }
        vec_holder[found]=V1;
        val_holder[found]=value;
        found++;
        error = 1;
    }
}

void inverse_power_iteration(const RealMatrix &A, RealVector &V0, RealVector &V1, const double tolerance, double &value, size_t n_values, RealMatrix*_vec_holder, RealVector *_val_holder){
    auto &vec_holder = *_vec_holder;
    auto &val_holder = *_val_holder;
    double error = 1E20;
    double old_val{0};
    size_t found{0};
    const unsigned max_iter = 2000;
    unsigned iter = 0;
    RealMatrix working_m=A;
    crout(working_m,working_m,working_m);
    for(size_t k=0; k<n_values; ++k){
        iter = 0;
        while(error > tolerance && iter < max_iter){
            iter ++;
            randomize(V0);
            normalize(V0);
            for (size_t i = 0; i < found; i++){
                V0 -= vec_holder[i] * (V0*vec_holder[i]);
            }
            solucion_crout(working_m,V1,V0);
            value = 1 / (V0*V1);
            error = std::abs(old_val-value);
            old_val = value;
            normalize(V1);
            V0 = V1;
        }
        if(_vec_holder == nullptr || _val_holder == nullptr || 
            k >= n_values || iter >= max_iter){ return; }
        vec_holder[found]=V1;
        val_holder[found]=value;
        found++;
        error = 1;
    }
}
void solve_cholesky(mymtx::RealMatrix &cholesky_factored,mymtx::RealVector &variables, mymtx::RealVector &solutions){
    mymtx::RealVector tmp(solutions.size);
    solucion_triangular_inf(cholesky_factored,tmp,solutions);
    solucion_triangular_sup(cholesky_factored,variables,tmp);
}

void jacobi_eigen(mymtx::RealMatrix &A, mymtx::RealVector &e, mymtx::RealMatrix  &U, const unsigned max_iter){
    size_t col, row;
    const size_t n = A.shape_y;
    unsigned iter = 0;
    auto IndexOfMax = [](const mymtx::RealMatrix &A, size_t &col, size_t &row){
        double max = row = col= 0;
        for(size_t k=0; k<A.shape_y; ++k)
        for(auto i = A.begin(k); i<A.end(k); ++i)
            if(*i > max){
                max = *i;
                col = i.get_col();
                row = i.get_row();
            }
    };
    auto rotate =[](mymtx::RealMatrix &A_mtx, mymtx::RealMatrix &U_mtx,const size_t row_lmb,const size_t col_lmb){
        double tan_2,tan_, cos_, sin_, theta;
        int n = A_mtx.shape_y;
        if (row_lmb == col_lmb) return;
        if(A_mtx[row_lmb][row_lmb]!=A_mtx[col_lmb][col_lmb]){
            tan_2 = (2*A_mtx[row_lmb][col_lmb])/(A_mtx[row_lmb][row_lmb]-A_mtx[col_lmb][col_lmb]);
            tan_ = tan_2*tan_2/(1+ std::sqrt( 1 + tan_2*tan_2 ) );
            cos_ = 1/sqrt(tan_*tan_ + 1);
            sin_ = cos_*tan_;
        }else{
            cos_ = std::cos( PI / 4 );
            sin_ = std::sin( PI / 4 );
        }
        auto R = mymtx::RealMatrix::identity(n);
        R[row_lmb][row_lmb] = R[col_lmb][col_lmb] = cos_; R[row_lmb][col_lmb] = sin_; R[col_lmb][row_lmb] = -sin_;
        A_mtx = (mymtx::RealMatrix::traspose(R)*A_mtx);
        A_mtx = (A_mtx * R);
        U_mtx = (R*U_mtx);
    };
    mymtx::RealMatrix M(n,n);
    U = mymtx::RealMatrix::identity(n);
    while(iter++ < max_iter){
        M=A;
        mymtx::abs(M);
        for(size_t i=0; i<n; ++i) M[i][i]=0;
        IndexOfMax(M,row,col);
        if(row==col){
            for(size_t i=0; i<n; ++i) e[i]=A[i][i];
            return;
        }
        double Amax = A[row][col];
        rotate(A, U, row, col);
        if (Amax < ZERO_UMBRAL) break;
    }
    for(size_t i=0; i<n; ++i) e[i]=A[i][i];
}
