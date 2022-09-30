#include "funcion_matriz.hpp"
#define PI 3.141592653589793
#define JACOBI_UMBRAL 10e-4
#define ZERO_UMBRAL 10e-5

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

void solucion_diagonal(const mymtx::RealMatrix &matriz, mymtx::RealVector &incognitas, mymtx::RealVector &result){
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
void inversa_diagonal(mymtx::RealMatrix &matriz_diagonal,mymtx::RealMatrix &inversa){
    const size_t size = matriz_diagonal.shape_x;
    for(int i=0; i<size; ++i){
        inversa[i][i] = 1 / matriz_diagonal[i][i];
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
        if(std::abs(divide_privote)<ZERO_UMBRAL) throw std::runtime_error("zero division");
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
void normalize(RealMatrix &mtx){
    double norm1{0},norm2{0};
    for(size_t i = 0; i < mtx.shape_y; ++i){
        norm1 += mtx[i][0]*mtx[i][0];
    }
    norm1 = std::sqrt(norm1);
    for(size_t k = 1; k < mtx.shape_x; ++k){
        for(size_t i = 0; i < mtx.shape_y; ++i){
            norm2 += mtx[i][k]*mtx[i][k];
            mtx[i][k-1] /= norm1;
        }
        norm1 = std::sqrt(norm2);
        norm2=0;
    }
    for(size_t i = 0; i < mtx.shape_y; ++i){
        norm1 += mtx[i][mtx.shape_x-1]*mtx[i][mtx.shape_x-1];
    }
}
void randomize(RealVector& vec){
    randgen &rng = randgen::get_randgen();
    for(auto i = vec.begin(); i != vec.end(); ++i){
        *i= rng.generate();
    }
}
void randomize(RealMatrix& mtx){
    randgen &rng = randgen::get_randgen();
    for(size_t i = 0; i < mtx.shape_y; ++i){
        for( auto j = mtx.begin(i); j < mtx.end(i); ++j){
            *j =rng.generate();
        }
    }
}
void power_iteration(const RealMatrix &A, RealVector &V1, 
        const double tolerance, double &value, size_t n_values, 
        RealMatrix *_vec_holder, RealVector *_val_holder,const size_t max_iter){
    RealVector V0(V1.size);
    auto &vec_holder = *_vec_holder;
    auto &val_holder = *_val_holder;
    double error = 1E20;
    double old_val{0};
    size_t found{0},iter{0};
    for(size_t k=0; k<n_values; ++k){
        iter = 0;
        randomize(V0);
        normalize(V0);
        while(error > tolerance && iter < max_iter){
            iter++;
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
        if(_vec_holder == nullptr || _val_holder == nullptr || 
            k >= n_values || iter >= max_iter){ return; }
        vec_holder[found]=V1;
        val_holder[found]=value;
        found++;
        error = 1;
    }
}
void inverse_power_iteration(const RealMatrix &A, RealVector &V1, const double tolerance, double &value, size_t n_values, RealMatrix*_vec_holder, RealVector *_val_holder, const size_t max_iter){
    RealVector V0(V1.size);
    auto &vec_holder = *_vec_holder;
    auto &val_holder = *_val_holder;
    double error = 1E20;
    double old_val{0};
    size_t found{0};
    unsigned iter = 0;
    RealMatrix working_m=A;
    crout(working_m,working_m,working_m);
    for(size_t k=0; k<n_values; ++k){
        iter = 0;
        randomize(V0);
        normalize(V0);
        while(error > tolerance && iter < max_iter){
            iter ++;
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
    mymtx::RealMatrix B = A;
    const mymtx::RealMatrix I = mymtx::RealMatrix::identity(n);
    mymtx::RealMatrix R(n,n);
    mymtx::RealVector new_row_i(n);
    mymtx::RealVector new_row_j(n);
    auto IndexOfMax = [](const mymtx::RealMatrix &inA, size_t &incol, size_t &inrow){
        double max = inrow = 0;
        incol = 1;
        for(size_t k=0; k<inA.shape_y; ++k)
        for(auto i = inA.begin(k); i<inA.end(k); ++i){
            if( i.get_col() == i.get_row() )continue;
            if(std::abs(*i) > std::abs(max)){
                max = *i;
                incol = i.get_col();
                inrow = i.get_row();
            }
        }
    };
    auto rotate =[&](mymtx::RealMatrix &U_mtx,const size_t row_lmb,const size_t col_lmb){
        double tan_2,tan_, cos_, sin_, theta;
        const size_t n = B.shape_y;
        if (row_lmb == col_lmb) return;
        if(B[row_lmb][row_lmb]!=B[col_lmb][col_lmb]){
            tan_2 = (2*B[row_lmb][col_lmb])/(B[row_lmb][row_lmb]-B[col_lmb][col_lmb]);
            tan_ = tan_2*tan_2/(1+ std::sqrt( 1 + tan_2*tan_2 ) );
            cos_ = 1/std::sqrt(tan_*tan_ + 1);
            sin_ = cos_*tan_;
        }else{
            cos_ = std::cos( PI / 4 );
            sin_ = std::sin( PI / 4 );
        }
        R = I;
        R[row_lmb][row_lmb] = R[col_lmb][col_lmb] = cos_; 
        R[row_lmb][col_lmb] = sin_; R[col_lmb][row_lmb] = -sin_;
        U_mtx *= R;  
        B = mymtx::MatrixTraspose(U_mtx)*A.prod_as_band(U_mtx,1,1);
    };
    U = mymtx::RealMatrix::identity(n);
    while(iter++<max_iter){
        IndexOfMax(B,row,col);
        if (std::abs(B[row][col]) < JACOBI_UMBRAL){ break;}
        rotate(U, row, col);
        // printf("%d\n",iter);
    }
    for(size_t i=0; i<n; ++i) e[i]=B[i][i];
}
void subspace_pow(const mymtx::RealMatrix &A,mymtx::RealMatrix &I_t,mymtx::RealVector &eig){
    size_t m = eig.size;
    mymtx::MatrixTraspose I(I_t);
    mymtx::RealMatrix B(m,m);
    mymtx::RealMatrix ro(m,m);
    mymtx::RealVector d(m);
    double val1;
    randomize(I_t);
    for (size_t i = 0; i < I_t.shape_y; i++){
        I_t[i]/=I_t[i].distance();
    }
    size_t iter,row,col;
    double error = 1;
    auto eig_old = eig;
    auto IndexOfMax = [](const mymtx::RealMatrix &inA){
        double max = 0;
        for(size_t k=0; k<inA.shape_y; ++k){
        for(auto i = inA.begin(k); i<inA.end(k); ++i){
            if( i.get_col() == i.get_row() )continue;
            if(std::abs(*i) > std::abs(max)){
                max = *i;
            }
        }}
        return max;
    };
    while (1){
        iter = 0;
        for(size_t k=0; k<m; ++k){
        iter = 0;
        mymtx::RealVector V0 = I_t[k];
        while(iter < 40){
            iter++;
            for (size_t i = 0; i < k; i++){
                V0 -= I_t[i] * (I_t[i]*V0);
            }
            V0 = A * V0;
            normalize(V0);
            }
            I_t[k] = V0;
        }
        B = I_t * A * I;
        eig_old = eig;
        jacobi_eigen(B,eig,ro,100000);
        if((eig-eig_old).distance()<ZERO_UMBRAL)break;
        I_t = mymtx::MatrixTraspose(ro) * I_t;
        for (size_t i = 0; i < I_t.shape_y; i++){
            I_t[i]/=I_t[i].distance();
        }
    }
    for (size_t i = 0; i < I_t.shape_y; i++){
        I_t[i]/=I_t[i].distance();
    }
}
void subspace_ipow(const mymtx::RealMatrix &A,mymtx::RealMatrix &I_t,mymtx::RealVector &eig){
    size_t m = eig.size;
    mymtx::MatrixTraspose I(I_t);
    mymtx::RealMatrix B(m,m);
    mymtx::RealMatrix Chol(A.shape_y,A.shape_x);
    mymtx::RealMatrix ro(m,m);
    mymtx::RealVector d(m);
    double val1;
    randomize(I_t);
    for (size_t i = 0; i < I_t.shape_y; i++){
        I_t[i]/=I_t[i].distance();
    }
    size_t iter,row,col;
    double error = 1;
    auto eig_old = eig;
    auto IndexOfMax = [](const mymtx::RealMatrix &inA){
        double max = 0;
        for(size_t k=0; k<inA.shape_y; ++k){
        for(auto i = inA.begin(k); i<inA.end(k); ++i){
            if( i.get_col() == i.get_row() )continue;
            if(std::abs(*i) > std::abs(max)){
                max = *i;
            }
        }}
        return max;
    };
    factor_cholesky(A,Chol);
    while (1){
        iter = 0;
        for(size_t k=0; k<m; ++k){
        iter = 0;
        mymtx::RealVector V0 = I_t[k];
        mymtx::RealVector V1(V0.size);
        while(iter < 40){
            iter++;
            for (size_t i = 0; i < k; i++){
                V0 -= I_t[i] * (I_t[i]*V0);
            }
            solve_cholesky(Chol,V1,V0);
            V0 = V1;
            normalize(V0);
            }
            I_t[k] = V0;
        }
        B = I_t * A * I;
        eig_old = eig;
        jacobi_eigen(B,eig,ro,100000);
        if((eig-eig_old).distance()<ZERO_UMBRAL)break;
        I_t = mymtx::MatrixTraspose(ro) * I_t;
        for (size_t i = 0; i < I_t.shape_y; i++){
            I_t[i]/=I_t[i].distance();
        }
    }
    for (size_t i = 0; i < I_t.shape_y; i++){
        I_t[i]/=I_t[i].distance();
    }
}


void rayleigh_method(const mymtx::RealMatrix &A,mymtx::RealVector &V1, double &val){
    const auto identity = mymtx::RealMatrix::identity(V1.size);
    if(((A - identity*val)*V1).distance()< ZERO_UMBRAL)return;
    normalize(V1);
    RealVector V0=V1;
    val = V0*(A * V0);
    mymtx::RealMatrix B=(A - identity*val);
    mymtx::RealMatrix Q(B.shape_y,B.shape_x),R(B.shape_y,B.shape_x);
    double old_val = val+1;
    while((B*V1).distance()>ZERO_UMBRAL && std::abs(old_val-val)>ZERO_UMBRAL){
        try{
            // crout(B,B,B);
            // solucion_crout(B,V1,V0);
            gauss(B,V1,V0);
            // qr_decomposition(B,Q,R);
            // solve_qr(Q,R,V1,V0);
        }catch(...){return;}
        normalize(V1);
        V0 = V1;
        old_val = val;
        val = V0*(A*V0);
        B = (A - identity*val);
    }
}

void qr_decomposition(const mymtx::RealMatrix& A, mymtx::RealMatrix&Q, mymtx::RealMatrix&R){
    mymtx::RealVector col0 = A.column(0).as_vector();
    R[0][0]=normalize(col0);
    Q.column(0) = col0;
    RealVector a_ux(A.shape_y);
    for (size_t j = 1; j < A.shape_x; j++){
        //calc R_ij
        for (size_t i = 0; i < j; i++){
            R[i][j] = Q.column(i) * A.column(j);
        }
        //calc a*
        a_ux = R[0][j]*Q.column(0);
        for (size_t i = 1; i < j; i++)
            a_ux += R[i][j]*Q.column(i);
        a_ux = A.column(j) - a_ux;
        //calc R_jj
        R[j][j] = normalize(a_ux);
        //assing q_j
        Q.column(j)=a_ux;
    }
}
void solve_qr(const mymtx::RealMatrix &Q, mymtx::RealMatrix &R, mymtx::RealVector &var, const mymtx::RealVector &res){
    mymtx::RealVector Qres = mymtx::MatrixTraspose(Q) * res;
    solucion_triangular_sup(R,var,Qres);
}
void conjugate_gradient(const mymtx::RealMatrix &A, mymtx::RealVector &x, mymtx::RealVector &b){
    auto r=b;
    auto p=r;
    unsigned k{0};
    mymtx::RealVector w(x.size);
    double alpha,betha;
    while (k < r.size * 2){
        w = A*p;
        alpha = (p*r)/(p*w);
        x = x + alpha*p;
        r = r - alpha*w;
        if(std::sqrt(r*r)< ZERO_UMBRAL) break;
        betha = (p*r)/(p*p);
        p = r + betha*p;
    }
    
}
void conjugate_gradient_jacobi(const mymtx::RealMatrix &A, mymtx::RealVector &x, mymtx::RealVector &b){
    double alpha,betha;
    auto r = b;
    mymtx::RealVector z(x.size);
    mymtx::RealVector w(x.size);
    solucion_diagonal(A,z,r);
    auto p = z;
    unsigned i{0};
    double rzold;
    while( i++< x.size * 3 ){
        rzold = r*z;
        w = A*p;
        alpha = (rzold)/(p*w);
        x   = x + alpha*p;
        r   = r - alpha*w;
        if(std::sqrt(r*r)< ZERO_UMBRAL) break;
        solucion_diagonal(A,z,r);
        betha =(r*z)/(rzold);
        p = r + (betha)*p;
    }
}
