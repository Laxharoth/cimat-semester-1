#include "interpolation.hpp"
#include "matrix_like/real_matrix.hpp"
PolyFunction interpolate_line(const mymtx::RealVector &X, const mymtx::RealVector &Y){
    return interpolate_poly(X,Y,1);
}
PolyFunction interpolate_poly(const mymtx::RealVector &X, const mymtx::RealVector &Y, unsigned int grade){
    mymtx::RealVector ab(grade+1);
    mymtx::RealVector yx(grade+1);
    mymtx::RealVector YX(Y.size);
    mymtx::RealVector Xpow(Y.size);
    mymtx::RealMatrix A(grade+1,grade+1);
    double xpow;

    for (size_t i = 0; i <= grade; i++){
        for(size_t j = 0; j < X.size; j++) YX[j] = Y[j]*pow(X[j], i);
        yx[i] = mymtx::reduce(YX, [](double acc, double cur){return acc+cur;},0);
        
    }
    for (size_t i = 0; i <= grade*2; i++){
        for (size_t j = 0; j < Xpow.size; j++){
            Xpow[j] = std::pow(X[j],i);
        }
        xpow = mymtx::reduce(Xpow, [](double acc, double cur){return acc+cur;},0);
        size_t current = i, qty = current + 1;
        size_t start = 0;
        if( current >= A.shape_x ){
            start = current - A.shape_x;
            current = A.shape_x - 1;
        }
        for (size_t j = start; j < qty && j < A.shape_y; ++j, --current){
            A(current,j) = xpow;
        }
    }
    gauss(A,ab,yx);
    // factor_cholesky(A,A);
    // solve_cholesky(A,ab,yx);
    return PolyFunction(ab);
}
MultiFunctionWrapper interpolate_funcs(const mymtx::RealVector &X, const mymtx::RealVector &Y, std::vector<FunctionWrapper *>fns){
    mymtx::RealVector cs(fns.size());
    mymtx::RealVector ys(cs.size);
    mymtx::RealMatrix A(cs.size, cs.size);
    for (size_t i = 0; i < fns.size(); i++){
        for (size_t j = 0; j < X.size; j++){
            A(i,0) += fns[i]->eval(X[j]) * fns[0]->eval(X[j]);
            ys[i]  += fns[i]->eval(X[j]);
        }
        for (size_t j = 1; j < fns.size(); j++){
            for (size_t k = 0; k < X.size; k++){
                A(i,j) += fns[i]->eval(X[k]) * fns[j]->eval(X[k]);
            }
        }    
    }
    gauss(A,cs,ys);
    return MultiFunctionWrapper(fns,cs);
}
PolyFunction interpolate_poly_2(const mymtx::RealVector &X, const mymtx::RealVector &Y){
    mymtx::RealMatrix A(X.size,X.size);
    mymtx::RealVector ys=Y;
    mymtx::RealVector as=Y;
    for (size_t i = 0; i < X.size; i++){
        for (size_t j = 0; j < X.size; j++){
            A(i,j) = std::pow(X[i],j);
        }
    }
    gauss(A,as,ys);
    return PolyFunction(as);
}
LagramFunction interpolate_lagram(const mymtx::RealVector &X, const mymtx::RealVector &Y){
    return LagramFunction(X,Y);
}
NewtonPolyFunction interpolate_newton(const mymtx::RealVector &X, const mymtx::RealVector &Y){
    return NewtonPolyFunction(X,Y);
}
double PolyFunction::eval(const double &x) {
    return const_cast<PolyFunction *>(this)->eval(x);
}
double PolyFunction::eval(const double &x) const{
    double sum = 0;
    for(auto a = A.begin(); a != A.end();++a)
        sum += *a * std::pow(x,a.get_col());
    
    return sum;
}
PolyFunction::PolyFunction(const mymtx::RealVector &v):A(mymtx::RealVector(v)){}
double MultiFunctionWrapper::eval(const double &x){
    return const_cast<MultiFunctionWrapper *>(this)->eval(x);
}
double MultiFunctionWrapper::eval(const double &x) const{
    double sum = 0;
    auto c = coef.begin();
    for(auto i = fns.begin(); i != fns.end(); ++i,++c){
        sum += (*i)->eval(x) * (*c);
    }
    return sum;
}
MultiFunctionWrapper::MultiFunctionWrapper(std::vector<FunctionWrapper *> fns, mymtx::RealVector coef):fns(fns),coef(coef){}
LagramFunction::LagramFunction(const mymtx::RealVector &X, const mymtx::RealVector &Y):X(mymtx::RealVector(X)),Y(mymtx::RealVector(Y)){}
double LagramFunction::eval(const double &x){
    return const_cast<LagramFunction *>(this)->eval(x);
}
double LagramFunction::eval(const double &x) const {
    auto Li = [&](const double x, const unsigned int i){
        double prod = 1;
        for (size_t j = 0; j < X.size; j++){
            if(i==j) continue;
            prod *= (x-X[j]) /(X[i]-X[j]);
        }
        return prod;
    };
    double sum = 0;
    for (size_t i = 0; i < X.size; i++){
        sum += Y[i]*Li(x,i);
    }
    return sum;
}
double NewtonPolyFunction::NewtonPolyFunction::eval(const double &x){
    return const_cast<NewtonPolyFunction *>(this)->eval(x);
}
double NewtonPolyFunction::NewtonPolyFunction::eval(const double &x) const{
    double acumulative_prod = 1;
    double sum = as[0];
    for (size_t i = 1; i < X.size; i++){
        acumulative_prod *= (x-X[i-1]);
        sum += acumulative_prod * as[i];
    }
    return sum;
}
double divided_difference(mymtx::RealMatrix &computed,const mymtx::RealVector &X,const mymtx::RealVector &Y, unsigned int from, unsigned int to){
    if(from == to) computed(from,to) = Y[from];
    if(!computed(from,to) && from != to){ 
        computed(from,to) = (divided_difference(computed,X,Y,from + 1, to) - divided_difference(computed,X,Y,from,to-1)) / (X[to] - X[from]);
    }
    return computed(from,to);
};
NewtonPolyFunction::NewtonPolyFunction(const mymtx::RealVector &X, const mymtx::RealVector &Y):X(mymtx::RealVector(X)),Y(mymtx::RealVector(Y)),as(mymtx::RealVector(X.size)){
    mymtx::RealMatrix computed(X.size,X.size);
    divided_difference(computed,X,Y,0,X.size-1);
    as = computed[0];
}

namespace mymtx{
    RealVector map(const RealVector &v, FunctionWrapper *fn){
        RealVector v_new(v.size);
        auto i = v.begin();
        for(auto j =v_new.begin(); j !=v_new.end();++i,++j)
            *j = fn->eval(*i);
        return v_new;
    }
    RealVector map(const RealVector &v, const FunctionWrapper &fn){
        RealVector v_new(v.size);
        auto i = v.begin();
        for(auto j =v_new.begin(); j !=v_new.end();++i,++j)
            *j = fn.eval(*i);
        return v_new;
    }
}
