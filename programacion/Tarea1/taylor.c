#include <stdlib.h>
#include <stdio.h>

double taylor_e_x(double x, unsigned int iteraciones, double *factorial_memo, unsigned int *memo_size);
double factorial(unsigned int n, double *factorial_memo, unsigned int *factorial_memo_size);
double mpow(double x, unsigned int p);

int main(int argc, char const *argv[]){
    unsigned int factorial_memo_size = 98;
    double *factorial_memo=NULL;
    double exponent;
    unsigned int iter;
    printf("valor de x:");
    scanf("%lf",&exponent);
    printf("iteraciones:");
    scanf("%ld",&iter);
    printf("e^%f = %f\n",-(exponent*exponent), taylor_e_x(exponent,iter,factorial_memo,&factorial_memo_size));
    if(factorial_memo != NULL)
        free(factorial_memo);
    return 0;
}

double factorial(unsigned int n, double *factorial_memo, unsigned int *factorial_memo_size){
    if(factorial_memo == NULL) factorial_memo=calloc(*factorial_memo_size , sizeof(double));
    if(n > *factorial_memo_size){
        *factorial_memo_size = n+100;
        factorial_memo = realloc(factorial_memo, sizeof(double)*(*factorial_memo_size) );
    }
    if(n == 0 || n == 1) return 1;
    if( factorial_memo[n-2] == 0 ) factorial_memo[n-2] = n * factorial(n-1,factorial_memo,factorial_memo_size);
    return factorial_memo[n-2];
}
double mpow(double x, unsigned int p){
    double result=1;
    for (size_t i = 0; i < p; i++){
        result *= x;
    }
    return result;
}
double taylor_e_x(double x, unsigned int iteraciones, double *factorial_memo, unsigned int *memo_size){
    double e=0;
    for(unsigned int i=0; i<iteraciones; i++){
        e+= mpow(-(x*x),i) / factorial(i,factorial_memo,memo_size);
    }
    return e;
}
