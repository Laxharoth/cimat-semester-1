#include <stdio.h>

#define scan_double_prompt(prompt, variable) printf(prompt);scanf("%lf",&variable);
#define scan_char_prompt(prompt, variable) printf(prompt);variable = getchar();

double suma_division(double a, double b, double c, double d, char _operator){
    switch (_operator){
    case '+':
        return a / b + c / d;
    case '-':
        return a / b - c / d;
    case '*':
        return a / b * c / d;
    case '/':
        return a / b / c / d;
    default:
        return -1;
    }
}

int main(int argc, char const *argv[]){
    double a,b,c,d;
    char oper;
    scan_double_prompt("Introduce a: ", a);
    scan_double_prompt("Introduce b: ", b);
    scan_double_prompt("Introduce c: ", c);
    scan_double_prompt("Introduce d: ", d);
    //consume line
    getchar();
    scan_char_prompt("Introduce operador: ", oper);
    printf("%f/%f%c%f/%f=%f\n",a,b,oper,c,d,suma_division(a,b,c,d,'+') );
}
