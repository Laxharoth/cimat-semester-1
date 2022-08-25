#ifndef FACTORIZACION_HPP
#define FACTORIZACION_HPP

#include <exception>

void metodo_de_crout(matrix_like<double> &matriz, matrix_like<double> &matriz_inferior, matrix_like<double> &matriz_superior, const int &size);
void metodo_de_doolittle(matrix_like<double> &matriz, matrix_like<double> &matriz_inferior, matrix_like<double> &matriz_superior, const int &size);
void factorizacion_LDU(matrix_like<double> &matriz, matrix_like<double> &matriz_inferior, matrix_like<double> &matriz_diagonal, matrix_like<double> &matriz_superior, const int &size);

class cant_factor_exception : public std::exception{
    using std::exception::exception;
    public:
    virtual const char* 
        what() 
        const throw();
};

class LDU_wrapper : public matrix_like<double>{
    matrix_like *data;
protected:
    LDU_wrapper(matrix_like *data, char type);
    class array_wrapper :public array_like<double>{
        public:
        matrix_like *data;
        static double default_val;
        size_t row;
        explicit array_wrapper(matrix_like<double> *data);
        virtual double &operator[](const size_t &col) = 0;
    };
    class diagonal_strategy : public array_wrapper{
        using array_wrapper::array_wrapper;
        double &operator[](const size_t &col);
    };
    class crout_strategy : public array_wrapper{
        using array_wrapper::array_wrapper;
        double &operator[](const size_t &col);
    };
    class doolittle_strategy : public array_wrapper{
        using array_wrapper::array_wrapper;
        double &operator[](const size_t &col);
    };
    array_wrapper *arr_wrapper;
    public:
    const static char DIAGONAL = 0b000;
    const static char CROUT = 0b001;
    const static char DOOLITTLE = 0b010;
    virtual array_like<double> &operator[](const size_t &row);
    static LDU_wrapper from(matrix_like *data, char type);
};


#endif /* FACTORIZACION_HPP */
