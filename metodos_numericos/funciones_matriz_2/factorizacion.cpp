#include "matrix_like/matrix_like.tcc"
#include "factorizacion.hpp"

void metodo_de_crout(matrix_like<double> &matriz, matrix_like<double> &matriz_inferior, matrix_like<double> &matriz_superior, const int &size){
    factorizacion_LDU(matriz, matriz_inferior, matriz_inferior, matriz_superior, size);
}
void metodo_de_doolittle(matrix_like<double> &matriz, matrix_like<double> &matriz_inferior, matrix_like<double> &matriz_superior, const int &size){
    auto calcular_factor_superior = [&](const int &i, const int &j){
        matriz_superior[i][j] = matriz[i][j];
        for(int k = 0; k <=i-1; ++k)
            matriz_superior[i][j] -= matriz_inferior[i][k] * matriz_superior[k][j]; 
    };
    auto calcular_factor_inferior = [&](const int &i, const int &j){
        matriz_inferior[i][j] = matriz[i][j];
        for(int k = 0; k <=j-1; ++k)
            matriz_inferior[i][j] -= matriz_inferior[i][k] * matriz_superior[k][j]; 
        matriz_inferior[i][j] /= matriz_superior[j][j];
    };
    for(int j=0; j<size; ++j){
        for(int i=0; i<j; ++i){
            calcular_factor_superior(i,j);
            calcular_factor_inferior(j,i);
        }
        calcular_factor_superior(j,j);
        if(matriz_inferior[j][j] == 0) throw cant_factor_exception();
    }
}
void factorizacion_LDU(matrix_like<double> &matriz, matrix_like<double> &matriz_inferior, matrix_like<double> &matriz_diagonal, matrix_like<double> &matriz_superior, const int &size){
    auto calcular_factor_inferior = [&](const int &i, const int &j){
        matriz_inferior[i][j] = matriz[i][j];
        for(int k = 0; k <= j-1; ++k)
            matriz_inferior[i][j] -= matriz_diagonal[k][k] * matriz_inferior[i][k] * matriz_superior[k][j]; 
        matriz_inferior[i][j]/= matriz_diagonal[j][j];
    };
    auto calcular_factor_superior = [&](const int &i, const int &j){
        matriz_superior[i][j] = matriz[i][j];
        for(int k = 0; k <= i-1; ++k)
            matriz_superior[i][j] -= matriz_diagonal[k][k] * matriz_inferior[i][k] * matriz_superior[k][j];
        matriz_superior[i][j] /= matriz_diagonal[i][i];
    };
    auto calcular_factor_diagonal = [&](const int &i, const int &j){
        matriz_diagonal[i][j] = matriz[i][j];
        for(int k = 0; k <= j-1; ++k)
            matriz_diagonal[i][j] -= matriz_diagonal[k][k] * matriz_inferior[i][k] * matriz_superior[k][j]; 
    };
    for(int i=0; i<size; ++i){
        for(int j=0; j<i; ++j){
            calcular_factor_inferior(i,j);
            calcular_factor_superior(j,i);
        }
        calcular_factor_diagonal(i,i);
        if(matriz_inferior[i][i] == 0) throw cant_factor_exception();
    }
}
const char* cant_factor_exception::what() const throw(){
    return "zero found in diagonal";
}

LDU_wrapper::LDU_wrapper(matrix_like *data,char type):matrix_like(data->get_shape_y(), data->get_shape_x()),data(data){
    switch (type){
    case LDU_wrapper::DIAGONAL:arr_wrapper = new LDU_wrapper::diagonal_strategy(data);
        break;
    case LDU_wrapper::CROUT:arr_wrapper = new LDU_wrapper::crout_strategy(data);
        break;
    case LDU_wrapper::DOOLITTLE:arr_wrapper = new LDU_wrapper::doolittle_strategy(data);
        break;
    default:
    //TODO: make this exception
        throw cant_factor_exception();
        break;
    }
}
LDU_wrapper LDU_wrapper::from(matrix_like *data,char type){
    return LDU_wrapper(data, type);
}
array_like<double> &LDU_wrapper::operator[](const size_t &row){
    (*(this->arr_wrapper)).row = row;
    return *(this->arr_wrapper);
}
LDU_wrapper::array_wrapper::array_wrapper(matrix_like *data):array_like<double>(data->get_shape_x()),data(data),row(0){}
double LDU_wrapper::array_wrapper::default_val = 0;
double &LDU_wrapper::diagonal_strategy::operator[](const size_t &col){
    if(row == col){
        return (*data)[row][col];
    }
    default_val = 0;
    return default_val;
}
size_t LDU_wrapper::diagonal_strategy::get_rbegin_n() const { return row; }
size_t LDU_wrapper::diagonal_strategy::get_rend_n() const { return row+1; }
double &LDU_wrapper::crout_strategy::operator[](const size_t &col){
    if(row == col){
        default_val = 1;
        return default_val;
    }
    if(row > col){
        default_val = 0;
        return default_val;
    }
    return (*data)[row][col];
}
size_t LDU_wrapper::crout_strategy::get_rbegin_n() const { return row; }
size_t LDU_wrapper::crout_strategy::get_rend_n() const { return get_size(); }
double &LDU_wrapper::doolittle_strategy::operator[](const size_t &col){
    if(row == col){
        default_val = 1;
        return default_val;
    }
    if(row < col){
        default_val = 0;
        return default_val;
    }
    return (*data)[row][col];
}
size_t LDU_wrapper::doolittle_strategy::get_rbegin_n() const { return 0; }
size_t LDU_wrapper::doolittle_strategy::get_rend_n() const { return row+1; }
