/***********************************
 *                                 *
 *                                 *
 *                                 *
 *                                 *
 *                                 *
***********************************/
#include "funcion_matriz.hpp"

#include <iostream>
#include <array>

#define matrix_to_ptr(matrix,ptr,size) for(int i = 0; i < size; ++i) ptr[i] = matrix[i];

template <typename T>
void print_array(T* array, int size){
  for(int i = 0; i < size; ++i){
    std::cout << array[i] << " ";
  }
  std::cout << std::endl;
}

template <typename T>
void print_matrix(T** matrix, int size){
  for(int i = 0; i < size; ++i){
    print_array(matrix[i],size);
  }
}

int main(int argc, char **argv){
  const int size = 5;
  std::array<double , size> matriz_diagonal{1.0,2.0,3.0,4.0,5.0};
  std::array<double , size> vector_solucion{1.0,4.0,2.0,5.0,7.0};
  std::array<double , size> vector_incognitas{};
  std::array<double , size> diagonal_inversa{};
  double matriz_inferior[5][5]{
  {1.0,0.0,0.0,0.0,0.0},
  {3.0,2.0,0.0,0.0,0.0},
  {4.0,5.0,1.0,0.0,0.0},
  {2.0,3.0,3.0,3.0,0.0},
  {7.0,2.0,5.0,2.0,6.0}
  };
  double matriz_superior[5][5]{
  {1.0,2.0,8.0,2.0,1.0},
  {0.0,3.0,1.0,5.0,5.0},
  {0.0,0.0,3.0,7.0,4.0},
  {0.0,0.0,0.0,2.0,5.0},
  {0.0,0.0,0.0,0.0,7.0}
  };
  std::cout << "matriz diagonal: ";
  print_array(matriz_diagonal.data(), matriz_diagonal.size());
  std::cout << "determinante de la matriz diagonal "
    <<determinante_diagonal( matriz_diagonal.data(), matriz_diagonal.size() )
    << std::endl;
  std::cout << "matriz diagonal inversa: ";
  inversa_diagonal(matriz_diagonal.data(), diagonal_inversa.data() ,matriz_diagonal.size());
  print_array( diagonal_inversa.data(),matriz_diagonal.size());
  std::cout << "solucion de la matriz diagonal: ";
  solucion_diagonal(matriz_diagonal.data(), vector_incognitas.data(), vector_solucion.data() , matriz_diagonal.size() );
  print_array( vector_incognitas.data(),vector_incognitas.size());
  double *ptr_matrix[size];
  matrix_to_ptr(matriz_inferior, ptr_matrix, size)
  std::cout << "matriz triangular inferior: " << std::endl;
  print_matrix(ptr_matrix, size);
  std::cout << "determinante de la triangular inferior "
    <<determinante_triangular( ptr_matrix, size )
    << std::endl; 
  std::cout << "solucion de la matriz diagonal: ";
  solucion_triangular_inf(ptr_matrix, vector_incognitas.data(), vector_solucion.data() , size );
  print_array( vector_incognitas.data(),vector_incognitas.size());
  matrix_to_ptr(matriz_superior, ptr_matrix, size)
  std::cout << "matriz triangular superior: " << std::endl;
  print_matrix(ptr_matrix, size);
  std::cout << "determinante de la triangular superior "
    <<determinante_triangular( ptr_matrix, size )
    << std::endl; 
  std::cout << "solucion de la matriz diagonal: ";
  solucion_triangular_sup(ptr_matrix, vector_incognitas.data(), vector_solucion.data() , size );
  print_array( vector_incognitas.data(),vector_incognitas.size());
}
