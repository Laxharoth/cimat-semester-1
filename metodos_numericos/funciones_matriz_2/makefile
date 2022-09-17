#******************************************#
#   Para ejecutar las pruebas hechas       #
#   como comprobación de la funcionalidad  #
#   utilizar "make test"                   #
#   Para ejecutar las pruebas hechas       #
#   para la matriz en el archivo A.txt     #
#   utilizar "make"                        #
#******************************************#
#	la funcion para encontrar maximos y minimos
# se encuentra en el archivo search_min_max.cpp
# 	las funciones para resolver ecuaciones 
# se encuentran en los arhivos funcion_matriz.*
# 	las funcions para factorizar matrices
# se encuentran en los archivos factorizacion.*
############################################
object_files =funcion_matriz.o matriz_banda.o factorizacion.o matrix_exception.o real_matrix.o

 matrix_like = matrix_like/matrix_like.tcc matrix_like/matrix_iterator.tcc

programa: programa.run
	./programa.run
programa.run: main.o $(object_files)
	g++ -O2 -pthread -o programa.run main.o $(object_files)
main.o: main.cpp funcion_matriz.hpp factorizacion.hpp matrix_like/matriz_banda.hpp
	g++ -O2 -c main.cpp
funcion_matriz.o: funcion_matriz.hpp funcion_matriz.cpp
	g++ -O2 -c funcion_matriz.cpp
factorizacion.o: factorizacion.hpp factorizacion.cpp
	g++ -O2 -c factorizacion.cpp
test: test.run
	clear 
	./test.run
test.run: test.o $(object_files)
	g++ -O2 -pthread  -o test.run test.o $(object_files)
test.o: test.cpp funcion_matriz.hpp factorizacion.hpp $(matrix_like) matrix_like/matriz_banda.hpp matrix_like/matriz_banda.cpp macros.hpp function_wrapper/function_wrapper.hpp search_min_max.cpp
	g++ -O2 -c test.cpp 
clean:
	rm *.o &
	rm *.run &
matrix_exception.o: matrix_like/matrix_exception.hpp matrix_like/matrix_exception.cpp
	$(MAKE) -C matrix_like matrix_exception.o
	mv matrix_like/matrix_exception.o .
matriz_banda.o: matrix_like/matriz_banda.hpp matrix_like/matriz_banda.cpp $(matrix_like)
	$(MAKE) -C matrix_like matriz_banda.o
	mv matrix_like/matriz_banda.o .
real_matrix.o: matrix_like/real_matrix.hpp matrix_like/real_matrix.cpp  matrix_like/matrix_like.tcc
	$(MAKE) -C matrix_like real_matrix.o
	mv matrix_like/real_matrix.o .
