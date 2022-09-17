###################################################
#
# Para ejecutar el programa con el metodo de jacobi
#  se ejecuta el comando 'make jacobi'
#  Es recomendable modificar el valor asignado al macro SIZE
#   en el archivo jacobi_main.cpp, 
#   el cual corresponde al tamaño de la matriz
#
# El método de jacobi se encuentra en el archivo
#  funcion_matriz.cpp en la linea 263
#
# La comparación entre valores y vectores propios 
#  se realiza en el notebook npeig.ipynb y requiere numpy
#
###############################
object_files = real_matrix.o funcion_matriz.o factorizacion.o

run: programa.run
	./programa.run
programa.run: main.o $(object_files)
	g++ -O2 -pthread -o programa.run main.o $(object_files)
jacobi: jacobi.run
	./jacobi.run
jacobi.run: jacobi_main.o $(object_files)
	g++ -O2 -pthread -o jacobi.run jacobi_main.o $(object_files)
jacobi_main.o: jacobi_main.cpp matrix_like/real_matrix.hpp funcion_matriz.hpp macros.hpp
	g++ -O2 -c jacobi_main.cpp
bathe: bathe.run
	./bathe.run
bathe.run: bathe_main.o $(object_files)
	g++ -O2 -pthread -o bathe.run bathe_main.o $(object_files)
bathe_main.o: bathe_main.cpp matrix_like/real_matrix.hpp funcion_matriz.hpp macros.hpp
	g++ -O2 -c bathe_main.cpp
main.o: main.cpp matrix_like/real_matrix.hpp funcion_matriz.hpp macros.hpp
	g++ -O2 -c main.cpp
test: test.run
	./test.run
test.run: test.o $(object_files)
	g++ -O2 -pthread -o test.run test.o $(object_files)
test.o: test.cpp matrix_like/real_matrix.hpp
	g++ -O2 -c test.cpp
funcion_matriz.o: funcion_matriz.hpp funcion_matriz.cpp matrix_like/real_matrix.hpp
	g++ -O2 -c funcion_matriz.cpp
real_matrix.o: matrix_like/real_matrix.hpp matrix_like/real_matrix.cpp matrix_like/real_vector.cpp matrix_like/real_vector_iterator.cpp
	$(MAKE) -C matrix_like real_matrix.o
	mv matrix_like/real_matrix.o .
factorizacion.o: factorizacion.hpp factorizacion.cpp
	g++ -O2 -c factorizacion.cpp
clean:
	rm -rf *.o
	rm -rf *.out
	rm -rf *.run