###################################################
#
# Para ejecutar las pruebas de los metodos
# utilizar el comando make subspace
#
# La salida en texto se guarda en el archivo subspace_out.txt
#
# Para obtener los valores esperados se utilizo numpy en el notebook subspace.ipynb
# Para comparar los valores esperados con los obtenidos se utilizo numpy en el notebook subspace_cmp.ipynb
#
# Las funciones implementadas se encuentran el el archivo funcion_matrix.cpp a partir de la linea 318
#
# Los vectores y matrices se guardan en el directorio vector.out
#
# Los vectores se guandan en archivos binario,
# 	Los primeros 8 bytes corresponden a un entero sin signo el cual es el tamaño del vector
#	El resto del archivo representa los valores del vector y estan guardados en paquetes
#		de 8 bytes como flotantes de doble precision
#
# Las matrices seguardan en archivos binarios,
#	Los primeros 8 bytes son un entero sin signo que corresponde al numero de filas
#	Los siguentes 8 bytes son un entero sin signo que corresponde al numero de columnas
#	El resto del archivo representa los valores de la matriz y estan guardados por filas
#		en paquetes de 8 bytes flotantes de doble precision
#
###############################
object_files = real_matrix.o funcion_matriz.o

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
subspace: subspace.run
	./subspace.run
subspace.run: subspace_main.o $(object_files)
	g++ -O2 -pthread -o subspace.run subspace_main.o $(object_files)
subspace_main.o: subspace_main.cpp matrix_like/real_matrix.hpp funcion_matriz.hpp macros.hpp
	g++ -O2 -c subspace_main.cpp
main.o: main.cpp matrix_like/real_matrix.hpp funcion_matriz.hpp macros.hpp
	g++ -O2 -c main.cpp
test: test.run
	./test.run
test.run: test.o $(object_files)
	g++ -O2 -pthread -o test.run test.o $(object_files)
test.o: test.cpp matrix_like/real_matrix.hpp
	g++ -O2 -c test.cpp
funcion_matriz.o: funcion_matriz.hpp funcion_matriz.cpp factorizacion.cpp matrix_like/real_matrix.hpp
	g++ -O2 -c funcion_matriz.cpp
real_matrix.o: matrix_like/real_matrix.hpp matrix_like/real_matrix.cpp matrix_like/real_vector.cpp matrix_like/real_vector_iterator.cpp
	$(MAKE) -C matrix_like real_matrix.o
	mv matrix_like/real_matrix.o .
clean:
	rm -rf *.o
	rm -rf *.run