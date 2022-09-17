##############################
#
#		Para ejecutar el programa utilizado para
#   realizar las pruebas para que asegurarse 
#   de que las funciones funcionen como esperado
# 	ejecutar 'make test'
#  		Para ejecutar el programa que resuelve la
#	ecuacionde calor usando cholesky ejecutar 'make run'
#  		Para ejecutar el programa que resuelve la
#	ecuacionde calor usando jacobi ejecutar 'make jacob'
#  		Para ejecutar el programa que resuelve la
#	ecuacionde calor usando gauss-seidel ejecutar 'make gauss'
#
##############################
object_files = functions.o funcion_matriz.o matrix_exception.o factorizacion.o real_matrix.o
matrix_like_dep = matrix_like/matrix_like.tcc matrix_like/matrix_iterator.tcc matrix_like/matrix_exception.hpp

run: programa.run
	./programa.run
mains: programa.run gauss.run jacob.run
	@echo e
programa.run: main.o $(object_files)
	g++ -O2 -pthread -o programa.run main.o $(object_files)
gauss: gauss.run
	./gauss.run
gauss.run: maing.o $(object_files)
	g++ -O2 -pthread -o gauss.run maing.o $(object_files)
jacob: jacob.run
	./jacob.run
jacob.run: mainj.o $(object_files)
	g++ -O2 -pthread -o jacob.run mainj.o $(object_files)
test: test.run
	./test.run
test.run: test.o $(object_files)
	g++ -O2 -pthread -o test.run test.o $(object_files)
test.o: test.cpp macros.hpp $(matrix_like_dep)
	g++ -O2 -c test.cpp
funcion_matriz.o: funcion_matriz.hpp funcion_matriz.cpp
	g++ -O2 -c funcion_matriz.cpp
factorizacion.o: factorizacion.hpp factorizacion.cpp
	g++ -O2 -c factorizacion.cpp
matrix_exception.o: matrix_like/matrix_exception.hpp matrix_like/matrix_exception.cpp
	$(MAKE) -C matrix_like matrix_exception.o
	mv matrix_like/matrix_exception.o .
real_matrix.o: matrix_like/real_matrix.hpp matrix_like/real_matrix.cpp  matrix_like/matrix_like.tcc
	$(MAKE) -C matrix_like real_matrix.o
	mv matrix_like/real_matrix.o .
clean:
	rm *.o &
	rm *.run &
	rm -rf __pycache__ &
	rm -rf *.py &
	rm -rf *.txt &