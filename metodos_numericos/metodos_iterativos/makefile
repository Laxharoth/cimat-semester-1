object_files = functions.o funcion_matriz.o matrix_exception.o factorizacion.o real_matrix.o
matrix_like_dep = matrix_like/matrix_like.tcc matrix_like/matrix_iterator.tcc matrix_like/matrix_exception.hpp

test: test.run
	./test.run
test.run: test.o $(object_files)
	g++ -o test.run test.o $(object_files)
test.o: test.cpp macros.hpp $(matrix_like_dep)
	g++ -c test.cpp
funcion_matriz.o: funcion_matriz.hpp funcion_matriz.cpp
	g++ -c funcion_matriz.cpp
factorizacion.o: factorizacion.hpp factorizacion.cpp
	g++ -c factorizacion.cpp
matrix_exception.o: matrix_like/matrix_exception.hpp matrix_like/matrix_exception.cpp
	$(MAKE) -C matrix_like matrix_exception.o
	mv matrix_like/matrix_exception.o .
real_matrix.o: matrix_like/real_matrix.hpp matrix_like/real_matrix.cpp  matrix_like/matrix_like.tcc
	$(MAKE) -C matrix_like real_matrix.o
	mv matrix_like/real_matrix.o .