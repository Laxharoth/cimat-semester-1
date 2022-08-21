object_files =gauss_jordan.o funcion_matriz.o matriz_banda.o factorizacion.o matrix_exception.o

programa.run: main.o $(object_files)
	g++ -o programa.run main.o $(object_files)
main.o: main.cpp funcion_matriz.hpp gauss_jordan.hpp factorizacion.hpp matrix_like/matriz_banda.hpp
	g++ -c main.cpp
gauss_jordan.o: gauss_jordan.hpp funcion_matriz.hpp gauss_jordan.cpp
	g++ -c gauss_jordan.cpp
funcion_matriz.o: funcion_matriz.hpp funcion_matriz.cpp
	g++ -c funcion_matriz.cpp
factorizacion.o: factorizacion.hpp factorizacion.cpp
	g++ -c factorizacion.cpp
test: test.run
	clear 
	./test.run
test.run: test.o $(object_files)
	g++ -o test.run test.o $(object_files)
test.o: test.cpp funcion_matriz.hpp gauss_jordan.hpp factorizacion.hpp matrix_like/matrix_like.tcc matrix_like/matriz_banda.hpp matrix_like/matriz_banda.cpp macros.hpp
	g++ -c test.cpp 
clean:
	rm *.o
matrix_exception.o: matrix_like/matrix_exception.hpp matrix_like/matrix_exception.cpp
	$(MAKE) -C matrix_like matrix_exception.o
	mv matrix_like/matrix_exception.o .
matriz_banda.o: matrix_like/matriz_banda.hpp matrix_like/matriz_banda.cpp matrix_like/matrix_like.tcc
	$(MAKE) -C matrix_like matriz_banda.o
	mv matrix_like/matriz_banda.o .
