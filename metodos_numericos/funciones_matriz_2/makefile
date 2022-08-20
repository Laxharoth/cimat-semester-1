object_files =gauss_jordan.o funcion_matriz.o matriz_banda.o factorizacion.o

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
test.o: test.cpp funcion_matriz.hpp gauss_jordan.hpp factorizacion.hpp matrix_like/matriz_banda.hpp macros.hpp
	g++ -c test.cpp 
clean:
	rm *.o
matriz_banda.o: matrix_like/matriz_banda.hpp matrix_like/matriz_banda.cpp
	$(MAKE) matrix_like
	mv matrix_like/matriz_banda.o .
