object_files =gauss_jordan.o funcion_matriz.o matriz_banda.o factorizacion.o

programa.run: main.o $(object_files)
	g++ -o programa.run main.o $(object_files)
main.o: main.cpp funcion_matriz.hpp gauss_jordan.hpp factorizacion.hpp matriz_banda.hpp
	g++ -c main.cpp
gauss_jordan.o: gauss_jordan.hpp funcion_matriz.hpp gauss_jordan.cpp
	g++ -c gauss_jordan.cpp
funcion_matriz.o: funcion_matriz.hpp funcion_matriz.cpp
	g++ -c funcion_matriz.cpp
factorizacion.o: factorizacion.hpp factorizacion.cpp
	g++ -c factorizacion.cpp
matriz_banda.o: matriz_banda.hpp matriz_banda.cpp
	g++ -c matriz_banda.cpp
test: test.run
	clear 
	./test.run
test.run: test.o $(object_files)
	g++ -o test.run test.o $(object_files)
test.o: test.cpp funcion_matriz.hpp gauss_jordan.hpp factorizacion.hpp matriz_banda.hpp
	g++ -c test.cpp 
clean:
	rm *.o

