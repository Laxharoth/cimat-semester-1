object_files = main.o funcion_matriz.o

exec: programa.run
	./programa.run

programa.run: $(object_files)
	g++ -o programa.run $(object_files)

main.o: main.cpp funcion_matriz.hpp
	g++ -c main.cpp

funcion_matriz.o: funcion_matriz.hpp funcion_matriz.cpp
	g++ -c funcion_matriz.cpp