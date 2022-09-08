object_files = real_matrix.o funcion_matriz.o factorizacion.o

run: programa.run
	./programa.run
programa.run: main.o $(object_files)
	g++ -O2 -o programa.run main.o $(object_files)
main.o: main.cpp matrix_like/real_matrix.hpp funcion_matriz.hpp macros.hpp
	g++ -O2 -c main.cpp
test: test.run
	./test.run
test.run: test.o $(object_files)
	g++ -O2 -o test.run test.o $(object_files)
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