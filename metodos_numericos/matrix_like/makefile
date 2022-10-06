matrix.o: matrix.hpp matrix.cpp vector.cpp vector_iterator.cpp
	g++ -O2 -c matrix.cpp

funcion_matriz.o: funcion_matriz.hpp funcion_matriz.cpp matrix.o
	g++ -O2 -c funcion_matriz.cpp

test: test.run
	./test.run
test.run: test.cpp real_matrix.cpp real_vector.cpp real_vector_iterator.cpp
	g++ -g -pthread -o test.run test.cpp real_matrix.cpp
clean:
	rm -rf *.o
	rm -rf *.run