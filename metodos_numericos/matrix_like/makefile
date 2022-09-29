matrix_exception.o: matrix_exception.hpp matrix_exception.cpp 
	g++ -O2 -c matrix_exception.cpp

matriz_banda.o: matriz_banda.hpp matriz_banda.cpp matrix_like.tcc matrix_iterator.tcc matrix_exception.hpp
	g++ -O2 -c matriz_banda.cpp

real_matrix.o: real_matrix.hpp real_matrix.cpp real_vector.cpp real_vector_iterator.cpp
	g++ -O2 -c real_matrix.cpp

test: test.run
	./test.run
test.run: test.cpp real_matrix.cpp
	g++ -g -pthread -o test.run test.cpp real_matrix.cpp
clean:
	rm -rf *.o
	rm -rf *.run