matrix_exception.o: matrix_exception.hpp matrix_exception.cpp 
	g++ -c matrix_exception.cpp

matriz_banda.o: matriz_banda.hpp matriz_banda.cpp matrix_like.tcc matrix_iterator.tcc matrix_exception.hpp
	g++ -c matriz_banda.cpp

real_matrix.o: real_matrix.hpp real_matrix.cpp  matrix_like.tcc
	g++ -c real_matrix.cpp

test: test.run
test.run: test.o real_matrix.o matrix_exception.o
	g++ -o test.run test.o real_matrix.o
test.o: test.cpp real_matrix.hpp matrix_exception.hpp matrix_like.tcc
	g++ -c test.cpp