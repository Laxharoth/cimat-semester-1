object_files = real_matrix.o function_wrapper.o funcion_matriz.o factorizacion.o 

main: main.run
	./main.run
main.run :main_interpolacion.cpp interpolation.cpp $(object_files)
	g++ -O2 -pthread -g -o main.run main_interpolacion.cpp interpolation.cpp $(object_files)
test: test.run
	./test.run
test.run :test.cpp interpolation.cpp $(object_files)
	g++ -O2 -pthread -g -o test.run test.cpp interpolation.cpp $(object_files)
real_matrix.o: matrix_like/real_matrix.hpp matrix_like/real_matrix.cpp matrix_like/real_vector.cpp matrix_like/real_vector_iterator.cpp
	$(MAKE) -C matrix_like real_matrix.o
	mv matrix_like/real_matrix.o .
function_wrapper.o: function_wrapper/function_wrapper.hpp function_wrapper/function_wrapper.cpp
	$(MAKE) -C function_wrapper function_wrapper.o
	mv function_wrapper/function_wrapper.o .
interpolation.o: interpolation.hpp interpolation.cpp matrix_like/real_matrix.hpp function_wrapper/function_wrapper.hpp funcion_matriz.hpp
	g++ -O2  -c interpolation.cpp
funcion_matriz.o: funcion_matriz.hpp funcion_matriz.cpp matrix_like/real_matrix.hpp factorizacion.hpp
	g++ -O2 -c funcion_matriz.cpp
factorizacion.o: factorizacion.hpp factorizacion.cpp matrix_like/real_matrix.hpp
	g++ -O2 -c factorizacion.cpp