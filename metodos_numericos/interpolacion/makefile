##################################
#
# El programa realiza integrales numericas utilizando el método
#   de Cuadratura Gaussiana
#   para ejecutar el programa utilizar el comando `make int`
#
###################################

object_files = matrix.o function_wrapper.o funcion_matriz.o

int: int.run
	@./int.run
int.run :main_integration.cpp interpolation.cpp $(object_files)
	g++ -Og -pthread -g -o int.run main_integration.cpp interpolation.cpp $(object_files)
int2: int2.run
	./int2.run
int2.run :interpolation_finite2.cpp interpolation.cpp $(object_files)
	g++ -Og -pthread -g -o int2.run interpolation_finite2.cpp interpolation.cpp $(object_files)
main: main.run
	./main.run
main.run :main_interpolacion.cpp interpolation.cpp $(object_files)
	g++ -O2 -pthread -g -o main.run main_interpolacion.cpp interpolation.cpp $(object_files)
main2: main2.run
	./main2.run
main2.run :main_interpolacion2.cpp interpolation.cpp $(object_files)
	g++ -O2 -pthread -g -o main2.run main_interpolacion2.cpp interpolation.cpp $(object_files)
main3: main3.run
	./main3.run
main3.run :main_interpolacion3.cpp interpolation.cpp $(object_files)
	g++ -O2 -pthread -g -o main3.run main_interpolacion3.cpp interpolation.cpp $(object_files) -lm
test: test.run
	./test.run
test.run :test.cpp interpolation.cpp $(object_files)
	g++ -O2 -pthread -g -o test.run test.cpp interpolation.cpp $(object_files)
matrix.o: matrix_like/matrix.hpp matrix_like/matrix.cpp matrix_like/vector.cpp matrix_like/vector_iterator.cpp
	$(MAKE) -C matrix_like matrix.o
	mv matrix_like/matrix.o .
funcion_matriz.o: matrix_like/funcion_matriz.hpp matrix_like/funcion_matriz.cpp
	$(MAKE) -C matrix_like funcion_matriz.o
	mv matrix_like/funcion_matriz.o .
function_wrapper.o: function_wrapper/function_wrapper.hpp function_wrapper/function_wrapper.cpp
	$(MAKE) -C function_wrapper function_wrapper.o
	mv function_wrapper/function_wrapper.o .
interpolation.o: interpolation.hpp interpolation.cpp matrix_like/real_matrix.hpp function_wrapper/function_wrapper.hpp funcion_matriz.hpp
	g++ -O2  -c interpolation.cpp
clean:
	rm -f *.o
	rm -f *.run