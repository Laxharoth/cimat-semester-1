##############################
#
# Para ejecutar el programa del primer problema utilizar make solucion
# 	una vez compilado el programa se puede ejecutar con ./solucion.run
#
# Para ejecutar el programa del segundo problema utilizar make inversa
# 	una vez compilado el programa se puede ejecutar con ./inversa.run
#
# Para ejecutar el programa del tercer problema utilizar make eigen
# 	una vez compilado el programa se puede ejecutar con ./eigen.run
#
####

real_m_src = matrix/real_matrix.hpp matrix/real_matrix.cpp matrix/real_vector.cpp matrix/real_vector_iterator.cpp
funcion_matrix_src = matrix/funcion_matriz.hpp matrix/funcion_matriz.cpp matrix/factorizacion.cpp

eigen: funcion_matriz.o real_matrix.o
	g++ -pthread -fopenmp -o eigen.run main_eigen.cpp funcion_matriz.o real_matrix.o
	./eigen.run

inversa: funcion_matriz.o real_matrix.o
	g++ -pthread -fopenmp -o inversa.run main_inversa.cpp funcion_matriz.o real_matrix.o
	./inversa.run

solucion: funcion_matriz.o real_matrix.o
	g++ -pthread -fopenmp -o solucion.run main_solucion.cpp funcion_matriz.o real_matrix.o
	./solucion.run
	
funcion_matriz.o: $(funcion_matrix_src)
	$(MAKE) -C matrix funcion_matriz.o
	mv matrix/funcion_matriz.o .

real_matrix.o : $(real_m_src)
	$(MAKE) -C matrix real_matrix.o
	mv matrix/real_matrix.o .

clean:
	rm -f *.o *.run
