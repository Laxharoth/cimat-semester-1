##############################
# Para Ejecutar el programa  #
#  utilizar el comando make  #
##############################
objetos = main.o function_wrapper.o fparser.o

exec: programa.run
	./programa.run x^3 -1 5
programa.run: $(objetos)
	g++ -o programa.run $(objetos)
main.o: main.cpp biseccion.cpp newton.cpp search_min_max.cpp function_wrapper/parser_wrapper.hpp fparser/fparser.hh
	g++ -c main.cpp
function_wrapper.o:
	$(MAKE) -C function_wrapper
	mv function_wrapper/function_wrapper.o .
fparser.o:
	$(MAKE) -C fparser
	mv fparser/fparser.o .
clean:
	rm *.o
	rm *.run
