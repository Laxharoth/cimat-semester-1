objetos = main.o function_wrapper/function_wrapper.o fparser/fparser.o

exec: programa.run
	programa.run x^3 -1 5
programa.run: $(objetos)
	g++ -o programa.run $(objetos)
main.o: main.cpp biseccion.cpp newton.cpp search_min_max.cpp function_wrapper/parser_wrapper.hpp fparser/fparser.hh
	g++ -c main.cpp
function_wrapper.o:
	$(MAKE) -C function_wrapper
fparser.o:
	$(MAKE) -C fparser