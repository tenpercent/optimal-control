CC = g++ $(CFLAGS)
CFLAGS = -O2
OBJFILES = main.o newton_method.o runge-kutta_method.o calculate_functional.o boundary_conditions.o

all: optcontrol
optcontrol: $(OBJFILES)
	$(CC) -o optcontrol $(OBJFILES)
main.o: main.cc defs.hh newton_method.hh calculate_functional.hh boundary_conditions.hh
	$(CC) -c main.cc
newton_method.o: newton_method.cc newton_method.hh defs.hh boundary_conditions.hh
	$(CC) -c newton_method.cc
runge-kutta_method.o: runge-kutta_method.cc runge-kutta_method.hh boundary_conditions.hh
	$(CC) -c runge-kutta_method.cc
calculate_functional.o: calculate_functional.cc calculate_functional.hh boundary_conditions.hh
	$(CC) -c calculate_functional.cc
boundary_conditions.o: boundary_conditions.cc boundary_conditions.hh
	$(CC) -c boundary_conditions.cc
clean:
	rm -f optcontrol *.o *.cc~ *.hh~
