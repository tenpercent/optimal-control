CC = g++ $(CFLAGS)
CFLAGS = -O2

all: task2

task2: 
	$(CC) -o task2 \
  main.cc \
  newton_method.cc \
  runge-kutta_method.cc \
  boundary_conditions.cc \
  calculate_functional.cc
clean:
	rm -f task2 *.o *.cc~ *.hh~
