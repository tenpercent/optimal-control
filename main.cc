#include <cstdio>

#include "newton_method.hh"
#include "calculate_functional.hh"
#include "defs.hh"

void print_boundary_conditions (
    boundary_conditions const & bc_in_0, 
    boundary_conditions const & bc_in_1
);
double get_alpha (void);
double get_precision (void);
void set_initial_boundary_conditions (boundary_conditions & bc_in_0, boundary_conditions & bc_in_1);


int main () {
  double const parameter = get_alpha ();
  double const precision = get_precision ();

  size_t const runge_kutta_steps = 1000;
  double const step_length = 1. / runge_kutta_steps;
  
  boundary_conditions bc_in_0, bc_in_1;
  set_initial_boundary_conditions (bc_in_0, bc_in_1);
  
  if (newton_method (step_length, runge_kutta_steps, bc_in_0, bc_in_1, parameter)) {
    printf ("Acceleration of Newton algorithm broke down.\n");
    return 1;
  }
  
  print_boundary_conditions (bc_in_0, bc_in_1);

  double result = calculate_functional (step_length, runge_kutta_steps, bc_in_0, parameter);  
  printf ("answer = %.14le\n", result);
  
  return 0;
}

double get_alpha (void) {
  double parameter;
  printf ("Enter alpha\n");
  scanf ("%lf", &parameter);
  return parameter;
}

double get_precision (void) {
  double precision;
  printf ("Enter precision\n");
  scanf ("%lf", &precision);
  return precision;
}

void print_boundary_conditions (
    boundary_conditions const & bc_in_0, 
    boundary_conditions const & bc_in_1) 
{ 
  printf ("Вoundary conditions:\n");
  printf ("\tx1(0) = %.4lf\tx1(1) = %.4lf\n", bc_in_0.x1, bc_in_1.x1);
  printf ("\tx2(0) = %.4lf\tx2(1) = %.4lf\n", bc_in_0.x2, bc_in_1.x2);
  printf ("\tp1(0) = %.4lf\tp1(1) = %.4lf\n", bc_in_0.p1, bc_in_1.p1);
  printf ("\tp2(0) = %.4lf\tp2(1) = %.4lf\n", bc_in_0.p2, bc_in_1.p2);
}

void set_initial_boundary_conditions (boundary_conditions & bc_in_0, boundary_conditions & bc_in_1) {
  bc_in_0 = boundary_conditions (X1_0, X2_0, P1_0, P2_0); // look in defs.hh
  bc_in_1 = boundary_conditions (X1_1, X2_1, 0, 0);
}
