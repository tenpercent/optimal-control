#include <cstdio>

#include "newton_method.hh"
#include "calculate_functional.hh"
#include "defs.hh"

int main () {
  double parameter;
  printf ("Enter alpha\n");
  scanf ("%lf", &parameter);
  
  size_t total_steps;
  printf ("Enter number of steps\n");
  scanf ("%ld", &total_steps);
  
  double step_length = 1. / total_steps;
  printf ("Step of Runge-Kutta algorithm = %.14le\n", step_length);
  
  boundary_conditions bc_in_0 (X1_0, X2_0, P1_0, P2_0); // look in defs.hh
  boundary_conditions bc_in_1 (X1_1, X2_1, 0, 0);
  
  if (newton_method (step_length, total_steps, bc_in_0, bc_in_1, parameter)) {
    printf ("Acceleration of Newton algorithm broke down.\n");
    return 1;
  }
  
  printf ("Вoundary conditions:\n");
  printf ("\tx1(0) = %.4lf\tx1(1) = %.4lf\n", bc_in_0.x1, bc_in_1.x1);
  printf ("\tx2(0) = %.4lf\tx2(1) = %.4lf\n", bc_in_0.x2, bc_in_1.x2);
  printf ("\tp1(0) = %.4lf\tp1(1) = %.4lf\n", bc_in_0.p1, bc_in_1.p1);
  printf ("\tp2(0) = %.4lf\tp2(1) = %.4lf\n", bc_in_0.p2, bc_in_1.p2);
  
  double result = calculate_functional (step_length, total_steps, bc_in_0, parameter);
  
  printf ("answer = %.14le\n", result);
  
  return 0;
}
