#include <cstdio>
#include <cmath>

#include "newton_method.hh"
#include "calculate_functional.hh"
#include "defs.hh"

static double MIN_DOUBLE = 1e-15;

double get_alpha (void);
double get_step_length (void);
void set_initial_boundary_conditions (boundary_conditions & bc_in_0, boundary_conditions & bc_in_1);

int main () {
  double const parameter = get_alpha ();
  double const step_length = get_step_length ();

  if (step_length < MIN_DOUBLE) {
    printf ("Incorrect step length\n");
    return 2;
  }

  size_t const runge_kutta_steps = 1. / step_length;
  EPS = 10 * step_length * step_length * step_length * step_length; 
  DELTA = 100 * EPS;
  
  boundary_conditions bc_in_0, bc_in_1;
  set_initial_boundary_conditions (bc_in_0, bc_in_1);
  
  if (newton_method (step_length, runge_kutta_steps, bc_in_0, bc_in_1, parameter)) {
    printf ("Newton algorithm acceleration failed.\n");
    return 1;
  }

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

double get_step_length (void) {
  double precision;
  printf ("Enter Runge-Kutta method step length\n");
  scanf ("%lf", &precision);
  return precision;
}

void set_initial_boundary_conditions (boundary_conditions & bc_in_0, boundary_conditions & bc_in_1) {
  bc_in_0 = boundary_conditions (X1_0, X2_0, P1_0, P2_0); // look in defs.hh
  bc_in_1 = boundary_conditions (X1_1, X2_1, 0, 0);       // look in defs.hh
}
