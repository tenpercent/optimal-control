#include "calculate_functional.hh"
#include "runge-kutta_method.hh"

double integrant_value (boundary_conditions const & bvalue) {
  return bvalue.p2 * bvalue.p2;
}

double simpson_formula (
    boundary_conditions const &   left,   // node 1
    boundary_conditions const &   middle,   // node 2
    boundary_conditions const &   right,   // node 3
    double const                  segment_length   // step
)
{
  return (segment_length / 6) * (
    integrant_value (left) + 
    integrant_value (middle) * 4 + 
    integrant_value (right)
  );
}

double calculate_functional (
    double const                  step_length,
    size_t const                  total_steps,
    boundary_conditions const &   bc_start,
    double const                  parameter
) 
{
  boundary_conditions bc_left, bc_middle, bc_right;
  double sum = 0.0;
  
  bc_left = bc_start;
  
  for (size_t i = 0; i < total_steps; ++i) {
    bc_middle = runge_kutta_method (0.5 * step_length, 1, bc_left, parameter);
    bc_right = runge_kutta_method (0.5 * step_length, 1, bc_middle, parameter);
    sum += simpson_formula (bc_left, bc_middle, bc_right, step_length);
    bc_left = bc_right;
  }
  
  return sum;
}