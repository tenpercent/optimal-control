#include "boundary_conditions.hh"

typedef unsigned long size_t;

double integrant_value (boundary_conditions const & b1);

double simpson_formula (
    boundary_conditions const &   left,   // node 1
    boundary_conditions const &   middle,   // node 2
    boundary_conditions const &   right,   // node 3
    double const                  segment_length   // step
);

double calculate_functional (
    double const                  step_length,
    size_t const                  total_steps,
    boundary_conditions const &   bc_start,
    double const                  parameter
);
