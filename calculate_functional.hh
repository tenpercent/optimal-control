#include <stdint.h>
#include "boundary_conditions.hh"

double integrant_value (boundary_conditions const & b1);

double simpson_formula (
    boundary_conditions const &   left,   // node 1
    boundary_conditions const &   middle,   // node 2
    boundary_conditions const &   right,   // node 3
    double const                  segment_length   // step
);

double calculate_functional (
    double const                  step_length,  // step
    uint32_t const                total_steps,  // step numbers
    boundary_conditions const &   bc_start, // start conditions
    double const                  parameter // parameter
);
