#include "boundary_conditions.hh"

int newton_method (
    double const                  step_length,
    size_t const                  total_steps,
    boundary_conditions &         bc_0,
    boundary_conditions &         bc_1,
    double const                  parameter
);
