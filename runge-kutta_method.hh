#include "boundary_conditions.hh"

typedef unsigned long size_t;

boundary_conditions get_ode_system_right_part_value (
    const boundary_conditions & a, 
    double const alpha
);

boundary_conditions runge_kutta_method (
    double const                  tau,  // step
    size_t const                  n,    // step numbers
    boundary_conditions const &   bc_0, // initial boundary conditions
    double const                  alpha // parametr
);