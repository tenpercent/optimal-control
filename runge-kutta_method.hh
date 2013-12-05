#include "boundary_conditions.hh"
#include <cstddef>

boundary_conditions get_ode_system_right_part_value (
    boundary_conditions const & bc_value, 
    double const parameter
);

boundary_conditions integrator (
    double const                    step_length,
    size_t const                    total_steps,
    boundary_conditions const &     bc_initial, 
    double const                    parameter 
);
