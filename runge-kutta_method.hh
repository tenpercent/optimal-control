#include <stdint.h>
#include "boundary_conditions.hh"

boundary_conditions get_ode_system_right_part_value (
    boundary_conditions const & bc_value, 
    double const parameter
);

boundary_conditions runge_kutta_method (
    double const                    step_length,
    uint32_t const                  total_steps,
    boundary_conditions const &     bc_initial, 
    double const                    parameter 
);
