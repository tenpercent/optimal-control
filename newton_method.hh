#include <stdint.h>
#include "boundary_conditions.hh"

int newton_method (
    double const                  tau,  // step
    uint32_t const                n,    // step numbers
    boundary_conditions &         bc_0, // conditions in 0
    boundary_conditions &         bc_1, // conditions in 1
    double const                  alpha // parametr
);

