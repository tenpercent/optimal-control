#include "boundary_conditions.hh"

int newton_method (
    double const                  tau,
    size_t const                  n,
    boundary_conditions &         bc_0,
    boundary_conditions &         bc_1,
    double const                  alpha
);

