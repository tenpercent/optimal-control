#include "runge-kutta_method.hh"
#include "newton_method.hh"
#include "defs.hh"

double calc_residual (boundary_conditions const & calculated, boundary_conditions const & wanted, double const * jacobi_matrix) {
  double const k1 = jacobi_matrix[0] * jacobi_matrix[0] + jacobi_matrix[1] * jacobi_matrix[1];
  double const k2 = jacobi_matrix[2] * jacobi_matrix[2] + jacobi_matrix[3] * jacobi_matrix[3];
  return sqrt (
      (calculated.x1 - wanted.x1) * (calculated.x1 - wanted.x1) / k1 + 
      (calculated.x2 - wanted.x2) * (calculated.x2 - wanted.x2) / k2
  );
}

boundary_conditions get_bc_at_right_end (
  double const step_length, 
  size_t const total_steps, 
  boundary_conditions const & bc_0, 
  double const parameter) 
{
  boundary_conditions start = bc_0;
  return integrator (step_length, total_steps, start, parameter);
}

void get_jacobi_matrix (
  double const step_length, 
  size_t const total_steps, 
  boundary_conditions const & bc_0, 
  boundary_conditions const & bc_at_right_end, 
  double const parameter, 
  double * jacobi_matrix
  ) {
  boundary_conditions start = bc_0;
  start.p1 = bc_0.p1 + DELTA;

  boundary_conditions bc_at_right_end_moved = integrator (step_length, total_steps, start, parameter);

  jacobi_matrix[0] = (bc_at_right_end_moved.x1 - bc_at_right_end.x1) / DELTA;
  jacobi_matrix[2] = (bc_at_right_end_moved.x2 - bc_at_right_end.x2) / DELTA;

  start.p1 = bc_0.p1;
  start.p2 = bc_0.p2 + DELTA;
  
  bc_at_right_end_moved = integrator (step_length, total_steps, start, parameter);

  jacobi_matrix[1] = (bc_at_right_end_moved.x1 - bc_at_right_end.x1) / DELTA;
  jacobi_matrix[3] = (bc_at_right_end_moved.x2 - bc_at_right_end.x2) / DELTA;
}

boundary_conditions get_next (
  boundary_conditions const & bc_0, 
  boundary_conditions const & bc_1, 
  boundary_conditions const & bc_at_right_end, 
  boundary_conditions const & start,
  double const gamma, 
  double const * jacobi_matrix
) {
  double jacobi_det = jacobi_matrix[0] * jacobi_matrix[3] - jacobi_matrix[1] * jacobi_matrix[2];
  boundary_conditions next = start;

  next.p1 = bc_0.p1 - gamma * (
    (bc_at_right_end.x1 - bc_1.x1) * jacobi_matrix[3] - 
    (bc_at_right_end.x2 - bc_1.x2) * jacobi_matrix[1]
    ) 
  / jacobi_det;

  next.p2 = bc_0.p2 - gamma * (
    (bc_at_right_end.x2 - bc_1.x2) * jacobi_matrix[0] - 
    (bc_at_right_end.x1 - bc_1.x1) * jacobi_matrix[2]
    ) 
  / jacobi_det;
  return next;
}

int newton_method (
    double const                 step_length,    // step length
    size_t const                 total_steps,    // step number
    boundary_conditions &        bc_0,           // conditions in 0
    boundary_conditions &        bc_1,           // conditions in 1
    double const                 parameter       // parameter
)
{
  boundary_conditions start, next, bc_at_right_end, bc_at_right_end_moved;
  size_t iterations = 0;
  double gamma, residual, next_residual;
  double jacobi_matrix[4] = {0};
  
  start = bc_0;
  
  while (iterations < MAX_ITERATIONS) {
    gamma = 1.0;
    bc_at_right_end = get_bc_at_right_end(step_length, total_steps, bc_0, parameter);
    get_jacobi_matrix (step_length, total_steps, bc_0, bc_at_right_end, parameter, jacobi_matrix);
    residual = calc_residual (bc_at_right_end, bc_1, jacobi_matrix);
    if (residual < EPS) {
      bc_1.p1 = bc_at_right_end.p1;
      bc_1.p2 = bc_at_right_end.p2;
      return 0;
    }
      
    next_residual = residual + 1;
    while (next_residual > residual) {
      next = get_next (bc_0, bc_1, bc_at_right_end, start, gamma, jacobi_matrix);
      bc_at_right_end = get_bc_at_right_end (step_length, total_steps, next, parameter);
      get_jacobi_matrix (step_length, total_steps, next, bc_at_right_end, parameter, jacobi_matrix);
      next_residual = calc_residual (bc_at_right_end, bc_1, jacobi_matrix);
      if (next_residual < EPS) {
        bc_0.p1 = next.p1;
        bc_0.p2 = next.p2;
        bc_1.p1 = bc_at_right_end.p1;
        bc_1.p2 = bc_at_right_end.p2;
        return 0;
      }
      gamma *= .5;
    }
    bc_0.p1 = next.p1;
    bc_0.p2 = next.p2;
    bc_1.p1 = bc_at_right_end.p1;
    bc_1.p2 = bc_at_right_end.p2;
    residual = next_residual;
    ++iterations;
  }  
  return 1;
}
