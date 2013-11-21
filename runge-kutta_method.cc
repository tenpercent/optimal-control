#include "runge-kutta_method.hh"
#include <cmath>

boundary_conditions get_ode_system_right_part_value (const boundary_conditions & a, double const alpha) {
  boundary_conditions ret;
  
  ret.x1 = a.x2;
  ret.x2 = a.p2 + a.x1 * cos (alpha * a.x2);
  ret.p1 = -a.p2 * cos (alpha * a.x2);
  ret.p2 = -a.p1 + alpha * a.p2 * a.x1 * sin (alpha * a.x2);
  
  return ret;
}

boundary_conditions runge_kutta_method (
    double const                    tau,  // step length
    size_t const                    n,    // total steps
    boundary_conditions const &     bc_0, // initial boundary conditions
    double const                    alpha // parameter
)
{
  boundary_conditions q1, q2, q3, q4, answer;  
  answer = bc_0;
  
  for (size_t i = 0; i < n; ++i) {
    q1 = get_ode_system_right_part_value (answer, alpha)           * tau;  
    q2 = get_ode_system_right_part_value (answer + q1 * .5, alpha) * tau;
    q3 = get_ode_system_right_part_value (answer + q2 * .5, alpha) * tau;  
    q4 = get_ode_system_right_part_value (answer + q3, alpha)      * tau;   
    answer += (q1 + q2 * 2 + q3 * 2 + q4) / 6;
  }
  return answer;
}