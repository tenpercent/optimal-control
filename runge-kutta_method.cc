#include "runge-kutta_method.hh"
#include <cmath>

boundary_conditions get_ode_system_right_part_value (
    boundary_conditions const & bc_value,
    double const parameter)
{
  boundary_conditions answer;

  answer.x1 = bc_value.x2;
  answer.x2 = bc_value.p2 + bc_value.x1 * cos (parameter * bc_value.x2);
  answer.p1 = -bc_value.p2 * cos (parameter * bc_value.x2);
  answer.p2 = -bc_value.p1 + parameter * bc_value.p2 * bc_value.x1 * sin (parameter * bc_value.x2);

  return answer;
}

boundary_conditions runge_kutta_method (
    double const                    step_length,  // step length
    uint32_t const                  total_steps,  // total steps
    boundary_conditions const &     bc_initial, // initial boundary conditions
    double const                    parameter // parameter
)
{
  boundary_conditions q1, q2, q3, q4, answer;
  answer = bc_initial;

  for (uint32_t i = 0; i < total_steps; ++i) {
    q1 = get_ode_system_right_part_value (answer, parameter)           * step_length;
    q2 = get_ode_system_right_part_value (answer + q1 * .5, parameter) * step_length;
    q3 = get_ode_system_right_part_value (answer + q2 * .5, parameter) * step_length;
    q4 = get_ode_system_right_part_value (answer + q3, parameter)      * step_length;
    answer += (q1 + q2 * 2 + q3 * 2 + q4) / 6;
  }
  return answer;
}
