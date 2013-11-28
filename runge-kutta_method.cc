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

boundary_conditions runge_kutta_method_4 (
    double const                    step_length,
    size_t const                    total_steps,
    boundary_conditions const &     bc_initial,
    double const                    parameter
)
{
  boundary_conditions q1, q2, q3, q4, answer;
  answer = bc_initial;
  
  for (size_t i = 0; i < total_steps; ++i) {
    q1 = get_ode_system_right_part_value (answer, parameter)           * step_length;
    q2 = get_ode_system_right_part_value (answer + q1 * .5, parameter) * step_length;
    q3 = get_ode_system_right_part_value (answer + q2 * .5, parameter) * step_length;
    q4 = get_ode_system_right_part_value (answer + q3, parameter)      * step_length;
    answer += (q1 + q2 * 2 + q3 * 2 + q4) / 6;
  }
  return answer;
}

/*
// it works, but say no to spaghetti

boundary_conditions runge_kutta_method_6 ( 
    double const                    step_length,
    size_t const                    total_steps,
    boundary_conditions const &     bc_initial,
    double const                    parameter
)
{
  boundary_conditions q1, q2, q3, q4, q5, q6, answer;
  answer = bc_initial;
  
  for (size_t i = 0; i < total_steps; ++i) {
    q1 = get_ode_system_right_part_value (answer, parameter)
      * step_length;
    q2 = get_ode_system_right_part_value (answer +
        q1 * 1./4,
      parameter)
      * step_length;
    q3 = get_ode_system_right_part_value (answer +
        q1 * 3./32 +
        q2 * 9./32,
      parameter)
      * step_length;
    q4 = get_ode_system_right_part_value (answer +
        q1 * 1932./2197 +
        q2 * -7200./2197 +
        q3 * 7296./2197,
      parameter)
      * step_length;
    q5 = get_ode_system_right_part_value (answer +
        q1 * 439./216 +
        q2 * -8 +
        q3 * 3680./513 +
        q4 * -845./4104,
      parameter)
      * step_length;
    q6 = get_ode_system_right_part_value (answer +
        q1 * -8./27 +
        q2 * 2 +
        q3 * -3544./2565 +
        q4 * 1859./4104 +
        q5 * -11./40,
      parameter)
      * step_length;
    answer += (
      q1 * 16./135 +
      q2 * 0. +
      q3 * 6656./12825 +
      q4 * 28561./56430 +
      q5 * -9./50 +
      q6 * 2./55
    );
  }
  return answer;
}
*/

boundary_conditions runge_kutta_method_6 (
    double const                    step_length,
    size_t const                    total_steps,
    boundary_conditions const &     bc_initial,
    double const                    parameter
)
{
  size_t const TOTAL_INNER_STEPS = 6;
  boundary_conditions inner_steps[TOTAL_INNER_STEPS];
  boundary_conditions inner_step_helper;
  boundary_conditions answer;

  double const sum_coef[TOTAL_INNER_STEPS] = {
    16./135, 0., 6656./12825, 28561./56430, -9./50, 2./55
  };

  double const inner_step_coef[TOTAL_INNER_STEPS - 1][TOTAL_INNER_STEPS - 1] = {
  {1./4,       0.,          0.,          0.,         0.     },
  {3./32,      9./32,       0.,          0.,         0.     },
  {1932./2197, -7200./2197, 7296./2197,  0.,         0.     },
  {439./216,   -8.,         3680./513,   -845./4104, 0.     },
  {-8./27,      2.,         -3544./2565, 1859./4104, -11./40}
  };

  answer = bc_initial;
  for (
      size_t tot_steps_cnt = 0; 
      tot_steps_cnt < total_steps; 
      ++tot_steps_cnt
      ) 
  {
    for (
      size_t in_steps_cnt = 0; 
      in_steps_cnt < TOTAL_INNER_STEPS; 
      ++in_steps_cnt) 
    {
      inner_step_helper = answer;
      for (size_t i = 0; i < in_steps_cnt; ++i) {
        inner_step_helper += inner_steps[i] * inner_step_coef[in_steps_cnt - 1][i];
      }
      inner_steps[in_steps_cnt] = get_ode_system_right_part_value (inner_step_helper, parameter) * step_length;
    }

    for (
      size_t in_steps_cnt = 0; 
      in_steps_cnt < TOTAL_INNER_STEPS; 
      ++in_steps_cnt) 
    {
      answer += inner_steps[in_steps_cnt] * sum_coef[in_steps_cnt];
    }
  }
  return answer;
}

boundary_conditions integrator (
    double const                    step_length,
    size_t const                    total_steps,
    boundary_conditions const &     bc_initial,
    double const                    parameter
)
{
  return runge_kutta_method_6 (step_length, total_steps, bc_initial, parameter);
  // return runge_kutta_method_4 (step_length, total_steps, bc_initial, parameter);
}
