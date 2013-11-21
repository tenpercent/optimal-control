#include "runge-kutta_method.hh"
#include "newton_method.hh"
#include "defs.hh"

int newton_method (
    double const                 tau,  // step
    size_t const                 n,    // step numbers
    boundary_conditions &        bc_0, // conditions in 0
    boundary_conditions &        bc_1, // conditions in 1
    double const                 alpha // parametr
)
{
  boundary_conditions start, ret1, ret2;
  size_t i = 0;
  double J11, J12, J21, J22, gamma, next_error, residual, k1, k2, det, next_p1, next_p2;
  
  start = bc_0;
  
  while (i < MAX_ITER)
    {
      // Инициализация итерации
      gamma = 1.0;

      // Найдём элементы матрицы Якоби
      start.p1 = bc_0.p1;
      start.p2 = bc_0.p2;

      ret1 = runge_kutta_method (tau, n, start, alpha);
      
      start.p1 = bc_0.p1 + DELTA;
      start.p2 = bc_0.p2;
      
      ret2 = runge_kutta_method (tau, n, start, alpha);

      J11 = (ret2.x1 - ret1.x1) / DELTA;
      J21 = (ret2.x2 - ret1.x2) / DELTA;

      start.p1 = bc_0.p1;
      start.p2 = bc_0.p2 + DELTA;
      
      ret2 = runge_kutta_method (tau, n, start, alpha);

      J12 = (ret2.x1 - ret1.x1) / DELTA;
      J22 = (ret2.x2 - ret1.x2) / DELTA;

      // Считаем невязку с нормировкой Федоренко
      k1 = J11 * J11 + J12 * J12;
      k2 = J21 * J21 + J22 * J22;

      residual = sqrt ((ret1.x1 - bc_1.x1) * (ret1.x1 - bc_1.x1) / k1 + (ret1.x2 - bc_1.x2) * (ret1.x2 - bc_1.x2) / k2);

      // Если норма достаточно маленькая, заканчиваем алгоритм
      if (residual < EPS)
        {
          bc_1.p1 = ret1.p1;
          bc_1.p2 = ret1.p2;
          return 0;
        }
        
      next_error = 1e+300;

      // Если нет, ищем начальные условия так, чтобы норма уменьшилась относительно предыдущего шага
      while (next_error > residual)
        {
          det = J11 * J22 - J12 * J21;
          
          if (fabs (det) < MINIMAL_FOR_COMPARE)
            return 100; // infernal error

          next_p1 = bc_0.p1 - gamma * ((ret1.x1 - bc_1.x1) * J22 - (ret1.x2 - bc_1.x2) * J12) / det;
          next_p2 = bc_0.p2 - gamma * ((ret1.x2 - bc_1.x2) * J11 - (ret1.x1 - bc_1.x1) * J21) / det;

          gamma *= .5;

          start.p1 = next_p1;
          start.p2 = next_p2;

          ret1 = runge_kutta_method (tau, n, start, alpha);

          start.p1 = next_p1 + DELTA;
          start.p2 = next_p2;
          
          ret2 = runge_kutta_method (tau, n, start, alpha);

          J11 = (ret2.x1 - ret1.x1) / DELTA;
          J21 = (ret2.x2 - ret1.x2) / DELTA;

          start.p1 = next_p1;
          start.p2 = next_p2 + DELTA;

          ret2 = runge_kutta_method (tau, n, start, alpha);

          J12 = (ret2.x1 - ret1.x1) / DELTA;
          J22 = (ret2.x2 - ret1.x2) / DELTA;

          // Считаем невязку с нормировкой Федоренко
          k1 = J11 * J11 + J12 * J12;
          k2 = J21 * J21 + J22 * J22;

          next_error = sqrt ((ret1.x1 - bc_1.x1) * (ret1.x1 - bc_1.x1) / k1 + 
                             (ret1.x2 - bc_1.x2) * (ret1.x2 - bc_1.x2) / k2);

          if (next_error < EPS)
            {
              bc_0.p1 = next_p1;
              bc_0.p2 = next_p2;
              bc_1.p1 = ret1.p1;
              bc_1.p2 = ret1.p2;
              return 0;
            }
        }
      // Заканчиваем итерацию
      bc_0.p1 = next_p1;
      bc_0.p2 = next_p2;
      bc_1.p1 = ret1.p1;
      bc_1.p1 = ret1.p1;

      residual = next_error;
      ++i;
    }  
  return 1;
}
