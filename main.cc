#include <cstdio>

#include "newton_method.hh"
#include "calculate_functional.hh"
#include "defs.hh"

int main () {
  double alpha;
  printf ("Enter alpha\n");
  scanf ("%lf", &alpha);
  
  size_t n;
  printf ("Enter number of steps\n");
  scanf ("%ld", &n);
  
  double tau = 1. / n;
  printf ("Step of Runge-Kutta algorithm = %.4lf(%le)\n", tau, tau);
  
  boundary_conditions bc_0 (X1_0, X2_0, P1_0, P2_0); 
  boundary_conditions bc_1 (X1_1, X2_1, 0, 0);
  
  if (newton_method (tau, n, bc_0, bc_1, alpha))
    {
      printf ("Acceleration of Newton algorithm broke down.\n");
      return 1;
    }
  
  printf ("Вoundary conditions:\n");
  printf ("\tx1(0) = %.4lf\tx1(1) = %.4lf\n", bc_0.x1, bc_1.x1);
  printf ("\tx2(0) = %.4lf\tx2(1) = %.4lf\n", bc_0.x2, bc_1.x2);
  printf ("\tp1(0) = %.4lf\tp1(1) = %.4lf\n", bc_0.p1, bc_1.p1);
  printf ("\tp2(0) = %.4lf\tp2(1) = %.4lf\n", bc_0.p2, bc_1.p2);
  
  double rez = calculate_functional (tau, n, bc_0, alpha);
  
  printf ("answer = %.5lf(%le)\n", rez, rez);
  
  return 0;
}
