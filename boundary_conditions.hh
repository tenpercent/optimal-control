#ifndef BOUNDARYCONDITIONS 
#define BOUNDARYCONDITIONS 

struct boundary_conditions {
  double x1;
  double x2;
  double p1;
  double p2;
  boundary_conditions (double const t1, double const t2, double const t3, double const t4):
    x1(t1),
    x2(t2),
    p1(t3),
    p2(t4)
  {}
  boundary_conditions ()
    {}
  boundary_conditions & operator+= (boundary_conditions const & second);
  boundary_conditions operator+    (boundary_conditions const & second) const;
  boundary_conditions & operator*= (double const multiplicator);
  boundary_conditions operator*    (double const multiplicator) const;
  boundary_conditions & operator/= (double const divisor);
  boundary_conditions operator/    (double const divisor) const;
};

#endif