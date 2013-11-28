#include "boundary_conditions.hh"

boundary_conditions & boundary_conditions::operator+= (boundary_conditions const & second) {
  x1 += second.x1;
  x2 += second.x2;
  p1 += second.p1;
  p2 += second.p2;
  return *this;
}

boundary_conditions boundary_conditions::operator+ (boundary_conditions const & second) const {
  boundary_conditions answer = *this;
  answer += second;
  return answer;
}

boundary_conditions & boundary_conditions::operator*= (double const multiplicator) {
  x1 *= multiplicator;
  x2 *= multiplicator;
  p1 *= multiplicator;
  p2 *= multiplicator;
  return *this;
}

boundary_conditions boundary_conditions::operator* (double const multiplicator) const {
  boundary_conditions answer = *this;
  answer *= multiplicator;
  return answer;
}

boundary_conditions & boundary_conditions::operator/= (double const divisor) {
  *this *= (1. / divisor);
  return *this;
}

boundary_conditions boundary_conditions::operator/ (double const divisor) const {
  boundary_conditions answer = *this;
  answer *= (1. / divisor);
  return answer;
}
