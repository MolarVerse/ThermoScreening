#include <stdlib.h>    // macros: EXIT_SUCCESS, EXIT_FAILURE
#include <stdio.h>     // function: printf, scanf
#include <string.h>    // function:
#include <strings.h>   // function:
#include <errno.h>     // function:
#include <math.h>      // function:
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


void eigenValues(double tensor[3][3], double eigVals[3])
{
  // store the tensor elements in different variables for better readability 
  /**
   *        [a b c]
   * T =    [d e f]
   *        [g h i]
   */

  double t00 = tensor[0][0];
  double t01 = tensor[0][1];
  double t02 = tensor[0][2];
  double t10 = tensor[1][0];
  double t11 = tensor[1][1];
  double t12 = tensor[1][2];
  double t20 = tensor[2][0];
  double t21 = tensor[2][1];
  double t22 = tensor[2][2];

  /**
   * Eigenvalue equation det(d*I-T) = 0
   * Transform into cubic equation:
   * A*x^3 + B*x^2 + C*x + D = 0
   */


  // Cubic coefficient: minus sign is pulled out
  double a = -1.0;

  // Quadratic coefficient: Sum of diagonal elements
  double b = t00 + t11 + t22;

  // Linear coefficient: 
  double c = -t11 * t22 - t00 * t22 - t00 * t11 + t20 * t02 + t01 * t10 + t21 * t12;

  // Constant Coefficient: 
  double d = t00 * t11 * t22 + t02 * t10 * t21 + t01 * t12 * t20 - t20 * t02 * t11 - t10 * t01 * t22 - t21 * t12 * t00;

   /**
   * Solve Cubic equations: see therefore
   * https://stackoverflow.com/questions/13328676/c-solving-cubic-equations
   *
  */

  /**
   * Cubic discrimant: D = B^2*C^2-4*A*C^3-4*B^3*D+18*ABCD-27*A^2*D^2
   * D = q^3+r^2 with q = (3*C - B^2)/9, r = (-(27*D) + B * (9*C - 2*B^2))) / 54
   * 
   */

  // Transform to normal form x^3 + B*x^2 + C*x^ + D = 0
  // double a = 1.0;
  // double b = -2.0;
  // double c = -5.0;
  // double d = 6.0;
  b /= a;
  c /= a;
  d /= a;

  /** 
   * The cubic equation can be transformed into a depressed form:
   * z^3 + p*z + q = 0, with p = (3*C - B^2 / 9 and q = (-(27*D) + B*(9*C - 2*B^2)) / 54)
   * The discriminant is then D = p^3 + q^2
   */

  double p = (3.0 * c - (b * b)) / 9.0;
  double q = (-(27.0 * d) + b * (9.0 * c - 2.0 * (b * b))) / 54.0;
  double disc = p * p * p  + q * q;

  double term1 = b / 3.0;

  if (disc > 0)
  {
    double s = q + sqrt(disc);
    s = (s < 0) ? -pow(-s, (1.0 / 3.0)) : pow(s, 1.0 / 3.0);
    double t = q - sqrt(disc);
    t = (t < 0) ? -pow(-t, (1.0 / 3.0)) : pow(t, 1.0 / 3.0);
    eigVals[0] = -term1 + s + t;
    term1 += (s + t) / 2.0;
    eigVals[1] = -term1;
    eigVals[2] = -term1;
    return;
  }

  double q13;

  if (disc == 0)
  {
    q13 = (q < 0) ? -pow(-q, (1.0 / 3.0)) : pow(q, (1.0 / 3.0));
    eigVals[0] = -term1 + 2.0 * q13;
    eigVals[1] = -(q13 + term1);
    eigVals[2] = -(q13 + term1);
    return;
  }

  p *= -1.0;
  double dum1 = p * p * p;
  dum1 = acos(q / sqrt(dum1));
  q13 = 2.0 * sqrt(p);
  eigVals[0] = -term1 + q13 * cos(dum1 / 3.0);
  eigVals[1] = -term1 + q13 * cos((dum1 + 2.0 * M_PI) / 3.0);
  eigVals[2] = -term1 + q13 * cos((dum1 + 4.0 * M_PI) / 3.0);

  printf("Eigenvalues in Angstrom: %.12lf %.12lf %.12lf\n\n", eigVals[0], eigVals[1], eigVals[2]);

  return;
}
int main(void){
  double matrix[3][3] = {{11.0,-11.0,-7.0},{7.0,-7.0,-5.0},{3.0,-3.0,-1.0}};
  double eigVal[3]={0};
  eigenValues(matrix, eigVal);
  for (int i = 0; i<3; i++){
      printf("%f, ",eigVal[i]);
  }
  return EXIT_SUCCESS;
}