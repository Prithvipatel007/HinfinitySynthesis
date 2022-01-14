/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Hinfinity_rtwutil.c
 *
 * Code generation for function 'Hinfinity_rtwutil'
 *
 */

/* Include files */
#include "Hinfinity_rtwutil.h"
#include "Hinfinity.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
double rt_hypotd_snf(double u0, double u1)
{
  double y;
  double a;
  a = fabs(u0);
  y = fabs(u1);
  if (a < y) {
    a /= y;
    y *= sqrt(a * a + 1.0);
  } else if (a > y) {
    y /= a;
    y = a * sqrt(y * y + 1.0);
  } else {
    if (!rtIsNaN(y)) {
      y = a * 1.4142135623730951;
    }
  }

  return y;
}

/* End of code generation (Hinfinity_rtwutil.c) */
