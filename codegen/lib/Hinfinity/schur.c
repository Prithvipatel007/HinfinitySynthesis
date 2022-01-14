/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * schur.c
 *
 * Code generation for function 'schur'
 *
 */

/* Include files */
#include "schur.h"
#include "Hinfinity.h"
#include "Hinfinity_rtwutil.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
void schur(double A[4])
{
  double htmp2;
  boolean_T p;
  int i;
  int k;
  boolean_T exitg1;
  double ba;
  double bb;
  int b_i;
  int i1;
  double d;
  double ab;
  double d1;
  double a;
  double cs;
  double aa;
  double tau;
  double sn;
  double b_p;
  int b_A;
  int c_A;
  htmp2 = A[0];
  p = ((!rtIsInf(htmp2)) && (!rtIsNaN(htmp2)));
  if (p) {
    htmp2 = A[1];
    if (rtIsInf(htmp2) || rtIsNaN(htmp2)) {
      p = false;
    }
  } else {
    p = false;
  }

  if (p) {
    htmp2 = A[2];
    if (rtIsInf(htmp2) || rtIsNaN(htmp2)) {
      p = false;
    }
  } else {
    p = false;
  }

  if (p) {
    htmp2 = A[3];
    if (rtIsInf(htmp2) || rtIsNaN(htmp2)) {
      p = false;
    }
  } else {
    p = false;
  }

  if (!p) {
    A[0] = rtNaN;
    A[1] = rtNaN;
    A[2] = rtNaN;
    A[3] = rtNaN;
    A[1] = 0.0;
  } else {
    for (i = 1; i + 1 >= 1; i = k - 2) {
      k = i + 1;
      exitg1 = false;
      while ((!exitg1) && (k > 1)) {
        ba = fabs(A[1]);
        if (ba <= 2.0041683600089728E-292) {
          exitg1 = true;
        } else {
          bb = fabs(A[3]);
          if (ba <= 2.2204460492503131E-16 * (fabs(A[0]) + bb)) {
            htmp2 = fabs(A[2]);
            if (ba > htmp2) {
              ab = ba;
              ba = htmp2;
            } else {
              ab = htmp2;
            }

            htmp2 = fabs(A[0] - A[3]);
            if (bb > htmp2) {
              aa = bb;
              bb = htmp2;
            } else {
              aa = htmp2;
            }

            htmp2 = aa + ab;
            aa = 2.2204460492503131E-16 * (bb * (aa / htmp2));
            if ((2.0041683600089728E-292 > aa) || rtIsNaN(aa)) {
              aa = 2.0041683600089728E-292;
            }

            if (ba * (ab / htmp2) <= aa) {
              exitg1 = true;
            } else {
              k = 1;
            }
          } else {
            k = 1;
          }
        }
      }

      if (k > 1) {
        A[1] = 0.0;
      }

      if ((k != i + 1) && (k == i)) {
        b_i = i << 1;
        i1 = i + b_i;
        bb = A[i1];
        d = A[i];
        d1 = A[b_i];
        a = A[0];
        if (A[i] == 0.0) {
          cs = 1.0;
          sn = 0.0;
        } else if (A[b_i] == 0.0) {
          cs = 0.0;
          sn = 1.0;
          bb = A[0];
          a = A[i1];
          d1 = -A[i];
          d = 0.0;
        } else {
          tau = A[0] - A[i1];
          if ((tau == 0.0) && ((A[b_i] < 0.0) != (A[i] < 0.0))) {
            cs = 1.0;
            sn = 0.0;
          } else {
            b_p = 0.5 * tau;
            ab = fabs(A[b_i]);
            aa = fabs(A[i]);
            if ((ab > aa) || rtIsNaN(aa)) {
              bb = ab;
            } else {
              bb = aa;
            }

            aa = fabs(A[i]);
            if ((ab < aa) || rtIsNaN(aa)) {
              aa = ab;
            }

            if (!(A[i << 1] < 0.0)) {
              b_A = 1;
            } else {
              b_A = -1;
            }

            if (!(A[i] < 0.0)) {
              c_A = 1;
            } else {
              c_A = -1;
            }

            ab = aa * (double)b_A * (double)c_A;
            htmp2 = fabs(b_p);
            if ((!(htmp2 > bb)) && (!rtIsNaN(bb))) {
              htmp2 = bb;
            }

            ba = b_p / htmp2 * b_p + bb / htmp2 * ab;
            if (ba >= 8.8817841970012523E-16) {
              a = sqrt(htmp2) * sqrt(ba);
              if (b_p < 0.0) {
                a = -a;
              }

              ba = b_p + a;
              a = A[i1] + ba;
              bb = A[i1] - bb / ba * ab;
              tau = rt_hypotd_snf(A[i], ba);
              cs = ba / tau;
              sn = A[i] / tau;
              d1 = A[b_i] - A[i];
              d = 0.0;
            } else {
              htmp2 = A[b_i] + A[i];
              tau = rt_hypotd_snf(htmp2, tau);
              cs = sqrt(0.5 * (fabs(htmp2) / tau + 1.0));
              if (!(htmp2 < 0.0)) {
                b_A = 1;
              } else {
                b_A = -1;
              }

              sn = -(b_p / (tau * cs)) * (double)b_A;
              aa = A[0] * cs + A[b_i] * sn;
              bb = -A[0] * sn + A[b_i] * cs;
              htmp2 = A[i] * cs + A[i1] * sn;
              ab = -A[i] * sn + A[i1] * cs;
              d1 = bb * cs + ab * sn;
              d = -aa * sn + htmp2 * cs;
              aa = 0.5 * ((aa * cs + htmp2 * sn) + (-bb * sn + ab * cs));
              a = aa;
              bb = aa;
              if (d != 0.0) {
                if (d1 != 0.0) {
                  if ((d1 < 0.0) == (d < 0.0)) {
                    htmp2 = sqrt(fabs(d1));
                    ba = sqrt(fabs(d));
                    a = htmp2 * ba;
                    if (!(d < 0.0)) {
                      b_p = a;
                    } else {
                      b_p = -a;
                    }

                    tau = 1.0 / sqrt(fabs(d1 + d));
                    a = aa + b_p;
                    bb = aa - b_p;
                    d1 -= d;
                    d = 0.0;
                    ab = htmp2 * tau;
                    htmp2 = ba * tau;
                    aa = cs * ab - sn * htmp2;
                    sn = cs * htmp2 + sn * ab;
                    cs = aa;
                  }
                } else {
                  d1 = -d;
                  d = 0.0;
                  aa = cs;
                  cs = -sn;
                  sn = aa;
                }
              }
            }
          }
        }

        A[0] = a;
        A[b_i] = d1;
        A[i] = d;
        A[i1] = bb;
        if (2 > i + 1) {
          ab = cs * A[2];
          htmp2 = sn * A[2];
          A[2] = ab - htmp2;
          A[2] = ab + htmp2;
        }
      }
    }
  }
}

/* End of code generation (schur.c) */
