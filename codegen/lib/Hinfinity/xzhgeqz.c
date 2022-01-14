/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzhgeqz.c
 *
 * Code generation for function 'xzhgeqz'
 *
 */

/* Include files */
#include "xzhgeqz.h"
#include "Hinfinity.h"
#include "Hinfinity_rtwutil.h"
#include "rt_nonfinite.h"
#include "xzlartg.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
void xzhgeqz(const creal_T A[4], int ilo, int ihi, int *info, creal_T alpha1[2],
             creal_T beta1[2])
{
  creal_T b_A[4];
  double eshift_re;
  double eshift_im;
  creal_T ctemp;
  double anorm;
  double scale;
  double b_atol;
  boolean_T firstNonZero;
  int j;
  int jp1;
  double ascale;
  int i;
  int re_tmp;
  boolean_T guard1 = false;
  double temp1;
  boolean_T guard2 = false;
  int ifirst;
  double temp2;
  int istart;
  int ilast;
  int ilastm1;
  int ilastm;
  int iiter;
  boolean_T goto60;
  boolean_T goto70;
  boolean_T goto90;
  int jiter;
  int exitg1;
  boolean_T b_guard1 = false;
  boolean_T guard3 = false;
  boolean_T exitg2;
  creal_T shift;
  double ascale_im;
  double ad22_re;
  double ad22_im;
  double t1_re;
  double t1_im;
  double t1_im_tmp;
  creal_T b_ascale;
  double shift_re;
  double shift_im;
  memcpy(&b_A[0], &A[0], 4U * sizeof(creal_T));
  *info = 0;
  alpha1[0].re = 0.0;
  alpha1[0].im = 0.0;
  beta1[0].re = 1.0;
  beta1[0].im = 0.0;
  alpha1[1].re = 0.0;
  alpha1[1].im = 0.0;
  beta1[1].re = 1.0;
  beta1[1].im = 0.0;
  eshift_re = 0.0;
  eshift_im = 0.0;
  ctemp.re = 0.0;
  ctemp.im = 0.0;
  anorm = 0.0;
  if (ilo <= ihi) {
    scale = 0.0;
    anorm = 0.0;
    firstNonZero = true;
    for (j = ilo; j <= ihi; j++) {
      jp1 = j + 1;
      if (ihi < j + 1) {
        jp1 = ihi;
      }

      for (i = ilo; i <= jp1; i++) {
        re_tmp = (i + ((j - 1) << 1)) - 1;
        if (A[re_tmp].re != 0.0) {
          temp1 = fabs(A[re_tmp].re);
          if (firstNonZero) {
            anorm = 1.0;
            scale = temp1;
            firstNonZero = false;
          } else if (scale < temp1) {
            temp2 = scale / temp1;
            anorm = anorm * temp2 * temp2 + 1.0;
            scale = temp1;
          } else {
            temp2 = temp1 / scale;
            anorm += temp2 * temp2;
          }
        }

        if (A[re_tmp].im != 0.0) {
          temp1 = fabs(A[re_tmp].im);
          if (firstNonZero) {
            anorm = 1.0;
            scale = temp1;
            firstNonZero = false;
          } else if (scale < temp1) {
            temp2 = scale / temp1;
            anorm = anorm * temp2 * temp2 + 1.0;
            scale = temp1;
          } else {
            temp2 = temp1 / scale;
            anorm += temp2 * temp2;
          }
        }
      }
    }

    anorm = scale * sqrt(anorm);
  }

  scale = 2.2204460492503131E-16 * anorm;
  b_atol = 2.2250738585072014E-308;
  if (scale > 2.2250738585072014E-308) {
    b_atol = scale;
  }

  scale = 2.2250738585072014E-308;
  if (anorm > 2.2250738585072014E-308) {
    scale = anorm;
  }

  ascale = 1.0 / scale;
  firstNonZero = true;
  if (ihi + 1 <= 2) {
    alpha1[1] = A[3];
  }

  guard1 = false;
  guard2 = false;
  if (ihi >= ilo) {
    ifirst = ilo;
    istart = ilo;
    ilast = ihi - 1;
    ilastm1 = ihi - 2;
    ilastm = ihi;
    iiter = 0;
    goto60 = false;
    goto70 = false;
    goto90 = false;
    jiter = 0;
    do {
      exitg1 = 0;
      if (jiter <= 30 * ((ihi - ilo) + 1) - 1) {
        b_guard1 = false;
        if (ilast + 1 == ilo) {
          goto60 = true;
          b_guard1 = true;
        } else {
          jp1 = ilast + (ilastm1 << 1);
          if (fabs(b_A[jp1].re) + fabs(b_A[jp1].im) <= b_atol) {
            b_A[jp1].re = 0.0;
            b_A[jp1].im = 0.0;
            goto60 = true;
            b_guard1 = true;
          } else {
            j = ilastm1;
            guard3 = false;
            exitg2 = false;
            while ((!exitg2) && (j + 1 >= ilo)) {
              if (j + 1 == ilo) {
                guard3 = true;
                exitg2 = true;
              } else if (fabs(b_A[j].re) + fabs(b_A[j].im) <= b_atol) {
                b_A[j].re = 0.0;
                b_A[j].im = 0.0;
                guard3 = true;
                exitg2 = true;
              } else {
                j--;
                guard3 = false;
              }
            }

            if (guard3) {
              ifirst = j + 1;
              goto70 = true;
            }

            if (goto70) {
              b_guard1 = true;
            } else {
              alpha1[0].re = rtNaN;
              alpha1[0].im = 0.0;
              beta1[0].re = rtNaN;
              beta1[0].im = 0.0;
              alpha1[1].re = rtNaN;
              alpha1[1].im = 0.0;
              beta1[1].re = rtNaN;
              beta1[1].im = 0.0;
              *info = 1;
              exitg1 = 1;
            }
          }
        }

        if (b_guard1) {
          if (goto60) {
            goto60 = false;
            alpha1[ilast] = b_A[ilast + (ilast << 1)];
            ilast = ilastm1;
            ilastm1--;
            if (ilast + 1 < ilo) {
              firstNonZero = false;
              guard2 = true;
              exitg1 = 1;
            } else {
              iiter = 0;
              eshift_re = 0.0;
              eshift_im = 0.0;
              ilastm = ilast + 1;
              jiter++;
            }
          } else {
            if (goto70) {
              goto70 = false;
              iiter++;
              if (iiter - iiter / 10 * 10 != 0) {
                jp1 = ilastm1 + (ilastm1 << 1);
                anorm = ascale * b_A[jp1].re;
                scale = ascale * b_A[jp1].im;
                if (scale == 0.0) {
                  shift.re = anorm / 0.70710678118654746;
                  shift.im = 0.0;
                } else if (anorm == 0.0) {
                  shift.re = 0.0;
                  shift.im = scale / 0.70710678118654746;
                } else {
                  shift.re = anorm / 0.70710678118654746;
                  shift.im = scale / 0.70710678118654746;
                }

                jp1 = ilast + (ilast << 1);
                anorm = ascale * b_A[jp1].re;
                scale = ascale * b_A[jp1].im;
                if (scale == 0.0) {
                  ad22_re = anorm / 0.70710678118654746;
                  ad22_im = 0.0;
                } else if (anorm == 0.0) {
                  ad22_re = 0.0;
                  ad22_im = scale / 0.70710678118654746;
                } else {
                  ad22_re = anorm / 0.70710678118654746;
                  ad22_im = scale / 0.70710678118654746;
                }

                t1_re = 0.5 * (shift.re + ad22_re);
                t1_im = 0.5 * (shift.im + ad22_im);
                t1_im_tmp = t1_re * t1_im;
                jp1 = ilastm1 + (ilast << 1);
                anorm = ascale * b_A[jp1].re;
                scale = ascale * b_A[jp1].im;
                if (scale == 0.0) {
                  temp2 = anorm / 0.70710678118654746;
                  ascale_im = 0.0;
                } else if (anorm == 0.0) {
                  temp2 = 0.0;
                  ascale_im = scale / 0.70710678118654746;
                } else {
                  temp2 = anorm / 0.70710678118654746;
                  ascale_im = scale / 0.70710678118654746;
                }

                jp1 = ilast + (ilastm1 << 1);
                anorm = ascale * b_A[jp1].re;
                scale = ascale * b_A[jp1].im;
                if (scale == 0.0) {
                  temp1 = anorm / 0.70710678118654746;
                  anorm = 0.0;
                } else if (anorm == 0.0) {
                  temp1 = 0.0;
                  anorm = scale / 0.70710678118654746;
                } else {
                  temp1 = anorm / 0.70710678118654746;
                  anorm = scale / 0.70710678118654746;
                }

                shift_re = shift.re * ad22_re - shift.im * ad22_im;
                shift_im = shift.re * ad22_im + shift.im * ad22_re;
                shift.re = ((t1_re * t1_re - t1_im * t1_im) + (temp2 * temp1 -
                  ascale_im * anorm)) - shift_re;
                shift.im = ((t1_im_tmp + t1_im_tmp) + (temp2 * anorm + ascale_im
                  * temp1)) - shift_im;
                if (shift.im == 0.0) {
                  if (shift.re < 0.0) {
                    anorm = 0.0;
                    scale = sqrt(-shift.re);
                  } else {
                    anorm = sqrt(shift.re);
                    scale = 0.0;
                  }
                } else if (shift.re == 0.0) {
                  if (shift.im < 0.0) {
                    anorm = sqrt(-shift.im / 2.0);
                    scale = -anorm;
                  } else {
                    anorm = sqrt(shift.im / 2.0);
                    scale = anorm;
                  }
                } else if (rtIsNaN(shift.re)) {
                  anorm = shift.re;
                  scale = shift.re;
                } else if (rtIsNaN(shift.im)) {
                  anorm = shift.im;
                  scale = shift.im;
                } else if (rtIsInf(shift.im)) {
                  anorm = fabs(shift.im);
                  scale = shift.im;
                } else if (rtIsInf(shift.re)) {
                  if (shift.re < 0.0) {
                    anorm = 0.0;
                    scale = shift.im * -shift.re;
                  } else {
                    anorm = shift.re;
                    scale = 0.0;
                  }
                } else {
                  scale = fabs(shift.re);
                  anorm = fabs(shift.im);
                  if ((scale > 4.4942328371557893E+307) || (anorm >
                       4.4942328371557893E+307)) {
                    scale *= 0.5;
                    anorm = rt_hypotd_snf(scale, anorm * 0.5);
                    if (anorm > scale) {
                      anorm = sqrt(anorm) * sqrt(scale / anorm + 1.0);
                    } else {
                      anorm = sqrt(anorm) * 1.4142135623730951;
                    }
                  } else {
                    anorm = sqrt((rt_hypotd_snf(scale, anorm) + scale) * 0.5);
                  }

                  if (shift.re > 0.0) {
                    scale = 0.5 * (shift.im / anorm);
                  } else {
                    if (shift.im < 0.0) {
                      scale = -anorm;
                    } else {
                      scale = anorm;
                    }

                    anorm = 0.5 * (shift.im / scale);
                  }
                }

                if ((t1_re - ad22_re) * anorm + (t1_im - ad22_im) * scale <= 0.0)
                {
                  shift.re = t1_re + anorm;
                  shift.im = t1_im + scale;
                } else {
                  shift.re = t1_re - anorm;
                  shift.im = t1_im - scale;
                }
              } else {
                jp1 = ilast + (ilastm1 << 1);
                anorm = ascale * b_A[jp1].re;
                scale = ascale * b_A[jp1].im;
                if (scale == 0.0) {
                  temp2 = anorm / 0.70710678118654746;
                  ascale_im = 0.0;
                } else if (anorm == 0.0) {
                  temp2 = 0.0;
                  ascale_im = scale / 0.70710678118654746;
                } else {
                  temp2 = anorm / 0.70710678118654746;
                  ascale_im = scale / 0.70710678118654746;
                }

                eshift_re += temp2;
                eshift_im += ascale_im;
                shift.re = eshift_re;
                shift.im = eshift_im;
              }

              j = ilastm1;
              jp1 = ilastm1 + 1;
              exitg2 = false;
              while ((!exitg2) && (j + 1 > ifirst)) {
                istart = 2;
                ctemp.re = ascale * b_A[3].re - shift.re * 0.70710678118654746;
                ctemp.im = ascale * b_A[3].im - shift.im * 0.70710678118654746;
                anorm = fabs(ctemp.re) + fabs(ctemp.im);
                temp2 = ascale * (fabs(b_A[jp1 + 2].re) + fabs(b_A[jp1 + 2].im));
                scale = anorm;
                if (temp2 > anorm) {
                  scale = temp2;
                }

                if ((scale < 1.0) && (scale != 0.0)) {
                  anorm /= scale;
                  temp2 /= scale;
                }

                if ((fabs(b_A[1].re) + fabs(b_A[1].im)) * temp2 <= anorm *
                    b_atol) {
                  goto90 = true;
                  exitg2 = true;
                } else {
                  jp1 = 1;
                  j = 0;
                }
              }

              if (!goto90) {
                istart = ifirst;
                jp1 = (ifirst + ((ifirst - 1) << 1)) - 1;
                ctemp.re = ascale * b_A[jp1].re - shift.re * 0.70710678118654746;
                ctemp.im = ascale * b_A[jp1].im - shift.im * 0.70710678118654746;
              }

              goto90 = false;
              jp1 = ((istart - 1) << 1) + 1;
              b_ascale.re = ascale * b_A[jp1].re;
              b_ascale.im = ascale * b_A[jp1].im;
              xzlartg(ctemp, b_ascale, &temp2, &shift);
              j = istart;
              while (j < ilast + 1) {
                for (j = 1; j <= ilastm; j++) {
                  re_tmp = (j - 1) << 1;
                  anorm = b_A[re_tmp].re;
                  scale = b_A[re_tmp].im;
                  jp1 = re_tmp + 1;
                  shift_re = shift.re * b_A[jp1].re - shift.im * b_A[jp1].im;
                  shift_im = shift.re * b_A[jp1].im + shift.im * b_A[jp1].re;
                  temp1 = temp2 * b_A[jp1].im - (shift.re * b_A[re_tmp].im -
                    shift.im * b_A[re_tmp].re);
                  b_A[jp1].re = temp2 * b_A[jp1].re - (shift.re * b_A[re_tmp].re
                    + shift.im * b_A[re_tmp].im);
                  b_A[jp1].im = temp1;
                  anorm = temp2 * anorm + shift_re;
                  b_A[re_tmp].re = anorm;
                  scale = temp2 * scale + shift_im;
                  b_A[re_tmp].im = scale;
                }

                shift.re = -shift.re;
                shift.im = -shift.im;
                for (i = ifirst; i < 3; i++) {
                  anorm = b_A[i + 1].re;
                  scale = b_A[i - 1].re;
                  temp1 = b_A[i - 1].im;
                  ad22_im = temp2 * b_A[i + 1].im + (shift.re * temp1 + shift.im
                    * scale);
                  b_A[i - 1].re = temp2 * b_A[i - 1].re - (shift.re * b_A[i + 1]
                    .re + shift.im * b_A[i + 1].im);
                  b_A[i - 1].im = temp2 * b_A[i - 1].im - (shift.re * b_A[i + 1]
                    .im - shift.im * anorm);
                  b_A[i + 1].re = temp2 * anorm + (shift.re * scale - shift.im *
                    temp1);
                  b_A[i + 1].im = ad22_im;
                }

                j = 2;
              }
            }

            jiter++;
          }
        }
      } else {
        guard2 = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  } else {
    guard1 = true;
  }

  if (guard2) {
    if (firstNonZero) {
      *info = ilast + 1;
      for (jp1 = 0; jp1 <= ilast; jp1++) {
        alpha1[jp1].re = rtNaN;
        alpha1[jp1].im = 0.0;
        beta1[jp1].re = rtNaN;
        beta1[jp1].im = 0.0;
      }
    } else {
      guard1 = true;
    }
  }

  if (guard1) {
    for (j = 0; j <= ilo - 2; j++) {
      alpha1[j] = b_A[j + (j << 1)];
    }
  }
}

/* End of code generation (xzhgeqz.c) */
