/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzgghrd.c
 *
 * Code generation for function 'xzgghrd'
 *
 */

/* Include files */
#include "xzgghrd.h"
#include "Hinfinity.h"
#include "Hinfinity_rtwutil.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
void xzgghrd(int ilo, int ihi, creal_T A[4])
{
  int jcol;
  int jrow;
  int f_re_tmp;
  double g2s;
  double f_im;
  double scale_tmp;
  double f2s;
  int y_tmp;
  double scale;
  double fs_re;
  double fs_im;
  double gs_re;
  double gs_im;
  int count;
  int rescaledir;
  boolean_T guard1 = false;
  double c;
  double g2;
  double b_gs_re;
  if (ihi >= ilo + 2) {
    for (jcol = ilo - 1; jcol + 1 < ihi - 1; jcol++) {
      for (jrow = ihi - 2; jrow + 2 > jcol + 2; jrow--) {
        f_re_tmp = jrow + (jcol << 1);
        g2s = A[f_re_tmp].re;
        f_im = A[f_re_tmp].im;
        scale_tmp = fabs(A[f_re_tmp].re);
        f2s = fabs(A[f_re_tmp].im);
        if (f2s > scale_tmp) {
          scale_tmp = f2s;
        }

        y_tmp = f_re_tmp + 1;
        f2s = fabs(A[y_tmp].re);
        scale = fabs(A[y_tmp].im);
        if (scale > f2s) {
          f2s = scale;
        }

        scale = scale_tmp;
        if (f2s > scale_tmp) {
          scale = f2s;
        }

        fs_re = A[f_re_tmp].re;
        fs_im = A[f_re_tmp].im;
        gs_re = A[y_tmp].re;
        gs_im = A[y_tmp].im;
        count = -1;
        rescaledir = 0;
        guard1 = false;
        if (scale >= 7.4428285367870146E+137) {
          do {
            count++;
            fs_re *= 1.3435752215134178E-138;
            fs_im *= 1.3435752215134178E-138;
            gs_re *= 1.3435752215134178E-138;
            gs_im *= 1.3435752215134178E-138;
            scale *= 1.3435752215134178E-138;
          } while (!(scale < 7.4428285367870146E+137));

          rescaledir = 1;
          guard1 = true;
        } else if (scale <= 1.3435752215134178E-138) {
          if ((A[y_tmp].re == 0.0) && (A[y_tmp].im == 0.0)) {
            c = 1.0;
            gs_re = 0.0;
            gs_im = 0.0;
          } else {
            do {
              count++;
              fs_re *= 7.4428285367870146E+137;
              fs_im *= 7.4428285367870146E+137;
              gs_re *= 7.4428285367870146E+137;
              gs_im *= 7.4428285367870146E+137;
              scale *= 7.4428285367870146E+137;
            } while (!(scale > 1.3435752215134178E-138));

            rescaledir = -1;
            guard1 = true;
          }
        } else {
          guard1 = true;
        }

        if (guard1) {
          scale = fs_re * fs_re + fs_im * fs_im;
          g2 = gs_re * gs_re + gs_im * gs_im;
          f2s = g2;
          if (1.0 > g2) {
            f2s = 1.0;
          }

          if (scale <= f2s * 2.0041683600089728E-292) {
            if ((A[f_re_tmp].re == 0.0) && (A[f_re_tmp].im == 0.0)) {
              c = 0.0;
              g2s = rt_hypotd_snf(A[y_tmp].re, A[y_tmp].im);
              f_im = 0.0;
              g2 = rt_hypotd_snf(gs_re, gs_im);
              gs_re /= g2;
              gs_im = -gs_im / g2;
            } else {
              g2s = sqrt(g2);
              c = rt_hypotd_snf(fs_re, fs_im) / g2s;
              if (scale_tmp > 1.0) {
                g2 = rt_hypotd_snf(A[f_re_tmp].re, A[f_re_tmp].im);
                fs_re = A[f_re_tmp].re / g2;
                fs_im = A[f_re_tmp].im / g2;
              } else {
                scale = 7.4428285367870146E+137 * A[f_re_tmp].re;
                f2s = 7.4428285367870146E+137 * A[f_re_tmp].im;
                g2 = rt_hypotd_snf(scale, f2s);
                fs_re = scale / g2;
                fs_im = f2s / g2;
              }

              b_gs_re = gs_re / g2s;
              gs_im = -gs_im / g2s;
              gs_re = fs_re * b_gs_re - fs_im * gs_im;
              gs_im = fs_re * gs_im + fs_im * b_gs_re;
              g2s = c * A[f_re_tmp].re + (gs_re * A[y_tmp].re - gs_im * A[y_tmp]
                .im);
              f_im = c * A[f_re_tmp].im + (gs_re * A[y_tmp].im + gs_im * A[y_tmp]
                .re);
            }
          } else {
            f2s = sqrt(g2 / scale + 1.0);
            g2s = f2s * fs_re;
            f_im = f2s * fs_im;
            c = 1.0 / f2s;
            g2 += scale;
            scale = g2s / g2;
            f2s = f_im / g2;
            b_gs_re = gs_re;
            gs_im = -gs_im;
            gs_re = scale * b_gs_re - f2s * gs_im;
            gs_im = scale * gs_im + f2s * b_gs_re;
            if (rescaledir > 0) {
              for (rescaledir = 0; rescaledir <= count; rescaledir++) {
                g2s *= 7.4428285367870146E+137;
                f_im *= 7.4428285367870146E+137;
              }
            } else {
              if (rescaledir < 0) {
                for (rescaledir = 0; rescaledir <= count; rescaledir++) {
                  g2s *= 1.3435752215134178E-138;
                  f_im *= 1.3435752215134178E-138;
                }
              }
            }
          }
        }

        A[f_re_tmp].re = g2s;
        A[f_re_tmp].im = f_im;
        A[y_tmp].re = 0.0;
        A[y_tmp].im = 0.0;
        b_gs_re = gs_re * A[1].re - gs_im * A[1].im;
        g2 = gs_re * A[1].im + gs_im * A[1].re;
        scale = c * A[0].re;
        f2s = c * A[0].im;
        g2s = A[0].im;
        scale_tmp = A[0].re;
        A[1].re = c * A[1].re - (gs_re * A[0].re + gs_im * A[0].im);
        A[1].im = c * A[1].im - (gs_re * g2s - gs_im * scale_tmp);
        A[0].re = scale + b_gs_re;
        A[0].im = f2s + g2;
        fs_re = c * A[2].re + (gs_re * A[3].re - gs_im * A[3].im);
        fs_im = c * A[2].im + (gs_re * A[3].im + gs_im * A[3].re);
        g2s = A[2].im;
        scale_tmp = A[2].re;
        A[3].re = c * A[3].re - (gs_re * A[2].re + gs_im * A[2].im);
        A[3].im = c * A[3].im - (gs_re * g2s - gs_im * scale_tmp);
        A[2].re = fs_re;
        A[2].im = fs_im;
        gs_re = -gs_re;
        gs_im = -gs_im;
        fs_re = c * A[2].re + (gs_re * A[0].re - gs_im * A[0].im);
        fs_im = c * A[2].im + (gs_re * A[0].im + gs_im * A[0].re);
        g2s = A[2].im;
        scale_tmp = A[2].re;
        A[0].re = c * A[0].re - (gs_re * A[2].re + gs_im * A[2].im);
        A[0].im = c * A[0].im - (gs_re * g2s - gs_im * scale_tmp);
        A[2].re = fs_re;
        A[2].im = fs_im;
        b_gs_re = gs_re * A[1].re - gs_im * A[1].im;
        g2 = gs_re * A[1].im + gs_im * A[1].re;
        scale = c * A[3].re;
        f2s = c * A[3].im;
        g2s = A[3].im;
        scale_tmp = A[3].re;
        A[1].re = c * A[1].re - (gs_re * A[3].re + gs_im * A[3].im);
        A[1].im = c * A[1].im - (gs_re * g2s - gs_im * scale_tmp);
        A[3].re = scale + b_gs_re;
        A[3].im = f2s + g2;
      }
    }
  }
}

/* End of code generation (xzgghrd.c) */
