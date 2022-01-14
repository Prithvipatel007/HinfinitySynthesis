/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzgeev.c
 *
 * Code generation for function 'xzgeev'
 *
 */

/* Include files */
#include "xzgeev.h"
#include "Hinfinity.h"
#include "Hinfinity_rtwutil.h"
#include "rt_nonfinite.h"
#include "xzgghrd.h"
#include "xzhgeqz.h"

/* Function Definitions */
void xzgeev(const double A[4], int *info, creal_T alpha1[2], creal_T beta1[2])
{
  creal_T At[4];
  double anrm;
  int k;
  boolean_T exitg1;
  double absxk;
  boolean_T ilascl;
  double anrmto;
  boolean_T guard1 = false;
  int ilo;
  double ctoc;
  int ihi;
  boolean_T notdone;
  int j;
  double cfrom1;
  int ii;
  double cto1;
  double a;
  int nzcount;
  int jj;
  boolean_T exitg2;
  int At_tmp;
  At[0].re = A[0];
  At[0].im = 0.0;
  At[1].re = A[1];
  At[1].im = 0.0;
  At[2].re = A[2];
  At[2].im = 0.0;
  At[3].re = A[3];
  At[3].im = 0.0;
  *info = 0;
  anrm = 0.0;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 4)) {
    absxk = rt_hypotd_snf(At[k].re, 0.0);
    if (rtIsNaN(absxk)) {
      anrm = rtNaN;
      exitg1 = true;
    } else {
      if (absxk > anrm) {
        anrm = absxk;
      }

      k++;
    }
  }

  if (rtIsInf(anrm) || rtIsNaN(anrm)) {
    alpha1[0].re = rtNaN;
    alpha1[0].im = 0.0;
    beta1[0].re = rtNaN;
    beta1[0].im = 0.0;
    alpha1[1].re = rtNaN;
    alpha1[1].im = 0.0;
    beta1[1].re = rtNaN;
    beta1[1].im = 0.0;
  } else {
    ilascl = false;
    anrmto = anrm;
    guard1 = false;
    if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
      anrmto = 6.7178761075670888E-139;
      ilascl = true;
      guard1 = true;
    } else {
      if (anrm > 1.4885657073574029E+138) {
        anrmto = 1.4885657073574029E+138;
        ilascl = true;
        guard1 = true;
      }
    }

    if (guard1) {
      absxk = anrm;
      ctoc = anrmto;
      notdone = true;
      while (notdone) {
        cfrom1 = absxk * 2.0041683600089728E-292;
        cto1 = ctoc / 4.9896007738368E+291;
        if ((cfrom1 > ctoc) && (ctoc != 0.0)) {
          a = 2.0041683600089728E-292;
          absxk = cfrom1;
        } else if (cto1 > absxk) {
          a = 4.9896007738368E+291;
          ctoc = cto1;
        } else {
          a = ctoc / absxk;
          notdone = false;
        }

        At[0].re *= a;
        At[0].im *= a;
        At[1].re *= a;
        At[1].im *= a;
        At[2].re *= a;
        At[2].im *= a;
        At[3].re *= a;
        At[3].im *= a;
      }
    }

    ilo = 1;
    ihi = 2;
    k = 0;
    j = 0;
    notdone = false;
    ii = 2;
    exitg1 = false;
    while ((!exitg1) && (ii > 0)) {
      nzcount = 0;
      k = ii;
      j = 2;
      jj = 0;
      exitg2 = false;
      while ((!exitg2) && (jj < 2)) {
        At_tmp = (ii + (jj << 1)) - 1;
        if ((At[At_tmp].re != 0.0) || (At[At_tmp].im != 0.0) || (ii == jj + 1))
        {
          if (nzcount == 0) {
            j = jj + 1;
            nzcount = 1;
            jj++;
          } else {
            nzcount = 2;
            exitg2 = true;
          }
        } else {
          jj++;
        }
      }

      if (nzcount < 2) {
        notdone = true;
        exitg1 = true;
      } else {
        ii--;
      }
    }

    if (!notdone) {
      k = 0;
      j = 0;
      notdone = false;
      jj = 1;
      exitg1 = false;
      while ((!exitg1) && (jj < 3)) {
        nzcount = 0;
        k = 2;
        j = jj;
        ii = 1;
        exitg2 = false;
        while ((!exitg2) && (ii < 3)) {
          At_tmp = (ii + ((jj - 1) << 1)) - 1;
          if ((At[At_tmp].re != 0.0) || (At[At_tmp].im != 0.0) || (ii == jj)) {
            if (nzcount == 0) {
              k = ii;
              nzcount = 1;
              ii++;
            } else {
              nzcount = 2;
              exitg2 = true;
            }
          } else {
            ii++;
          }
        }

        if (nzcount < 2) {
          notdone = true;
          exitg1 = true;
        } else {
          jj++;
        }
      }

      if (notdone) {
        if (k != 1) {
          absxk = At[k - 1].re;
          ctoc = At[k - 1].im;
          At[k - 1] = At[0];
          At[0].re = absxk;
          At[0].im = ctoc;
          absxk = At[k + 1].re;
          ctoc = At[k + 1].im;
          At[k + 1] = At[2];
          At[2].re = absxk;
          At[2].im = ctoc;
        }

        if (j != 1) {
          k = (j - 1) << 1;
          absxk = At[k].re;
          ctoc = At[k].im;
          At[k] = At[0];
          At[0].re = absxk;
          At[0].im = ctoc;
          k++;
          absxk = At[k].re;
          ctoc = At[k].im;
          At[k] = At[1];
          At[1].re = absxk;
          At[1].im = ctoc;
        }

        ilo = 2;
      }
    } else {
      if (k != 2) {
        absxk = At[k - 1].re;
        ctoc = At[k - 1].im;
        At[k - 1] = At[1];
        At[1].re = absxk;
        At[1].im = ctoc;
        absxk = At[k + 1].re;
        ctoc = At[k + 1].im;
        At[k + 1] = At[3];
        At[3].re = absxk;
        At[3].im = ctoc;
      }

      if (j != 2) {
        k = (j - 1) << 1;
        absxk = At[k].re;
        ctoc = At[k].im;
        At[k] = At[2];
        At[2].re = absxk;
        At[2].im = ctoc;
        k++;
        absxk = At[k].re;
        ctoc = At[k].im;
        At[k] = At[3];
        At[3].re = absxk;
        At[3].im = ctoc;
      }

      ihi = 1;
    }

    xzgghrd(ilo, ihi, At);
    xzhgeqz(At, ilo, ihi, info, alpha1, beta1);
    if ((*info == 0) && ilascl) {
      notdone = true;
      while (notdone) {
        cfrom1 = anrmto * 2.0041683600089728E-292;
        cto1 = anrm / 4.9896007738368E+291;
        if ((cfrom1 > anrm) && (anrm != 0.0)) {
          a = 2.0041683600089728E-292;
          anrmto = cfrom1;
        } else if (cto1 > anrmto) {
          a = 4.9896007738368E+291;
          anrm = cto1;
        } else {
          a = anrm / anrmto;
          notdone = false;
        }

        alpha1[0].re *= a;
        alpha1[0].im *= a;
        alpha1[1].re *= a;
        alpha1[1].im *= a;
      }
    }
  }
}

/* End of code generation (xzgeev.c) */
