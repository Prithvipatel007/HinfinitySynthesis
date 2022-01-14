/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Hinfinity.c
 *
 * Code generation for function 'Hinfinity'
 *
 */

/* Include files */
#include "Hinfinity.h"
#include "Hinfinity_data.h"
#include "Hinfinity_emxutil.h"
#include "Hinfinity_initialize.h"
#include "Hinfinity_rtwutil.h"
#include "rand.h"
#include "randn.h"
#include "rt_nonfinite.h"
#include "schur.h"
#include "xzgeev.h"
#include <math.h>

/* Function Declarations */
static double rt_powd_snf(double u0, double u1);

/* Function Definitions */
static double rt_powd_snf(double u0, double u1)
{
  double y;
  double d;
  double d1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d = fabs(u0);
    d1 = fabs(u1);
    if (rtIsInf(u1)) {
      if (d == 1.0) {
        y = 1.0;
      } else if (d > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

void Hinfinity(double g, double duration, double dt, double SteadyState,
               emxArray_real_T *pos, emxArray_real_T *vel, emxArray_real_T
               *poshat, emxArray_real_T *velhat, emxArray_real_T *poshatinf,
               emxArray_real_T *velhatinf, emxArray_real_T *HinfGains,
               emxArray_real_T *KalmanGains)
{
  double a[4];
  double r;
  double b_idx_0;
  double x_idx_0;
  double xhat_idx_0;
  double x_idx_1;
  double xhat_idx_1;
  double y;
  double Sw_idx_0;
  double Sw_idx_2_tmp;
  double Sw_idx_3;
  double P[4];
  double xhatinf_idx_0;
  double xhatinf_idx_1;
  double Pinf[4];
  int i;
  int t;
  emxArray_real_T *b_HinfGains;
  boolean_T exitg1;
  boolean_T guard1 = false;
  double T[4];
  double K[2];
  double d;
  double d1;
  double d2;
  double b_t;
  double d3;
  double L_idx_0;
  int info;
  double L_idx_1;
  double L_idx_2;
  double b_T[4];
  double b_a[4];
  static const double b[4] = { 0.01, 0.0, 0.0, 0.01 };

  int b_i;
  static const signed char c_a[4] = { 0, 0, 0, 100 };

  int loop_ub;
  static const double W[4] = { 3.0E-7, 5.0E-6, 5.0E-6, 0.0001 };

  boolean_T p;
  double b_d[2];
  creal_T lambda[2];
  boolean_T exitg3;
  creal_T beta1[2];
  int exitg2;
  if (isInitialized_Hinfinity == false) {
    Hinfinity_initialize();
  }

  /*  function Hinfinity(g, duration, dt, SteadyState) */
  /*  */
  /*  H-infinity filter simulation for a vehicle travelling along a road. */
  /*  This code also simulates the Kalman filter. */
  /*  INPUTS */
  /*    g = gamma */
  /*    if g is too large the program will terminate with an error message. */
  /*    In my simulation I set g = 0.01. */
  /*    duration = length of simulation (seconds). I used duration = 60. */
  /*    dt = step size (seconds). I used dt = 0.1. */
  /*    SteadyState = flag indicating use of steady state filter. */
  /*                  1 = steady state, 0 = time-varying */
  /*  nominal velocity measurement noise (feet/sec) */
  /*  nominal acceleration noise (feet/sec^2) */
  a[0] = 1.0;
  a[2] = dt;

  /*  transition matrix */
  r = dt * dt;
  b_idx_0 = r / 2.0;

  /*  input matrix */
  /*  measurement matrix */
  /*  initial state vector */
  a[1] = 0.0;
  x_idx_0 = 0.0;
  xhat_idx_0 = 0.0;
  a[3] = 1.0;
  x_idx_1 = 0.0;
  xhat_idx_1 = 0.0;
  y = 0.0;

  /*  initial measurement */
  /*  Initialize Kalman filter variables */
  /*  initial Kalman filter state estimate */
  /*  measurement error covariance */
  Sw_idx_0 = 0.040000000000000008 * (rt_powd_snf(dt, 4.0) / 4.0);
  Sw_idx_2_tmp = 0.040000000000000008 * (rt_powd_snf(dt, 3.0) / 2.0);
  Sw_idx_3 = 0.040000000000000008 * r;

  /*  process noise cov */
  P[0] = Sw_idx_0;
  P[1] = Sw_idx_2_tmp;
  P[2] = Sw_idx_2_tmp;
  P[3] = Sw_idx_3;

  /*  initial Kalman filter estimation covariance */
  /*  Initialize H-infinity filter variables */
  xhatinf_idx_0 = 0.0;
  xhatinf_idx_1 = 0.0;

  /*  initial H-infinity filter state estimate */
  Pinf[0] = 0.01;
  Pinf[1] = 0.0;
  Pinf[2] = 0.0;
  Pinf[3] = 0.01;

  /*  Initialize arrays for later plotting. */
  i = pos->size[0];
  pos->size[0] = 1;
  emxEnsureCapacity_real_T(pos, i);
  pos->data[0] = 0.0;

  /*  true position array */
  i = vel->size[0];
  vel->size[0] = 1;
  emxEnsureCapacity_real_T(vel, i);
  vel->data[0] = 0.0;

  /*  true velocity array */
  i = poshat->size[0];
  poshat->size[0] = 1;
  emxEnsureCapacity_real_T(poshat, i);
  poshat->data[0] = 0.0;

  /*  estimated position array (Kalman filter) */
  i = velhat->size[0];
  velhat->size[0] = 1;
  emxEnsureCapacity_real_T(velhat, i);
  velhat->data[0] = 0.0;

  /*  estimated velocity array (Kalman filter) */
  i = poshatinf->size[0];
  poshatinf->size[0] = 1;
  emxEnsureCapacity_real_T(poshatinf, i);
  poshatinf->data[0] = 0.0;

  /*  estimated position array (H-infinity) */
  i = velhatinf->size[0];
  velhatinf->size[0] = 1;
  emxEnsureCapacity_real_T(velhatinf, i);
  velhatinf->data[0] = 0.0;

  /*  estimated velocity array (H-infinity) */
  i = HinfGains->size[0] * HinfGains->size[1];
  HinfGains->size[0] = 2;
  HinfGains->size[1] = 1;
  emxEnsureCapacity_real_T(HinfGains, i);

  /*  H-infinity filter gains */
  i = KalmanGains->size[0] * KalmanGains->size[1];
  KalmanGains->size[0] = 2;
  KalmanGains->size[1] = 1;
  emxEnsureCapacity_real_T(KalmanGains, i);
  HinfGains->data[0] = 0.0;
  KalmanGains->data[0] = 0.0;
  HinfGains->data[1] = 0.0;
  KalmanGains->data[1] = 0.0;

  /*  Kalman filter gains */
  i = (int)(((duration - dt) + dt) / dt);
  t = 0;
  emxInit_real_T(&b_HinfGains, 2);
  exitg1 = false;
  while ((!exitg1) && (t <= i - 1)) {
    /*  Use a constant commanded acceleration of 1 foot/sec^2. */
    /*  Figure out the H-infinity estimate. */
    guard1 = false;
    if (SteadyState == 1.0) {
      /*  Use steady-state H-infinity gains */
      K[0] = 0.11;
      K[1] = 0.09;
      guard1 = true;
    } else {
      T[1] = 0.0;
      T[2] = 0.0;
      T[0] = 1.0;
      T[3] = 1.0;
      d = Pinf[0];
      d1 = Pinf[1];
      d2 = Pinf[2];
      d3 = Pinf[3];
      for (info = 0; info < 2; info++) {
        r = g * b[info];
        b_t = g * b[info + 2];
        b_T[info] = T[info] - (r * d + b_t * d1);
        b_i = c_a[info + 2];
        b_a[info] = 0.0 * d + (double)b_i * d1;
        b_T[info + 2] = T[info + 2] - (r * d2 + b_t * d3);
        b_a[info + 2] = 0.0 * d2 + (double)b_i * d3;
      }

      b_T[0] += b_a[0];
      b_T[1] += b_a[1];
      b_T[2] += b_a[2];
      b_T[3] += b_a[3];
      if (fabs(b_T[1]) > fabs(b_T[0])) {
        r = b_T[0] / b_T[1];
        b_t = 1.0 / (r * b_T[3] - b_T[2]);
        L_idx_0 = b_T[3] / b_T[1] * b_t;
        L_idx_1 = -b_t;
        L_idx_2 = -b_T[2] / b_T[1] * b_t;
        b_t *= r;
      } else {
        r = b_T[1] / b_T[0];
        b_t = 1.0 / (b_T[3] - r * b_T[2]);
        L_idx_0 = b_T[3] / b_T[0] * b_t;
        L_idx_1 = -r * b_t;
        L_idx_2 = -b_T[2] / b_T[0] * b_t;
      }

      for (info = 0; info < 2; info++) {
        d = a[info + 2];
        d1 = (double)(int)a[info] * Pinf[0] + d * Pinf[1];
        d = (double)(int)a[info] * Pinf[2] + d * Pinf[3];
        d2 = d1 * L_idx_0 + d * L_idx_1;
        T[info] = d2;
        d3 = d2 * 0.0;
        d2 = d1 * L_idx_2 + d * b_t;
        T[info + 2] = d2;
        d3 += d2;
        d3 *= 100.0;
        K[info] = d3;
      }

      for (info = 0; info < 2; info++) {
        d = T[info + 2];
        Pinf[info] = (T[info] + d * dt) + W[info];
        Pinf[info + 2] = (T[info] * 0.0 + d) + W[info + 2];
      }

      /*  Force Pinf to be symmetric. */
      b_T[0] = (Pinf[0] + Pinf[0]) / 2.0;
      b_t = (Pinf[1] + Pinf[2]) / 2.0;
      b_T[3] = (Pinf[3] + Pinf[3]) / 2.0;
      Pinf[0] = b_T[0];
      Pinf[1] = b_t;
      Pinf[2] = b_t;
      Pinf[3] = b_T[3];

      /*  Make sure the eigenvalues of Pinf are less than 1 in magnitude. */
      p = ((!rtIsInf(b_T[0])) && (!rtIsNaN(b_T[0])));
      if ((!p) || (rtIsInf(b_t) || rtIsNaN(b_t))) {
        p = false;
      }

      if ((!p) || (rtIsInf(b_t) || rtIsNaN(b_t))) {
        p = false;
      }

      if ((!p) || (rtIsInf(b_T[3]) || rtIsNaN(b_T[3]))) {
        p = false;
      }

      if (!p) {
        lambda[0].re = rtNaN;
        lambda[0].im = 0.0;
        lambda[1].re = rtNaN;
        lambda[1].im = 0.0;
      } else {
        p = true;
        info = 0;
        exitg3 = false;
        while ((!exitg3) && (info < 2)) {
          b_i = 0;
          do {
            exitg2 = 0;
            if (b_i <= info) {
              if (!(Pinf[b_i + (info << 1)] == Pinf[info + (b_i << 1)])) {
                p = false;
                exitg2 = 1;
              } else {
                b_i++;
              }
            } else {
              info++;
              exitg2 = 2;
            }
          } while (exitg2 == 0);

          if (exitg2 == 1) {
            exitg3 = true;
          }
        }

        if (p) {
          T[0] = b_T[0];
          T[1] = b_t;
          T[2] = b_t;
          T[3] = b_T[3];
          schur(T);
          lambda[0].re = T[0];
          lambda[0].im = 0.0;
          lambda[1].re = T[3];
          lambda[1].im = 0.0;
        } else {
          xzgeev(Pinf, &info, lambda, beta1);
          if (beta1[0].im == 0.0) {
            if (lambda[0].im == 0.0) {
              L_idx_1 = lambda[0].re / beta1[0].re;
              r = 0.0;
            } else if (lambda[0].re == 0.0) {
              L_idx_1 = 0.0;
              r = lambda[0].im / beta1[0].re;
            } else {
              L_idx_1 = lambda[0].re / beta1[0].re;
              r = lambda[0].im / beta1[0].re;
            }
          } else if (beta1[0].re == 0.0) {
            if (lambda[0].re == 0.0) {
              L_idx_1 = lambda[0].im / beta1[0].im;
              r = 0.0;
            } else if (lambda[0].im == 0.0) {
              L_idx_1 = 0.0;
              r = -(lambda[0].re / beta1[0].im);
            } else {
              L_idx_1 = lambda[0].im / beta1[0].im;
              r = -(lambda[0].re / beta1[0].im);
            }
          } else {
            L_idx_0 = fabs(beta1[0].re);
            r = fabs(beta1[0].im);
            if (L_idx_0 > r) {
              r = beta1[0].im / beta1[0].re;
              b_t = beta1[0].re + r * beta1[0].im;
              L_idx_1 = (lambda[0].re + r * lambda[0].im) / b_t;
              r = (lambda[0].im - r * lambda[0].re) / b_t;
            } else if (r == L_idx_0) {
              if (beta1[0].re > 0.0) {
                r = 0.5;
              } else {
                r = -0.5;
              }

              if (beta1[0].im > 0.0) {
                b_t = 0.5;
              } else {
                b_t = -0.5;
              }

              L_idx_1 = (lambda[0].re * r + lambda[0].im * b_t) / L_idx_0;
              r = (lambda[0].im * r - lambda[0].re * b_t) / L_idx_0;
            } else {
              r = beta1[0].re / beta1[0].im;
              b_t = beta1[0].im + r * beta1[0].re;
              L_idx_1 = (r * lambda[0].re + lambda[0].im) / b_t;
              r = (r * lambda[0].im - lambda[0].re) / b_t;
            }
          }

          lambda[0].re = L_idx_1;
          lambda[0].im = r;
          if (beta1[1].im == 0.0) {
            if (lambda[1].im == 0.0) {
              L_idx_1 = lambda[1].re / beta1[1].re;
              r = 0.0;
            } else if (lambda[1].re == 0.0) {
              L_idx_1 = 0.0;
              r = lambda[1].im / beta1[1].re;
            } else {
              L_idx_1 = lambda[1].re / beta1[1].re;
              r = lambda[1].im / beta1[1].re;
            }
          } else if (beta1[1].re == 0.0) {
            if (lambda[1].re == 0.0) {
              L_idx_1 = lambda[1].im / beta1[1].im;
              r = 0.0;
            } else if (lambda[1].im == 0.0) {
              L_idx_1 = 0.0;
              r = -(lambda[1].re / beta1[1].im);
            } else {
              L_idx_1 = lambda[1].im / beta1[1].im;
              r = -(lambda[1].re / beta1[1].im);
            }
          } else {
            L_idx_0 = fabs(beta1[1].re);
            r = fabs(beta1[1].im);
            if (L_idx_0 > r) {
              r = beta1[1].im / beta1[1].re;
              b_t = beta1[1].re + r * beta1[1].im;
              L_idx_1 = (lambda[1].re + r * lambda[1].im) / b_t;
              r = (lambda[1].im - r * lambda[1].re) / b_t;
            } else if (r == L_idx_0) {
              if (beta1[1].re > 0.0) {
                r = 0.5;
              } else {
                r = -0.5;
              }

              if (beta1[1].im > 0.0) {
                b_t = 0.5;
              } else {
                b_t = -0.5;
              }

              L_idx_1 = (lambda[1].re * r + lambda[1].im * b_t) / L_idx_0;
              r = (lambda[1].im * r - lambda[1].re * b_t) / L_idx_0;
            } else {
              r = beta1[1].re / beta1[1].im;
              b_t = beta1[1].im + r * beta1[1].re;
              L_idx_1 = (r * lambda[1].re + lambda[1].im) / b_t;
              r = (r * lambda[1].im - lambda[1].re) / b_t;
            }
          }

          lambda[1].re = L_idx_1;
          lambda[1].im = r;
        }
      }

      if ((rt_hypotd_snf(lambda[0].re, lambda[0].im) >= 1.0) || (rt_hypotd_snf
           (lambda[1].re, lambda[1].im) >= 1.0)) {
        exitg1 = true;
      } else {
        guard1 = true;
      }
    }

    if (guard1) {
      b_t = 0.0 * xhatinf_idx_0 + xhatinf_idx_1;
      L_idx_0 = y - b_t;
      L_idx_1 = ((xhatinf_idx_0 + dt * xhatinf_idx_1) + b_idx_0) + K[0] *
        L_idx_0;
      L_idx_2 = (b_t + dt) + K[1] * L_idx_0;
      xhatinf_idx_0 = L_idx_1;
      xhatinf_idx_1 = L_idx_2;
      info = b_HinfGains->size[0] * b_HinfGains->size[1];
      b_HinfGains->size[0] = 2;
      b_HinfGains->size[1] = HinfGains->size[1] + 1;
      emxEnsureCapacity_real_T(b_HinfGains, info);
      loop_ub = HinfGains->size[1];
      for (info = 0; info < loop_ub; info++) {
        b_HinfGains->data[2 * info] = HinfGains->data[2 * info];
        b_i = 2 * info + 1;
        b_HinfGains->data[b_i] = HinfGains->data[b_i];
      }

      b_HinfGains->data[2 * HinfGains->size[1]] = K[0];
      b_HinfGains->data[2 * HinfGains->size[1] + 1] = K[1];
      info = HinfGains->size[0] * HinfGains->size[1];
      HinfGains->size[0] = 2;
      HinfGains->size[1] = b_HinfGains->size[1];
      emxEnsureCapacity_real_T(HinfGains, info);
      loop_ub = b_HinfGains->size[0] * b_HinfGains->size[1];
      for (info = 0; info < loop_ub; info++) {
        HinfGains->data[info] = b_HinfGains->data[info];
      }

      info = poshatinf->size[0];
      b_i = poshatinf->size[0];
      poshatinf->size[0]++;
      emxEnsureCapacity_real_T(poshatinf, b_i);
      poshatinf->data[info] = L_idx_1;
      info = velhatinf->size[0];
      b_i = velhatinf->size[0];
      velhatinf->size[0]++;
      emxEnsureCapacity_real_T(velhatinf, b_i);
      velhatinf->data[info] = L_idx_2;

      /*  Simulate the linear system and noisy measurement. */
      /*  Note that randn is Matlab's Gaussian (normal) random number */
      /*  generator; rand is Matlab's uniform random number generator. */
      b_d[0] = 0.4 * b_idx_0;
      b_d[1] = 0.4 * dt;
      d = randn();
      d1 = randn();
      L_idx_1 = ((x_idx_0 + dt * x_idx_1) + b_idx_0) + b_d[0] * d;
      L_idx_2 = ((0.0 * x_idx_0 + x_idx_1) + dt) + b_d[1] * d1;
      x_idx_0 = L_idx_1;
      x_idx_1 = L_idx_2;
      r = 2.0 * (b_rand() - 0.5);
      y = (0.0 * L_idx_1 + L_idx_2) + r;

      /*  Compute the Kalman filter estimate. */
      /*  Extrapolate the most recent state estimate to the present time. */
      L_idx_1 = (xhat_idx_0 + dt * xhat_idx_1) + b_idx_0;
      L_idx_2 = (0.0 * xhat_idx_0 + xhat_idx_1) + dt;
      info = poshat->size[0];
      b_i = poshat->size[0];
      poshat->size[0]++;
      emxEnsureCapacity_real_T(poshat, b_i);
      poshat->data[info] = L_idx_1;
      info = velhat->size[0];
      b_i = velhat->size[0];
      velhat->size[0]++;
      emxEnsureCapacity_real_T(velhat, b_i);
      velhat->data[info] = L_idx_2;

      /*  Form the Innovation vector. */
      L_idx_0 = y - (0.0 * L_idx_1 + L_idx_2);
      if (SteadyState == 1.0) {
        K[0] = 0.1;
        K[1] = 0.01;
      } else {
        /*  Compute the covariance of the Innovation. */
        /*  Form the Kalman Gain matrix. */
        d = 0.0;
        for (info = 0; info < 2; info++) {
          d1 = a[info + 2];
          d2 = (double)(int)a[info] * P[0] + d1 * P[1];
          T[info] = d2;
          d3 = d2 * 0.0;
          b_i = info << 1;
          d2 = (double)(int)a[info] * P[2] + d1 * P[3];
          T[info + 2] = d2;
          d3 += d2;
          b_d[info] = d3;
          d += (0.0 * P[b_i] + P[b_i + 1]) * (double)info;
        }

        b_t = 1.0 / (d + 4.0);

        /*  Compute the covariance of the estimation error. */
        r = b_d[0] * b_t;
        K[0] = r;
        b_T[0] = r * 0.0;
        b_t *= b_d[1];
        b_T[1] = b_t * 0.0;
        K[1] = b_t;
        b_T[2] = r;
        b_T[3] = b_t;
        for (info = 0; info < 2; info++) {
          d = b_T[info + 2];
          b_a[info] = b_T[info] * P[0] + d * P[1];
          b_a[info + 2] = b_T[info] * P[2] + d * P[3];
        }

        for (info = 0; info < 2; info++) {
          d = T[info + 2];
          d1 = b_a[info + 2];
          b_T[info] = b_a[info] + d1 * dt;
          P[info] = T[info] + d * dt;
          b_T[info + 2] = b_a[info] * 0.0 + d1;
          P[info + 2] = T[info] * 0.0 + d;
        }

        P[0] = (P[0] - b_T[0]) + Sw_idx_0;
        P[1] = (P[1] - b_T[1]) + Sw_idx_2_tmp;
        P[2] = (P[2] - b_T[2]) + Sw_idx_2_tmp;
        P[3] = (P[3] - b_T[3]) + Sw_idx_3;

        /*  Force P to be symmetric. */
        b_T[0] = (P[0] + P[0]) / 2.0;
        b_t = (P[1] + P[2]) / 2.0;
        b_T[3] = (P[3] + P[3]) / 2.0;
        P[0] = b_T[0];
        P[1] = b_t;
        P[2] = b_t;
        P[3] = b_T[3];
      }

      /*  Update the Kalman filter state estimate. */
      xhat_idx_0 = L_idx_1 + K[0] * L_idx_0;
      xhat_idx_1 = L_idx_2 + K[1] * L_idx_0;

      /*  Save some parameters for plotting later. */
      info = b_HinfGains->size[0] * b_HinfGains->size[1];
      b_HinfGains->size[0] = 2;
      b_HinfGains->size[1] = KalmanGains->size[1] + 1;
      emxEnsureCapacity_real_T(b_HinfGains, info);
      loop_ub = KalmanGains->size[1];
      for (info = 0; info < loop_ub; info++) {
        b_HinfGains->data[2 * info] = KalmanGains->data[2 * info];
        b_i = 2 * info + 1;
        b_HinfGains->data[b_i] = KalmanGains->data[b_i];
      }

      b_HinfGains->data[2 * KalmanGains->size[1]] = K[0];
      b_HinfGains->data[2 * KalmanGains->size[1] + 1] = K[1];
      info = KalmanGains->size[0] * KalmanGains->size[1];
      KalmanGains->size[0] = 2;
      KalmanGains->size[1] = b_HinfGains->size[1];
      emxEnsureCapacity_real_T(KalmanGains, info);
      loop_ub = b_HinfGains->size[0] * b_HinfGains->size[1];
      for (info = 0; info < loop_ub; info++) {
        KalmanGains->data[info] = b_HinfGains->data[info];
      }

      info = pos->size[0];
      b_i = pos->size[0];
      pos->size[0]++;
      emxEnsureCapacity_real_T(pos, b_i);
      pos->data[info] = x_idx_0;
      info = vel->size[0];
      b_i = vel->size[0];
      vel->size[0]++;
      emxEnsureCapacity_real_T(vel, b_i);
      vel->data[info] = x_idx_1;
      t++;
    }
  }

  emxFree_real_T(&b_HinfGains);
}

/* End of code generation (Hinfinity.c) */
