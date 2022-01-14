/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Hinfinity.h
 *
 * Code generation for function 'Hinfinity'
 *
 */

#ifndef HINFINITY_H
#define HINFINITY_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "Hinfinity_types.h"

/* Function Declarations */
extern void Hinfinity(double g, double duration, double dt, double SteadyState,
                      emxArray_real_T *pos, emxArray_real_T *vel,
                      emxArray_real_T *poshat, emxArray_real_T *velhat,
                      emxArray_real_T *poshatinf, emxArray_real_T *velhatinf,
                      emxArray_real_T *HinfGains, emxArray_real_T *KalmanGains);

#endif

/* End of code generation (Hinfinity.h) */
