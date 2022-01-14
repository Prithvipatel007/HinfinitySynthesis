/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_Hinfinity_api.h
 *
 * Code generation for function '_coder_Hinfinity_api'
 *
 */

#ifndef _CODER_HINFINITY_API_H
#define _CODER_HINFINITY_API_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"

/* Type Definitions */
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  real_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_real_T*/

#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T

typedef struct emxArray_real_T emxArray_real_T;

#endif                                 /*typedef_emxArray_real_T*/

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void Hinfinity(real_T g, real_T duration, real_T dt, real_T SteadyState,
                      emxArray_real_T *pos, emxArray_real_T *vel,
                      emxArray_real_T *poshat, emxArray_real_T *velhat,
                      emxArray_real_T *poshatinf, emxArray_real_T *velhatinf,
                      emxArray_real_T *HinfGains, emxArray_real_T *KalmanGains);
extern void Hinfinity_api(const mxArray * const prhs[4], int32_T nlhs, const
  mxArray *plhs[8]);
extern void Hinfinity_atexit(void);
extern void Hinfinity_initialize(void);
extern void Hinfinity_terminate(void);
extern void Hinfinity_xil_shutdown(void);
extern void Hinfinity_xil_terminate(void);

#endif

/* End of code generation (_coder_Hinfinity_api.h) */
