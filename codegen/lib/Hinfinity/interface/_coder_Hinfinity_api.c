/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_Hinfinity_api.c
 *
 * Code generation for function '_coder_Hinfinity_api'
 *
 */

/* Include files */
#include "_coder_Hinfinity_api.h"
#include "_coder_Hinfinity_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131483U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "Hinfinity",                         /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

/* Function Declarations */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static const mxArray *b_emlrt_marshallOut(const emxArray_real_T *u);
static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);
static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *g, const
  char_T *identifier);
static const mxArray *emlrt_marshallOut(const emxArray_real_T *u);
static void emxFree_real_T(emxArray_real_T **pEmxArray);
static void emxInit_real_T(const emlrtStack *sp, emxArray_real_T **pEmxArray,
  int32_T numDimensions, boolean_T doPush);

/* Function Definitions */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = c_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static const mxArray *b_emlrt_marshallOut(const emxArray_real_T *u)
{
  const mxArray *y;
  const mxArray *m;
  static const int32_T iv[2] = { 0, 0 };

  y = NULL;
  m = emlrtCreateNumericArray(2, iv, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, &u->data[0]);
  emlrtSetDimensions((mxArray *)m, u->size, 2);
  emlrtAssign(&y, m);
  return y;
}

static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *g, const
  char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(g), &thisId);
  emlrtDestroyArray(&g);
  return y;
}

static const mxArray *emlrt_marshallOut(const emxArray_real_T *u)
{
  const mxArray *y;
  const mxArray *m;
  static const int32_T iv[1] = { 0 };

  y = NULL;
  m = emlrtCreateNumericArray(1, iv, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, &u->data[0]);
  emlrtSetDimensions((mxArray *)m, u->size, 1);
  emlrtAssign(&y, m);
  return y;
}

static void emxFree_real_T(emxArray_real_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_real_T *)NULL) {
    if (((*pEmxArray)->data != (real_T *)NULL) && (*pEmxArray)->canFreeData) {
      emlrtFreeMex((*pEmxArray)->data);
    }

    emlrtFreeMex((*pEmxArray)->size);
    emlrtFreeMex(*pEmxArray);
    *pEmxArray = (emxArray_real_T *)NULL;
  }
}

static void emxInit_real_T(const emlrtStack *sp, emxArray_real_T **pEmxArray,
  int32_T numDimensions, boolean_T doPush)
{
  emxArray_real_T *emxArray;
  int32_T i;
  *pEmxArray = (emxArray_real_T *)emlrtMallocMex(sizeof(emxArray_real_T));
  if (doPush) {
    emlrtPushHeapReferenceStackR2012b(sp, (void *)pEmxArray, (void (*)(void *))
      emxFree_real_T);
  }

  emxArray = *pEmxArray;
  emxArray->data = (real_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)emlrtMallocMex(sizeof(int32_T) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

void Hinfinity_api(const mxArray * const prhs[4], int32_T nlhs, const mxArray
                   *plhs[8])
{
  emxArray_real_T *pos;
  emxArray_real_T *vel;
  emxArray_real_T *poshat;
  emxArray_real_T *velhat;
  emxArray_real_T *poshatinf;
  emxArray_real_T *velhatinf;
  emxArray_real_T *HinfGains;
  emxArray_real_T *KalmanGains;
  real_T g;
  real_T duration;
  real_T dt;
  real_T SteadyState;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  emxInit_real_T(&st, &pos, 1, true);
  emxInit_real_T(&st, &vel, 1, true);
  emxInit_real_T(&st, &poshat, 1, true);
  emxInit_real_T(&st, &velhat, 1, true);
  emxInit_real_T(&st, &poshatinf, 1, true);
  emxInit_real_T(&st, &velhatinf, 1, true);
  emxInit_real_T(&st, &HinfGains, 2, true);
  emxInit_real_T(&st, &KalmanGains, 2, true);

  /* Marshall function inputs */
  g = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "g");
  duration = emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "duration");
  dt = emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "dt");
  SteadyState = emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "SteadyState");

  /* Invoke the target function */
  Hinfinity(g, duration, dt, SteadyState, pos, vel, poshat, velhat, poshatinf,
            velhatinf, HinfGains, KalmanGains);

  /* Marshall function outputs */
  pos->canFreeData = false;
  plhs[0] = emlrt_marshallOut(pos);
  emxFree_real_T(&pos);
  if (nlhs > 1) {
    vel->canFreeData = false;
    plhs[1] = emlrt_marshallOut(vel);
  }

  emxFree_real_T(&vel);
  if (nlhs > 2) {
    poshat->canFreeData = false;
    plhs[2] = emlrt_marshallOut(poshat);
  }

  emxFree_real_T(&poshat);
  if (nlhs > 3) {
    velhat->canFreeData = false;
    plhs[3] = emlrt_marshallOut(velhat);
  }

  emxFree_real_T(&velhat);
  if (nlhs > 4) {
    poshatinf->canFreeData = false;
    plhs[4] = emlrt_marshallOut(poshatinf);
  }

  emxFree_real_T(&poshatinf);
  if (nlhs > 5) {
    velhatinf->canFreeData = false;
    plhs[5] = emlrt_marshallOut(velhatinf);
  }

  emxFree_real_T(&velhatinf);
  if (nlhs > 6) {
    HinfGains->canFreeData = false;
    plhs[6] = b_emlrt_marshallOut(HinfGains);
  }

  emxFree_real_T(&HinfGains);
  if (nlhs > 7) {
    KalmanGains->canFreeData = false;
    plhs[7] = b_emlrt_marshallOut(KalmanGains);
  }

  emxFree_real_T(&KalmanGains);
  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

void Hinfinity_atexit(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  Hinfinity_xil_terminate();
  Hinfinity_xil_shutdown();
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void Hinfinity_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

void Hinfinity_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (_coder_Hinfinity_api.c) */
