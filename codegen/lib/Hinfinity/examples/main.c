/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * main.c
 *
 * Code generation for function 'main'
 *
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/

/* Include files */
#include "main.h"
#include "Hinfinity.h"
#include "Hinfinity_emxAPI.h"
#include "Hinfinity_terminate.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static double argInit_real_T(void);
static void main_Hinfinity(void);

/* Function Definitions */
static double argInit_real_T(void)
{
  return 0.0;
}

static void main_Hinfinity(void)
{
  emxArray_real_T *pos;
  emxArray_real_T *vel;
  emxArray_real_T *poshat;
  emxArray_real_T *velhat;
  emxArray_real_T *poshatinf;
  emxArray_real_T *velhatinf;
  emxArray_real_T *HinfGains;
  emxArray_real_T *KalmanGains;
  double g_tmp_tmp;
  emxInitArray_real_T(&pos, 1);
  emxInitArray_real_T(&vel, 1);
  emxInitArray_real_T(&poshat, 1);
  emxInitArray_real_T(&velhat, 1);
  emxInitArray_real_T(&poshatinf, 1);
  emxInitArray_real_T(&velhatinf, 1);
  emxInitArray_real_T(&HinfGains, 2);
  emxInitArray_real_T(&KalmanGains, 2);

  /* Initialize function 'Hinfinity' input arguments. */
  g_tmp_tmp = argInit_real_T();

  /* Call the entry-point 'Hinfinity'. */
  Hinfinity(g_tmp_tmp, g_tmp_tmp, g_tmp_tmp, argInit_real_T(), pos, vel, poshat,
            velhat, poshatinf, velhatinf, HinfGains, KalmanGains);
  emxDestroyArray_real_T(KalmanGains);
  emxDestroyArray_real_T(HinfGains);
  emxDestroyArray_real_T(velhatinf);
  emxDestroyArray_real_T(poshatinf);
  emxDestroyArray_real_T(velhat);
  emxDestroyArray_real_T(poshat);
  emxDestroyArray_real_T(vel);
  emxDestroyArray_real_T(pos);
}

int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* The initialize function is being called automatically from your entry-point function. So, a call to initialize is not included here. */
  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_Hinfinity();

  /* Terminate the application.
     You do not need to do this more than one time. */
  Hinfinity_terminate();
  return 0;
}

/* End of code generation (main.c) */
