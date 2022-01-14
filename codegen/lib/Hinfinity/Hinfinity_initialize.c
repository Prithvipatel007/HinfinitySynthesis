/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Hinfinity_initialize.c
 *
 * Code generation for function 'Hinfinity_initialize'
 *
 */

/* Include files */
#include "Hinfinity_initialize.h"
#include "Hinfinity.h"
#include "Hinfinity_data.h"
#include "eml_rand_mt19937ar_stateful.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void Hinfinity_initialize(void)
{
  rt_InitInfAndNaN();
  c_eml_rand_mt19937ar_stateful_i();
  isInitialized_Hinfinity = true;
}

/* End of code generation (Hinfinity_initialize.c) */
