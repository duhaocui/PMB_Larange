/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_AuctionJacobi_api.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 02-Jan-2018 11:57:42
 */

#ifndef _CODER_AUCTIONJACOBI_API_H
#define _CODER_AUCTIONJACOBI_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_AuctionJacobi_api.h"

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
extern void AuctionJacobi(emxArray_real_T *c, emxArray_real_T *assignment,
  emxArray_real_T *prices);
extern void AuctionJacobi_api(const mxArray * const prhs[1], const mxArray *
  plhs[2]);
extern void AuctionJacobi_atexit(void);
extern void AuctionJacobi_initialize(void);
extern void AuctionJacobi_terminate(void);
extern void AuctionJacobi_xil_terminate(void);

#endif

/*
 * File trailer for _coder_AuctionJacobi_api.h
 *
 * [EOF]
 */
