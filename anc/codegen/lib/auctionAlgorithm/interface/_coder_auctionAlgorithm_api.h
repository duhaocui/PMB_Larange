/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_auctionAlgorithm_api.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 02-Jan-2018 11:22:03
 */

#ifndef _CODER_AUCTIONALGORITHM_API_H
#define _CODER_AUCTIONALGORITHM_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_auctionAlgorithm_api.h"

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
extern void auctionAlgorithm(emxArray_real_T *RewardMatrix, emxArray_real_T
  *Customer2Item, emxArray_real_T *Item2Customer);
extern void auctionAlgorithm_api(const mxArray *prhs[1], const mxArray *plhs[2]);
extern void auctionAlgorithm_atexit(void);
extern void auctionAlgorithm_initialize(void);
extern void auctionAlgorithm_terminate(void);
extern void auctionAlgorithm_xil_terminate(void);

#endif

/*
 * File trailer for _coder_auctionAlgorithm_api.h
 *
 * [EOF]
 */
