//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: AuctionJacobi.cpp
//
// MATLAB Coder version            : 3.4
// C/C++ source code generated on  : 02-Jan-2018 11:57:42
//

// Include Files
#include "rt_nonfinite.h"
#include "AuctionJacobi.h"
#include "AuctionJacobi_emxutil.h"

// Function Declarations
static void PerformRoundAuction(emxArray_real_T *assignment, const
  emxArray_real_T *prices, const emxArray_real_T *c, double epsilon,
  emxArray_real_T *v);

// Function Definitions

//
// Arguments    : emxArray_real_T *assignment
//                const emxArray_real_T *prices
//                const emxArray_real_T *c
//                double epsilon
//                emxArray_real_T *v
// Return Type  : void
//
static void PerformRoundAuction(emxArray_real_T *assignment, const
  emxArray_real_T *prices, const emxArray_real_T *c, double epsilon,
  emxArray_real_T *v)
{
  int i1;
  int nx;
  emxArray_boolean_T *b;
  emxArray_int32_T *ii;
  int idx;
  int nxin;
  boolean_T exitg1;
  emxArray_int32_T *unAssignedPeople;
  emxArray_real_T *temp;
  int i;
  emxArray_real_T *b_index;
  emxArray_int32_T *indices;
  emxArray_int32_T *r0;
  emxArray_int32_T *b_indices;
  int n;
  double mtmp;
  int itmp;
  double b_mtmp;
  i1 = v->size[0] * v->size[1];
  v->size[0] = 1;
  v->size[1] = prices->size[1];
  emxEnsureCapacity_real_T(v, i1);
  nx = prices->size[0] * prices->size[1];
  for (i1 = 0; i1 < nx; i1++) {
    v->data[i1] = prices->data[i1];
  }

  emxInit_boolean_T(&b, 2);
  i1 = b->size[0] * b->size[1];
  b->size[0] = 1;
  b->size[1] = assignment->size[1];
  emxEnsureCapacity_boolean_T(b, i1);
  nx = assignment->size[0] * assignment->size[1];
  for (i1 = 0; i1 < nx; i1++) {
    b->data[i1] = rtIsInf(assignment->data[i1]);
  }

  emxInit_int32_T(&ii, 2);
  nx = b->size[1];
  idx = 0;
  i1 = ii->size[0] * ii->size[1];
  ii->size[0] = 1;
  ii->size[1] = b->size[1];
  emxEnsureCapacity_int32_T(ii, i1);
  nxin = 1;
  exitg1 = false;
  while ((!exitg1) && (nxin <= nx)) {
    if (b->data[nxin - 1]) {
      idx++;
      ii->data[idx - 1] = nxin;
      if (idx >= nx) {
        exitg1 = true;
      } else {
        nxin++;
      }
    } else {
      nxin++;
    }
  }

  if (b->size[1] == 1) {
    if (idx == 0) {
      i1 = ii->size[0] * ii->size[1];
      ii->size[0] = 1;
      ii->size[1] = 0;
      emxEnsureCapacity_int32_T(ii, i1);
    }
  } else {
    i1 = ii->size[0] * ii->size[1];
    if (1 > idx) {
      ii->size[1] = 0;
    } else {
      ii->size[1] = idx;
    }

    emxEnsureCapacity_int32_T(ii, i1);
  }

  emxInit_int32_T(&unAssignedPeople, 2);
  i1 = unAssignedPeople->size[0] * unAssignedPeople->size[1];
  unAssignedPeople->size[0] = 1;
  unAssignedPeople->size[1] = ii->size[1];
  emxEnsureCapacity_int32_T(unAssignedPeople, i1);
  nx = ii->size[0] * ii->size[1];
  for (i1 = 0; i1 < nx; i1++) {
    unAssignedPeople->data[i1] = ii->data[i1];
  }

  emxInit_real_T(&temp, 2);
  i1 = temp->size[0] * temp->size[1];
  temp->size[0] = 2;
  temp->size[1] = unAssignedPeople->size[1];
  emxEnsureCapacity_real_T(temp, i1);
  nx = unAssignedPeople->size[1] << 1;
  for (i1 = 0; i1 < nx; i1++) {
    temp->data[i1] = 0.0;
  }

  // compute and store the bids of each unsassigned individual in temp
  i = 0;
  emxInit_real_T(&b_index, 2);
  while (i <= unAssignedPeople->size[1] - 1) {
    nx = c->size[1];
    nxin = unAssignedPeople->data[i];
    i1 = b_index->size[0] * b_index->size[1];
    b_index->size[0] = 1;
    b_index->size[1] = nx;
    emxEnsureCapacity_real_T(b_index, i1);
    for (i1 = 0; i1 < nx; i1++) {
      b_index->data[b_index->size[0] * i1] = c->data[(nxin + c->size[0] * i1) -
        1] - prices->data[prices->size[0] * i1];
    }

    idx = 1;
    n = b_index->size[1];
    mtmp = b_index->data[0];
    itmp = 1;
    if (b_index->size[1] > 1) {
      if (rtIsNaN(b_index->data[0])) {
        nxin = 2;
        exitg1 = false;
        while ((!exitg1) && (nxin <= n)) {
          idx = nxin;
          if (!rtIsNaN(b_index->data[nxin - 1])) {
            mtmp = b_index->data[nxin - 1];
            itmp = nxin;
            exitg1 = true;
          } else {
            nxin++;
          }
        }
      }

      if (idx < b_index->size[1]) {
        while (idx + 1 <= n) {
          if (b_index->data[idx] > mtmp) {
            mtmp = b_index->data[idx];
            itmp = idx + 1;
          }

          idx++;
        }
      }
    }

    nxin = b_index->size[1] - 1;
    for (nx = itmp; nx <= nxin; nx++) {
      b_index->data[nx - 1] = b_index->data[nx];
    }

    if (1 > nxin) {
      i1 = 0;
    } else {
      i1 = nxin;
    }

    nxin = b_index->size[0] * b_index->size[1];
    b_index->size[1] = i1;
    emxEnsureCapacity_real_T(b_index, nxin);
    idx = 1;
    b_mtmp = b_index->data[0];
    if (i1 > 1) {
      if (rtIsNaN(b_index->data[0])) {
        nxin = 2;
        exitg1 = false;
        while ((!exitg1) && (nxin <= i1)) {
          idx = nxin;
          if (!rtIsNaN(b_index->data[nxin - 1])) {
            b_mtmp = b_index->data[nxin - 1];
            exitg1 = true;
          } else {
            nxin++;
          }
        }
      }

      if (idx < i1) {
        while (idx + 1 <= i1) {
          if (b_index->data[idx] > b_mtmp) {
            b_mtmp = b_index->data[idx];
          }

          idx++;
        }
      }
    }

    temp->data[temp->size[0] * i] = itmp;
    temp->data[1 + temp->size[0] * i] = (mtmp - b_mtmp) + epsilon;
    i++;
  }

  // each object which has received a bid determines the highest bidder and
  // update its price accordingly
  i = 0;
  emxInit_int32_T(&indices, 2);
  emxInit_int32_T1(&r0, 1);
  emxInit_int32_T1(&b_indices, 1);
  while (i <= prices->size[1] - 1) {
    nx = temp->size[1];
    i1 = b->size[0] * b->size[1];
    b->size[0] = 1;
    b->size[1] = nx;
    emxEnsureCapacity_boolean_T(b, i1);
    for (i1 = 0; i1 < nx; i1++) {
      b->data[b->size[0] * i1] = (temp->data[temp->size[0] * i1] == 1.0 +
        (double)i);
    }

    nx = b->size[1];
    idx = 0;
    i1 = ii->size[0] * ii->size[1];
    ii->size[0] = 1;
    ii->size[1] = b->size[1];
    emxEnsureCapacity_int32_T(ii, i1);
    nxin = 1;
    exitg1 = false;
    while ((!exitg1) && (nxin <= nx)) {
      if (b->data[nxin - 1]) {
        idx++;
        ii->data[idx - 1] = nxin;
        if (idx >= nx) {
          exitg1 = true;
        } else {
          nxin++;
        }
      } else {
        nxin++;
      }
    }

    if (b->size[1] == 1) {
      if (idx == 0) {
        i1 = ii->size[0] * ii->size[1];
        ii->size[0] = 1;
        ii->size[1] = 0;
        emxEnsureCapacity_int32_T(ii, i1);
      }
    } else {
      i1 = ii->size[0] * ii->size[1];
      if (1 > idx) {
        ii->size[1] = 0;
      } else {
        ii->size[1] = idx;
      }

      emxEnsureCapacity_int32_T(ii, i1);
    }

    i1 = indices->size[0] * indices->size[1];
    indices->size[0] = 1;
    indices->size[1] = ii->size[1];
    emxEnsureCapacity_int32_T(indices, i1);
    nx = ii->size[0] * ii->size[1];
    for (i1 = 0; i1 < nx; i1++) {
      indices->data[i1] = ii->data[i1];
    }

    if (!(indices->size[1] == 0)) {
      i1 = b_index->size[0] * b_index->size[1];
      b_index->size[0] = 1;
      b_index->size[1] = indices->size[1];
      emxEnsureCapacity_real_T(b_index, i1);
      nx = indices->size[1];
      for (i1 = 0; i1 < nx; i1++) {
        b_index->data[b_index->size[0] * i1] = temp->data[1 + temp->size[0] *
          (indices->data[indices->size[0] * i1] - 1)];
      }

      idx = 1;
      i1 = b_indices->size[0];
      b_indices->size[0] = indices->size[1];
      emxEnsureCapacity_int32_T1(b_indices, i1);
      nx = indices->size[1];
      for (i1 = 0; i1 < nx; i1++) {
        b_indices->data[i1] = indices->data[indices->size[0] * i1];
      }

      n = b_indices->size[0];
      i1 = r0->size[0];
      r0->size[0] = indices->size[1];
      emxEnsureCapacity_int32_T1(r0, i1);
      nx = indices->size[1];
      for (i1 = 0; i1 < nx; i1++) {
        r0->data[i1] = indices->data[indices->size[0] * i1];
      }

      mtmp = temp->data[1 + temp->size[0] * (r0->data[0] - 1)];
      itmp = 0;
      i1 = b_indices->size[0];
      b_indices->size[0] = indices->size[1];
      emxEnsureCapacity_int32_T1(b_indices, i1);
      nx = indices->size[1];
      for (i1 = 0; i1 < nx; i1++) {
        b_indices->data[i1] = indices->data[indices->size[0] * i1];
      }

      if (b_indices->size[0] > 1) {
        if (rtIsNaN(mtmp)) {
          nxin = 2;
          exitg1 = false;
          while ((!exitg1) && (nxin <= n)) {
            idx = nxin;
            if (!rtIsNaN(b_index->data[nxin - 1])) {
              mtmp = b_index->data[nxin - 1];
              itmp = nxin - 1;
              exitg1 = true;
            } else {
              nxin++;
            }
          }
        }

        i1 = b_indices->size[0];
        b_indices->size[0] = indices->size[1];
        emxEnsureCapacity_int32_T1(b_indices, i1);
        nx = indices->size[1];
        for (i1 = 0; i1 < nx; i1++) {
          b_indices->data[i1] = indices->data[indices->size[0] * i1];
        }

        if (idx < b_indices->size[0]) {
          while (idx + 1 <= n) {
            if (b_index->data[idx] > mtmp) {
              mtmp = b_index->data[idx];
              itmp = idx;
            }

            idx++;
          }
        }
      }

      i1 = b->size[0] * b->size[1];
      b->size[0] = 1;
      b->size[1] = assignment->size[1];
      emxEnsureCapacity_boolean_T(b, i1);
      nx = assignment->size[0] * assignment->size[1];
      for (i1 = 0; i1 < nx; i1++) {
        b->data[i1] = (assignment->data[i1] == 1.0 + (double)i);
      }

      nx = b->size[1];
      idx = 0;
      i1 = ii->size[0] * ii->size[1];
      ii->size[0] = 1;
      ii->size[1] = b->size[1];
      emxEnsureCapacity_int32_T(ii, i1);
      nxin = 1;
      exitg1 = false;
      while ((!exitg1) && (nxin <= nx)) {
        if (b->data[nxin - 1]) {
          idx++;
          ii->data[idx - 1] = nxin;
          if (idx >= nx) {
            exitg1 = true;
          } else {
            nxin++;
          }
        } else {
          nxin++;
        }
      }

      if (b->size[1] == 1) {
        if (idx == 0) {
          i1 = ii->size[0] * ii->size[1];
          ii->size[0] = 1;
          ii->size[1] = 0;
          emxEnsureCapacity_int32_T(ii, i1);
        }
      } else {
        i1 = ii->size[0] * ii->size[1];
        if (1 > idx) {
          ii->size[1] = 0;
        } else {
          ii->size[1] = idx;
        }

        emxEnsureCapacity_int32_T(ii, i1);
      }

      i1 = b_index->size[0] * b_index->size[1];
      b_index->size[0] = 1;
      b_index->size[1] = ii->size[1];
      emxEnsureCapacity_real_T(b_index, i1);
      nx = ii->size[0] * ii->size[1];
      for (i1 = 0; i1 < nx; i1++) {
        b_index->data[i1] = ii->data[i1];
      }

      if (!(b_index->size[1] == 0)) {
        assignment->data[(int)b_index->data[0] - 1] = rtInf;
      }

      assignment->data[unAssignedPeople->data[indices->data[itmp] - 1] - 1] =
        1.0 + (double)i;
      v->data[i] += mtmp;
    }

    i++;
  }

  emxFree_int32_T(&b_indices);
  emxFree_int32_T(&r0);
  emxFree_int32_T(&ii);
  emxFree_boolean_T(&b);
  emxFree_real_T(&b_index);
  emxFree_int32_T(&indices);
  emxFree_real_T(&temp);
  emxFree_int32_T(&unAssignedPeople);
}

//
// AUCTIONJACOBI Compute optimal assignement and optimal prices by Bertsekas
//  algorithm. c is a matrix ; assignement is a vector giving the assigned
//  column of a given row and the vector prices stores the prices of each
//  column.
// Arguments    : const emxArray_real_T *c
//                emxArray_real_T *assignment
//                emxArray_real_T *prices
// Return Type  : void
//
void AuctionJacobi(const emxArray_real_T *c, emxArray_real_T *assignment,
                   emxArray_real_T *prices)
{
  int varargin_1;
  int loop_ub;
  int b_loop_ub;
  int i0;
  double epsilon;
  emxArray_boolean_T *b;
  emxArray_real_T *b_prices;
  int exitg1;
  double x;

  //  This function returns altogether the optimal assignement and the dual
  //  prices relative to the cost function c, i.e it solves
  //  max c(1,sigma(1)) + ... + c(N,sigma(N)) where sigma is a permutation of
  //  {1,..,N} and c is a N-by-N matrix. The dual prices are the solution of
  //  min sum_j(v_j) + sum_i(Max_j(c_ij - v(i)))
  //  This is one of the numerous implementation of Pr. Dimitri Bertsekas' auction 
  //  algorithm. Reference papers can be found on his page
  //  http://web.mit.edu/dimitrib/www/home.html
  //  Implemented by Damien Bosc (Ecole Polytechnique, France), last modified
  //  9/7/09
  varargin_1 = c->size[0];
  loop_ub = c->size[0];
  b_loop_ub = c->size[0];
  i0 = assignment->size[0] * assignment->size[1];
  assignment->size[0] = 1;
  assignment->size[1] = b_loop_ub;
  emxEnsureCapacity_real_T(assignment, i0);
  for (i0 = 0; i0 < b_loop_ub; i0++) {
    assignment->data[i0] = rtInf;
  }

  i0 = prices->size[0] * prices->size[1];
  prices->size[0] = 1;
  prices->size[1] = varargin_1;
  emxEnsureCapacity_real_T(prices, i0);
  for (i0 = 0; i0 < varargin_1; i0++) {
    prices->data[i0] = 1.0;
  }

  epsilon = 1.0;
  emxInit_boolean_T(&b, 2);
  emxInit_real_T(&b_prices, 2);
  while (epsilon > 1.0 / (double)loop_ub) {
    i0 = assignment->size[0] * assignment->size[1];
    assignment->size[0] = 1;
    assignment->size[1] = loop_ub;
    emxEnsureCapacity_real_T(assignment, i0);
    for (i0 = 0; i0 < loop_ub; i0++) {
      assignment->data[i0] = rtInf;
    }

    do {
      exitg1 = 0;
      i0 = b->size[0] * b->size[1];
      b->size[0] = 1;
      b->size[1] = assignment->size[1];
      emxEnsureCapacity_boolean_T(b, i0);
      b_loop_ub = assignment->size[0] * assignment->size[1];
      for (i0 = 0; i0 < b_loop_ub; i0++) {
        b->data[i0] = rtIsInf(assignment->data[i0]);
      }

      if (b->size[1] == 0) {
        x = 0.0;
      } else {
        x = b->data[0];
        for (varargin_1 = 2; varargin_1 <= b->size[1]; varargin_1++) {
          x += (double)b->data[varargin_1 - 1];
        }
      }

      if (x != 0.0) {
        i0 = b_prices->size[0] * b_prices->size[1];
        b_prices->size[0] = 1;
        b_prices->size[1] = prices->size[1];
        emxEnsureCapacity_real_T(b_prices, i0);
        b_loop_ub = prices->size[0] * prices->size[1];
        for (i0 = 0; i0 < b_loop_ub; i0++) {
          b_prices->data[i0] = prices->data[i0];
        }

        PerformRoundAuction(assignment, b_prices, c, epsilon, prices);
      } else {
        exitg1 = 1;
      }
    } while (exitg1 == 0);

    // epsilon scaling as recommended by Bertsekas
    epsilon *= 0.25;
  }

  emxFree_real_T(&b_prices);
  emxFree_boolean_T(&b);
}

//
// File trailer for AuctionJacobi.cpp
//
// [EOF]
//
