//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: auctionAlgorithm.cpp
//
// MATLAB Coder version            : 3.4
// C/C++ source code generated on  : 02-Jan-2018 11:22:03
//

// Include Files
#include "rt_nonfinite.h"
#include "auctionAlgorithm.h"
#include "auctionAlgorithm_emxutil.h"

// Function Definitions

//
// Arguments    : emxArray_real_T *RewardMatrix
//                emxArray_real_T *Customer2Item
//                emxArray_real_T *Item2Customer
// Return Type  : void
//
void auctionAlgorithm(emxArray_real_T *RewardMatrix, emxArray_real_T
                      *Customer2Item, emxArray_real_T *Item2Customer)
{
  int NofItems;
  int NofCustomers;
  int idx;
  int k;
  emxArray_boolean_T *x;
  emxArray_real_T *varargin_1;
  emxArray_int32_T *iindx;
  emxArray_real_T *b_RewardMatrix;
  int exitg1;
  int cindx;
  boolean_T exitg2;
  int i;
  unsigned int unnamed_idx_1;
  int n;
  int ixstart;
  double mtmp;
  int itmp;
  double b_mtmp;

  // amount of deviation from the optimal reward
  NofItems = RewardMatrix->size[1];
  NofCustomers = RewardMatrix->size[0];
  idx = Item2Customer->size[0] * Item2Customer->size[1];
  Item2Customer->size[0] = 1;
  Item2Customer->size[1] = RewardMatrix->size[1];
  emxEnsureCapacity_real_T(Item2Customer, idx);
  k = RewardMatrix->size[1];
  for (idx = 0; idx < k; idx++) {
    Item2Customer->data[idx] = 0.0;
  }

  idx = Customer2Item->size[0] * Customer2Item->size[1];
  Customer2Item->size[0] = 1;
  Customer2Item->size[1] = RewardMatrix->size[0];
  emxEnsureCapacity_real_T(Customer2Item, idx);
  k = RewardMatrix->size[0];
  for (idx = 0; idx < k; idx++) {
    Customer2Item->data[idx] = 0.0;
  }

  emxInit_boolean_T(&x, 2);
  emxInit_real_T(&varargin_1, 2);
  emxInit_int32_T(&iindx, 2);
  emxInit_real_T1(&b_RewardMatrix, 1);
  do {
    exitg1 = 0;
    idx = x->size[0] * x->size[1];
    x->size[0] = 1;
    x->size[1] = Customer2Item->size[1];
    emxEnsureCapacity_boolean_T(x, idx);
    k = Customer2Item->size[0] * Customer2Item->size[1];
    for (idx = 0; idx < k; idx++) {
      x->data[idx] = (Customer2Item->data[idx] == 0.0);
    }

    k = x->size[1];
    if (1 < k) {
      k = 1;
    }

    idx = 0;
    cindx = 0;
    exitg2 = false;
    while ((!exitg2) && (cindx + 1 <= x->size[1])) {
      if (x->data[cindx]) {
        idx = 1;
        exitg2 = true;
      } else {
        cindx++;
      }
    }

    if (k == 1) {
      if (idx == 0) {
        k = 0;
      }
    } else {
      k = !(1 > idx);
    }

    if (!(k == 0)) {
      if (NofItems == 1) {
        // if there is only one item
        unnamed_idx_1 = (unsigned int)RewardMatrix->size[1];
        idx = iindx->size[0] * iindx->size[1];
        iindx->size[0] = 1;
        iindx->size[1] = (int)unnamed_idx_1;
        emxEnsureCapacity_int32_T(iindx, idx);
        k = (int)unnamed_idx_1;
        for (idx = 0; idx < k; idx++) {
          iindx->data[idx] = 1;
        }

        n = RewardMatrix->size[0];
        for (i = 0; i + 1 <= RewardMatrix->size[1]; i++) {
          k = i * n;
          ixstart = i * n + 1;
          idx = k + n;
          mtmp = RewardMatrix->data[k];
          itmp = 1;
          if (n > 1) {
            cindx = 1;
            if (rtIsNaN(RewardMatrix->data[k])) {
              k = ixstart;
              exitg2 = false;
              while ((!exitg2) && (k + 1 <= idx)) {
                cindx++;
                ixstart = k + 1;
                if (!rtIsNaN(RewardMatrix->data[k])) {
                  mtmp = RewardMatrix->data[k];
                  itmp = cindx;
                  exitg2 = true;
                } else {
                  k++;
                }
              }
            }

            if (ixstart < idx) {
              while (ixstart + 1 <= idx) {
                cindx++;
                if (RewardMatrix->data[ixstart] > mtmp) {
                  mtmp = RewardMatrix->data[ixstart];
                  itmp = cindx;
                }

                ixstart++;
              }
            }
          }

          iindx->data[i] = itmp;
        }

        idx = Item2Customer->size[0] * Item2Customer->size[1];
        Item2Customer->size[0] = 1;
        Item2Customer->size[1] = iindx->size[1];
        emxEnsureCapacity_real_T(Item2Customer, idx);
        k = iindx->size[0] * iindx->size[1];
        for (idx = 0; idx < k; idx++) {
          Item2Customer->data[idx] = iindx->data[idx];
        }

        // Assign the item to the best customer
        idx = iindx->size[0] * iindx->size[1];
        iindx->size[0] = 1;
        iindx->size[1] = Item2Customer->size[1];
        emxEnsureCapacity_int32_T(iindx, idx);
        k = Item2Customer->size[0] * Item2Customer->size[1];
        for (idx = 0; idx < k; idx++) {
          iindx->data[idx] = (int)Item2Customer->data[idx];
        }

        k = iindx->size[0] * iindx->size[1];
        for (idx = 0; idx < k; idx++) {
          Customer2Item->data[iindx->data[idx] - 1] = 1.0;
        }

        // Assign the corresponding customer to the item
      } else {
        for (i = 0; i < NofCustomers; i++) {
          if (!(Customer2Item->data[i] != 0.0)) {
            k = RewardMatrix->size[1];
            idx = varargin_1->size[0] * varargin_1->size[1];
            varargin_1->size[0] = 1;
            varargin_1->size[1] = k;
            emxEnsureCapacity_real_T(varargin_1, idx);
            for (idx = 0; idx < k; idx++) {
              varargin_1->data[varargin_1->size[0] * idx] = RewardMatrix->data[i
                + RewardMatrix->size[0] * idx];
            }

            ixstart = 1;
            n = RewardMatrix->size[1];
            mtmp = RewardMatrix->data[i];
            itmp = 0;
            idx = RewardMatrix->size[1];
            if (idx > 1) {
              if (rtIsNaN(mtmp)) {
                k = 1;
                exitg2 = false;
                while ((!exitg2) && (k + 1 <= n)) {
                  ixstart = k + 1;
                  if (!rtIsNaN(varargin_1->data[k])) {
                    mtmp = varargin_1->data[k];
                    itmp = k;
                    exitg2 = true;
                  } else {
                    k++;
                  }
                }
              }

              idx = RewardMatrix->size[1];
              if (ixstart < idx) {
                while (ixstart + 1 <= n) {
                  if (varargin_1->data[ixstart] > mtmp) {
                    mtmp = varargin_1->data[ixstart];
                    itmp = ixstart;
                  }

                  ixstart++;
                }
              }
            }

            // find maximum element value and its index
            k = RewardMatrix->size[1];
            idx = varargin_1->size[0] * varargin_1->size[1];
            varargin_1->size[0] = 1;
            varargin_1->size[1] = k;
            emxEnsureCapacity_real_T(varargin_1, idx);
            for (idx = 0; idx < k; idx++) {
              varargin_1->data[varargin_1->size[0] * idx] = RewardMatrix->data[i
                + RewardMatrix->size[0] * idx];
            }

            ixstart = 1;
            n = RewardMatrix->size[1];
            b_mtmp = RewardMatrix->data[i];
            idx = RewardMatrix->size[1];
            if (idx > 1) {
              if (rtIsNaN(b_mtmp)) {
                k = 2;
                exitg2 = false;
                while ((!exitg2) && (k <= n)) {
                  ixstart = k;
                  if (!rtIsNaN(varargin_1->data[k - 1])) {
                    b_mtmp = varargin_1->data[k - 1];
                    exitg2 = true;
                  } else {
                    k++;
                  }
                }
              }

              idx = RewardMatrix->size[1];
              if (ixstart < idx) {
                while (ixstart + 1 <= n) {
                  if (varargin_1->data[ixstart] < b_mtmp) {
                    b_mtmp = varargin_1->data[ixstart];
                  }

                  ixstart++;
                }
              }
            }

            RewardMatrix->data[i + RewardMatrix->size[0] * itmp] = b_mtmp - 1.0;

            // make the maximum minimum to find second maximum
            k = RewardMatrix->size[1];
            idx = varargin_1->size[0] * varargin_1->size[1];
            varargin_1->size[0] = 1;
            varargin_1->size[1] = k;
            emxEnsureCapacity_real_T(varargin_1, idx);
            for (idx = 0; idx < k; idx++) {
              varargin_1->data[varargin_1->size[0] * idx] = RewardMatrix->data[i
                + RewardMatrix->size[0] * idx];
            }

            ixstart = 1;
            n = RewardMatrix->size[1];
            b_mtmp = RewardMatrix->data[i];
            idx = RewardMatrix->size[1];
            if (idx > 1) {
              if (rtIsNaN(b_mtmp)) {
                k = 2;
                exitg2 = false;
                while ((!exitg2) && (k <= n)) {
                  ixstart = k;
                  if (!rtIsNaN(varargin_1->data[k - 1])) {
                    b_mtmp = varargin_1->data[k - 1];
                    exitg2 = true;
                  } else {
                    k++;
                  }
                }
              }

              idx = RewardMatrix->size[1];
              if (ixstart < idx) {
                while (ixstart + 1 <= n) {
                  if (varargin_1->data[ixstart] > b_mtmp) {
                    b_mtmp = varargin_1->data[ixstart];
                  }

                  ixstart++;
                }
              }
            }

            // find the second maximum value and its index
            RewardMatrix->data[i + RewardMatrix->size[0] * itmp] = mtmp;

            // restore the maximum value
            Customer2Item->data[i] = itmp + 1;

            // Assign the customer the item
            if ((int)(unsigned int)Item2Customer->data[itmp] != 0) {
              // if item is already assigned
              Customer2Item->data[(int)(unsigned int)Item2Customer->data[itmp] -
                1] = 0.0;

              // unassign the corresponding customer
            }

            Item2Customer->data[itmp] = 1.0 + (double)i;

            // Assign the item to the customer
            k = RewardMatrix->size[0];
            mtmp = (mtmp - b_mtmp) + 0.05;
            idx = b_RewardMatrix->size[0];
            b_RewardMatrix->size[0] = k;
            emxEnsureCapacity_real_T1(b_RewardMatrix, idx);
            for (idx = 0; idx < k; idx++) {
              b_RewardMatrix->data[idx] = RewardMatrix->data[idx +
                RewardMatrix->size[0] * itmp] - mtmp;
            }

            k = b_RewardMatrix->size[0];
            for (idx = 0; idx < k; idx++) {
              RewardMatrix->data[idx + RewardMatrix->size[0] * itmp] =
                b_RewardMatrix->data[idx];
            }

            // reduce the item's value
          }
        }
      }
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_real_T(&b_RewardMatrix);
  emxFree_int32_T(&iindx);
  emxFree_real_T(&varargin_1);
  emxFree_boolean_T(&x);
}

//
// File trailer for auctionAlgorithm.cpp
//
// [EOF]
//
