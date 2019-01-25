/* Copyright (c) 2018, Marvell Technology Group Ltd.
 * SPDX-License-Identifier: Apache-2.0
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/* Based on log2.c from https://github.com/ARM-software/optimized-routines */

#include <math.h>
#include <math_private.h>
#include <stdint.h>
#include <stdlib.h>
#include <libc-symbols.h>
#include "libmvec_util.h"
#include "math_config.h"

#define T __log2_data.tab
#define T2 __log2_data.tab2
#define B __log2_data.poly1
#define A __log2_data.poly
#define InvLn2hi __log2_data.invln2hi
#define InvLn2lo __log2_data.invln2lo
#define N (1 << LOG2_TABLE_BITS)
#define OFF 0x3fe6000000000000

__AARCH64_VECTOR_PCS_ATTR static __Float64x2_t
__scalar_log2(__Float64x2_t x)
{
  return (__Float64x2_t) { log2(x[0]), log2(x[1]) };
}

__AARCH64_VECTOR_PCS_ATTR __Float64x2_t
_ZGVnN2v_log2(__Float64x2_t x)
{
  double_t z_0, z_1;
  double_t invc_0, invc_1, logc_0, logc_1;
  double_t x_0, x_1, kd_0, kd_1;
  __Float64x2_t r_v, r2_v, r4_v, p_v, y_v, z_v, kd_v;
  __Float64x2_t hi_v, lo_v, invc_v, logc_v;
  __Float64x2_t negone_v;
  __Float64x2_t A0_v, A1_v, A2_v, A3_v, A4_v, A5_v;
  __Float64x2_t InvLn2hi_v, InvLn2lo_v, t1_v, t2_v, t3_v;
  uint64_t ix_0, ix_1, iz_0, iz_1, tmp_0, tmp_1;
  int i_0, i_1, k_0, k_1;

  x_0 = x[0];
  x_1 = x[1];

  if (__glibc_unlikely(!__builtin_isnormal (x_0) || !__builtin_isnormal (x_1)))
    return __scalar_log2 (x);
  if (x_0 <= 1.32 || x_1 <= 1.32)
    return __scalar_log2 (x);

  ix_0 = asuint64 (x_0);
  ix_1 = asuint64 (x_1);
  tmp_0 = ix_0 - OFF;
  tmp_1 = ix_1 - OFF;
  i_0 = (tmp_0 >> (52 - LOG2_TABLE_BITS)) % N;
  i_1 = (tmp_1 >> (52 - LOG2_TABLE_BITS)) % N;
  k_0 = (int64_t) tmp_0 >> 52;
  k_1 = (int64_t) tmp_1 >> 52;
  iz_0 = ix_0 - (tmp_0 & 0xfffULL << 52);
  iz_1 = ix_1 - (tmp_1 & 0xfffULL << 52);
  invc_0 = T[i_0].invc;
  invc_1 = T[i_1].invc;
  logc_0 = T[i_0].logc;
  logc_1 = T[i_1].logc;
  z_0 = asdouble (iz_0);
  z_1 = asdouble (iz_1);
  kd_0 = (double_t) k_0;
  kd_1 = (double_t) k_1;

  z_v = (__Float64x2_t) { z_0, z_1 };
  invc_v = (__Float64x2_t) { invc_0, invc_1 };
  negone_v = (__Float64x2_t) { -1.0, -1.0 };
  InvLn2hi_v = (__Float64x2_t) { InvLn2hi, InvLn2hi };
  InvLn2lo_v = (__Float64x2_t) { InvLn2lo, InvLn2lo };
  logc_v = (__Float64x2_t) { logc_0, logc_1 };
  kd_v = (__Float64x2_t) { kd_0, kd_1 };

  A0_v = (__Float64x2_t) { A[0], A[0] };
  A1_v = (__Float64x2_t) { A[1], A[1] };
  A2_v = (__Float64x2_t) { A[2], A[2] };
  A3_v = (__Float64x2_t) { A[3], A[3] };
  A4_v = (__Float64x2_t) { A[4], A[4] };
  A5_v = (__Float64x2_t) { A[5], A[5] };

  r_v = __builtin_aarch64_fmav2df (z_v, invc_v, negone_v);
  t1_v = r_v * InvLn2hi_v;
  t2_v = r_v * InvLn2lo_v + __builtin_aarch64_fmav2df (r_v, InvLn2hi_v, -t1_v);
  t3_v = kd_v + logc_v;
  hi_v = t3_v + t1_v;
  lo_v = t3_v - hi_v + t1_v + t2_v;
  r2_v = r_v * r_v;
  r4_v = r2_v * r2_v;
  p_v = A0_v + r_v * A1_v + r2_v * (A2_v + r_v * A3_v) + r4_v * (A4_v + r_v * A5_v);
  y_v = lo_v + r2_v * p_v + hi_v;
  return y_v;
}
weak_alias (_ZGVnN2v_log2, _ZGVnN2v___log2_finite)
