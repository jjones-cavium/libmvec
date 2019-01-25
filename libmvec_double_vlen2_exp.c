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

/* Based on exp.c from https://github.com/ARM-software/optimized-routines */

#include <math.h>
#include <float.h>
#include <stdint.h>
#include <stdlib.h>
#include <ieee754.h>
#include <math-narrow-eval.h>
#include "math_config.h"
#include "libmvec_util.h"

#define N (1 << EXP_TABLE_BITS)
#define InvLn2N __exp_data.invln2N
#define NegLn2hiN __exp_data.negln2hiN
#define NegLn2loN __exp_data.negln2loN
#define Shift __exp_data.shift
#define T __exp_data.tab
#define C2 __exp_data.poly[5 - EXP_POLY_ORDER]
#define C3 __exp_data.poly[6 - EXP_POLY_ORDER]
#define C4 __exp_data.poly[7 - EXP_POLY_ORDER]
#define C5 __exp_data.poly[8 - EXP_POLY_ORDER]

#define CUTOFF 700.0

__AARCH64_VECTOR_PCS_ATTR static __Float64x2_t
__scalar_exp(__Float64x2_t x)
{
  return (__Float64x2_t) { exp(x[0]), exp(x[1]) };
}

__AARCH64_VECTOR_PCS_ATTR __Float64x2_t
_ZGVnN2v_exp(__Float64x2_t x)
{
  double h, z_0, z_1;
  __Float64x2_t g, scale_v, tail_v, tmp_v, r_v, r2_v, kd_v;
  __Float64x2_t NegLn2hiN_v, NegLn2loN_v, C2_v, C3_v, C4_v, C5_v;
  uint64_t ki_0, ki_1, idx_0, idx_1;
  uint64_t top_0, top_1, sbits_0, sbits_1;

  if (__glibc_unlikely(!__builtin_isnormal (x[0]) || !__builtin_isnormal (x[1])))
    return __scalar_exp (x);
  g = __builtin_aarch64_absv2df (x);
  h = __builtin_aarch64_reduc_smax_scal_v2df (g);
  if (__glibc_unlikely(h > CUTOFF))
    return __scalar_exp (x);

  z_0 = InvLn2N * x[0];
  z_1 = InvLn2N * x[1];
  ki_0 = converttoint (z_0);
  ki_1 = converttoint (z_1);

  idx_0 = 2 * (ki_0 % N);
  idx_1 = 2 * (ki_1 % N);
  top_0 = ki_0 << (52 - EXP_TABLE_BITS);
  top_1 = ki_1 << (52 - EXP_TABLE_BITS);
  sbits_0 = T[idx_0 + 1] + top_0;
  sbits_1 = T[idx_1 + 1] + top_1;

  kd_v = (__Float64x2_t) { roundtoint (z_0), roundtoint (z_1) };
  scale_v = (__Float64x2_t) { asdouble (sbits_0), asdouble (sbits_1) };
  tail_v = (__Float64x2_t) { asdouble (T[idx_0]), asdouble (T[idx_1]) };
  NegLn2hiN_v = (__Float64x2_t) { NegLn2hiN, NegLn2hiN };
  NegLn2loN_v = (__Float64x2_t) { NegLn2loN, NegLn2loN };
  C2_v = (__Float64x2_t) { C2, C2 };
  C3_v = (__Float64x2_t) { C3, C3 };
  C4_v = (__Float64x2_t) { C4, C4 };
  C5_v = (__Float64x2_t) { C5, C5 };

  r_v = x + kd_v * NegLn2hiN_v + kd_v * NegLn2loN_v;
  r2_v = r_v * r_v;
  tmp_v = tail_v + r_v + r2_v * (C2_v + r_v * C3_v) + r2_v * r2_v * (C4_v + r_v * C5_v);
  return scale_v + scale_v * tmp_v;
}
weak_alias (_ZGVnN2v_exp, _ZGVnN2v___exp_finite)
