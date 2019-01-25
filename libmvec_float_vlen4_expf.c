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

/* Based on expf.c from https://github.com/ARM-software/optimized-routines */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <sysdeps/ieee754/flt-32/math_config.h>
#include "libmvec_util.h"

#define N (1 << EXP2F_TABLE_BITS)
#define LIMIT 80.0

#define InvLn2N __exp2f_data.invln2_scaled
#define T __exp2f_data.tab
#define C __exp2f_data.poly_scaled
#define SHIFT __exp2f_data.shift

__AARCH64_VECTOR_PCS_ATTR static __Float32x4_t
__scalar_expf (__Float32x4_t x)
{
  return (__Float32x4_t) { expf(x[0]), expf(x[1]), expf(x[2]), expf(x[3]) };
}

__AARCH64_VECTOR_PCS_ATTR __Float32x4_t
_ZGVnN4v_expf(__Float32x4_t x)
{
  __Float32x4_t g, result;
  __Float64x2_t xd_0, xd_1, vInvLn2N, z_0, z_1, vkd_0, vkd_1, r_0, r_1;
  __Float64x2_t vs_0, vs_1, c0, c1, c2, y_0, y_1, r2_0, r2_1, one;
  uint64_t ki_0, ki_1, ki_2, ki_3, t_0, t_1, t_2, t_3;
  double s_0, s_1, s_2, s_3, kd_0, kd_1, kd_2, kd_3;
  float f;

  g = __builtin_aarch64_absv4sf (x);
  f = __builtin_aarch64_reduc_smax_scal_v4sf (g);
  if (f > LIMIT)
    return __scalar_expf (x);
  
  if (!isnormal (x[0]) || !isnormal (x[1])
      || !isnormal (x[2]) || !isnormal (x[3]))
    return __scalar_expf (x);

  xd_0 = get_lo_and_extend (x);
  xd_1 = get_hi_and_extend (x);

  vInvLn2N = (__Float64x2_t) { InvLn2N, InvLn2N };
  /* x*N/Ln2 = k + r with r in [-1/2, 1/2] and int k.  */
  z_0 = vInvLn2N * xd_0;
  z_1 = vInvLn2N * xd_1;

    /* Round and convert z to int, the result is in [-150*N, 128*N] and
     ideally ties-to-even rule is used, otherwise the magnitude of r
     can be bigger which gives larger approximation error.  */
#if 1
  kd_0 = (double) (z_0[0] + SHIFT);
  kd_1 = (double) (z_0[1] + SHIFT);
  kd_2 = (double) (z_1[0] + SHIFT);
  kd_3 = (double) (z_1[1] + SHIFT);
  ki_0 = asuint64 (kd_0);
  ki_1 = asuint64 (kd_1);
  ki_2 = asuint64 (kd_2);
  ki_3 = asuint64 (kd_3);
  kd_0 -= SHIFT;
  kd_1 -= SHIFT;
  kd_2 -= SHIFT;
  kd_3 -= SHIFT;

  vkd_0 = (__Float64x2_t) {kd_0, kd_1 };
  vkd_1 = (__Float64x2_t) {kd_2, kd_3 };
  r_0 = z_0 - vkd_0;
  r_1 = z_1 - vkd_1;
#else
  kd_0 = __builtin_aarch64_roundv2df (z_0);
  kd_1 = __builtin_aarch64_roundv2df (z_1);
  r_0 = z_0 - kd_0;
  r_1 = z_1 - kd_1;

  ki_0 = (long) kd_0[0];
  ki_1 = (long) kd_0[1];
  ki_2 = (long) kd_1[0];
  ki_3 = (long) kd_1[1];
#endif

  /* exp(x) = 2^(k/N) * 2^(r/N) ~= s * (C0*r^3 + C1*r^2 + C2*r + 1) */
  t_0 = T[ki_0 % N];
  t_1 = T[ki_1 % N];
  t_2 = T[ki_2 % N];
  t_3 = T[ki_3 % N];
  t_0 += ki_0 << (52 - EXP2F_TABLE_BITS);
  t_1 += ki_1 << (52 - EXP2F_TABLE_BITS);
  t_2 += ki_2 << (52 - EXP2F_TABLE_BITS);
  t_3 += ki_3 << (52 - EXP2F_TABLE_BITS);
  s_0 = asdouble (t_0);
  s_1 = asdouble (t_1);
  s_2 = asdouble (t_2);
  s_3 = asdouble (t_3);

  vs_0 = (__Float64x2_t) { s_0, s_1 };
  vs_1 = (__Float64x2_t) { s_2, s_3 };
  c0 = (__Float64x2_t) { C[0], C[0] };
  c1 = (__Float64x2_t) { C[1], C[1] };
  c2 = (__Float64x2_t) { C[2], C[2] };
  one = (__Float64x2_t) { 1.0, 1.0 };

  z_0 = c0 * r_0 + c1;
  z_1 = c0 * r_1 + c1;
  r2_0 = r_0 * r_0;
  r2_1 = r_1 * r_1;
  y_0 = c2 * r_0 + one;
  y_1 = c2 * r_1 + one;
  y_0 = z_0 * r2_0 + y_0;
  y_1 = z_1 * r2_1 + y_1;
  y_0 = y_0 * vs_0;
  y_1 = y_1 * vs_1;
  result = pack_and_trunc (y_0, y_1);
  return result;
}
weak_alias (_ZGVnN4v_expf, _ZGVnN4v___expf_finite)
