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

/* Based on logf.c from https://github.com/ARM-software/optimized-routines */

#include <math.h>
#include <stdint.h>
#include <sysdeps/ieee754/flt-32/math_config.h>
#include "libmvec_util.h"

__AARCH64_VECTOR_PCS_ATTR __Float32x4_t
__scalar_logf (__Float32x4_t x)
{
  return (__Float32x4_t) { logf(x[0]), logf(x[1]), logf(x[2]), logf(x[3]) };
}

#define T __logf_data.tab
#define A __logf_data.poly
#define Ln2 __logf_data.ln2
#define N (1 << LOGF_TABLE_BITS)
#define OFF 0x3f330000

#define CALC_Z(__n) \
	ix_##__n = asuint (x_##__n); \
	tmp_##__n = ix_##__n - OFF; \
	i_##__n = (tmp_##__n >> (23 - LOGF_TABLE_BITS)) % N; \
	k_##__n = (int32_t) tmp_##__n >> 23; \
	iz_##__n = ix_##__n - (tmp_##__n & 0x1ff << 23); \
	invc_##__n = T[i_##__n].invc; \
	logc_##__n = T[i_##__n].logc; \
	z_##__n = (double_t) asfloat (iz_##__n)

__AARCH64_VECTOR_PCS_ATTR __Float32x4_t
_ZGVnN4v_logf(__Float32x4_t x)
{
  __Float64x2_t vz_0, vz_1, vinvc_0, vinvc_1, vlogc_0, vlogc_1;
  __Float64x2_t vk_0, vk_1, vLn2, va0, va1, va2, vone;
  __Float64x2_t r_0, r_1, y0_0, y0_1, r2_0, r2_1, y_0, y_1;
  __Float32x4_t result;
  double z_0, z_1, z_2, z_3;
  double invc_0, invc_1, invc_2, invc_3, logc_0, logc_1, logc_2, logc_3;

  uint32_t ix_0, ix_1, ix_2, ix_3, iz_0, iz_1, iz_2, iz_3;
  uint32_t tmp_0, tmp_1, tmp_2, tmp_3;
  int k_0, k_1, k_2, k_3, i_0, i_1, i_2, i_3;
  float f, x_0, x_1, x_2, x_3;

  x_0 = x[0];
  x_1 = x[1];
  x_2 = x[2];
  x_3 = x[3];

  f = __builtin_aarch64_reduc_smin_scal_v4sf (x);
  if (f < 1.3)
    return __scalar_logf (x);


  if (!isnormal (x_0) || !isnormal (x_1) || !isnormal (x_2) || !isnormal (x_3))
    return __scalar_logf (x);

  CALC_Z(0);
  CALC_Z(1);
  CALC_Z(2);
  CALC_Z(3);

  vz_0 = (__Float64x2_t) { z_0, z_1 };
  vz_1 = (__Float64x2_t) { z_2, z_3 };
  vinvc_0 = (__Float64x2_t) { invc_0, invc_1 };
  vinvc_1 = (__Float64x2_t) { invc_2, invc_3 };
  vone = (__Float64x2_t) { 1.0, 1.0 };
  vlogc_0 = (__Float64x2_t) { logc_0, logc_1 };
  vlogc_1 = (__Float64x2_t) { logc_2, logc_3 };
  vk_0 =  (__Float64x2_t) { (double) k_0, (double) k_1 };
  vk_1 =  (__Float64x2_t) { (double) k_2, (double) k_3 };
  vLn2 = (__Float64x2_t) { Ln2, Ln2 };
  va0 = (__Float64x2_t) { A[0], A[0] };
  va1 = (__Float64x2_t) { A[1], A[1] };
  va2 = (__Float64x2_t) { A[2], A[2] };

   /* log(x) = log1p(z/c-1) + log(c) + k*Ln2 */
  r_0 = vz_0 * vinvc_0 - vone;
  r_1 = vz_1 * vinvc_1 - vone;
  y0_0 = vlogc_0 + vk_0 * vLn2;
  y0_1 = vlogc_1 + vk_1 * vLn2;
  /* Pipelined polynomial evaluation to approximate log1p(r).  */
  r2_0 = r_0 * r_0;
  r2_1 = r_1 * r_1;
  y_0 = va1 * r_0 + va2;
  y_1 = va1 * r_1 + va2;
  y_0 = va0 * r2_0 + y_0;
  y_1 = va0 * r2_1 + y_1;
  y_0 = y_0 * r2_0 + (y0_0 + r_0);
  y_1 = y_1 * r2_1 + (y0_1 + r_1);
  result = pack_and_trunc (y_0, y_1);
  return result;
}
weak_alias (_ZGVnN4v_logf, _ZGVnN4v___logf_finite)
