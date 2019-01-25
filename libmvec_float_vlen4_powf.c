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

#include <math.h>
#include "libmvec_util.h"

extern __Float32x4_t _ZGVnN4v_exp2f(__Float32x4_t);
extern __Float32x4_t _ZGVnN4v_log2f(__Float32x4_t);

__AARCH64_VECTOR_PCS_ATTR static __Float32x4_t
__scalar_powf(__Float32x4_t x, __Float32x4_t y)
{
  return (__Float32x4_t) { powf(x[0],y[0]), powf(x[1],y[1]),
			   powf(x[2],y[2]), powf(x[3],y[3]) };
}

#define CUTOFF 80.0

__AARCH64_VECTOR_PCS_ATTR __Float32x4_t
_ZGVnN4vv_powf(__Float32x4_t x, __Float32x4_t y)
{
  float x_0, x_1, x_2, x_3, y_0, y_1, y_2, y_3;
  float f, g;

  x_0 = x[0];
  x_1 = x[1];
  x_2 = x[2];
  x_3 = x[3];
  y_0 = y[0];
  y_1 = y[1];
  y_2 = y[2];
  y_3 = y[3];

  if (__glibc_unlikely(!__builtin_isnormal (x_0) || !__builtin_isnormal (x_1)
		    || !__builtin_isnormal (x_2) || !__builtin_isnormal (x_3)))
    return __scalar_powf (x, y);

  if (__glibc_unlikely(!__builtin_isnormal (y_0) || !__builtin_isnormal (y_1)
		    || !__builtin_isnormal (y_2) || !__builtin_isnormal (y_3)))
    return __scalar_powf (x, y);

  f = __builtin_aarch64_reduc_smin_scal_v4sf (x);
  g = __builtin_aarch64_reduc_smin_scal_v4sf (y);

  if (f < 0.0 || g < 0.0)
    return __scalar_powf (x, y);

  f = __builtin_aarch64_reduc_smax_scal_v4sf (x);
  g = __builtin_aarch64_reduc_smax_scal_v4sf (y);

  if (f > CUTOFF || g > CUTOFF)
    return __scalar_powf (x, y);

  return (_ZGVnN4v_exp2f (y * _ZGVnN4v_log2f (x)));
}
weak_alias (_ZGVnN4vv_powf, _ZGVnN4vv___powf_finite)
