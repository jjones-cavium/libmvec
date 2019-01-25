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

extern __AARCH64_VECTOR_PCS_ATTR __Float64x2_t _ZGVnN2v_exp2 (__Float64x2_t);
extern __AARCH64_VECTOR_PCS_ATTR __Float64x2_t _ZGVnN2v_log2 (__Float64x2_t);

#define CUTOFF 125.0

  /* The log2 version seems more precise, the log version fails one test
     in the glibc testsuite.   (0x1.430d4cp+0 ** 0x5.0e462p+4)  */
  /* pow(x,y) = e^(y * log(x))  */
  /* pow(x,y) = 2^(y * log2(x))  */

__AARCH64_VECTOR_PCS_ATTR static __Float64x2_t
__scalar_pow(__Float64x2_t x, __Float64x2_t y)
{
  return (__Float64x2_t) { pow(x[0],y[0]), pow(x[1],y[1]) };
}

__AARCH64_VECTOR_PCS_ATTR __Float64x2_t
_ZGVnN2vv_pow(__Float64x2_t x, __Float64x2_t y)
{
  double c,d,e,f;

  c = x[0];
  d = x[1];
  e = y[0];
  f = y[1];
  if (__glibc_unlikely(!__builtin_isnormal (c) || !__builtin_isnormal (d)
		|| !__builtin_isnormal (e) || !__builtin_isnormal (f)))
    return __scalar_pow (x, y);

  if (c < 0.0 || d < 0.0 || e < 0.0 || f < 0.0)
    return __scalar_pow (x, y);

  if (c > CUTOFF || d > CUTOFF || e > CUTOFF || f > CUTOFF)
    return __scalar_pow (x, y);

  return (_ZGVnN2v_exp2 (y * _ZGVnN2v_log2 (x)));
}
weak_alias (_ZGVnN2vv_pow, _ZGVnN2vv___pow_finite)
