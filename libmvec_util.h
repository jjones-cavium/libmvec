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

#include <stdint.h>

#ifdef COMPILER_SUPPORTS_SIMD_ABI
#  define __AARCH64_VECTOR_PCS_ATTR __attribute__((aarch64_vector_pcs))
#else
#  define __AARCH64_VECTOR_PCS_ATTR
#endif

static __always_inline
__Float64x2_t get_lo_and_extend (__Float32x4_t x)
{
	__Uint64x2_t  tmp1 = (__Uint64x2_t) x;
#ifdef BIG_ENDI
	uint64_t      tmp2 = (uint64_t) tmp1[1];
#else
	uint64_t      tmp2 = (uint64_t) tmp1[0];
#endif
	return __builtin_aarch64_float_extend_lo_v2df ((__Float32x2_t) tmp2);
}

static __always_inline
__Float64x2_t get_hi_and_extend (__Float32x4_t x)
{
	return __builtin_aarch64_vec_unpacks_hi_v4sf (x);
}

static __always_inline
__Float32x4_t pack_and_trunc (__Float64x2_t x, __Float64x2_t y)
{
        __Float32x2_t xx = __builtin_aarch64_float_truncate_lo_v2sf (x);
        __Float32x2_t yy = __builtin_aarch64_float_truncate_lo_v2sf (y);
        return (__builtin_aarch64_combinev2sf (xx, yy));
}
