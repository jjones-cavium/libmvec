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

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "endian.h"
#include "libmvec_double_sinetable.h"
#include "libmvec_util.h"

//
// main body of routine
//

__AARCH64_VECTOR_PCS_ATTR
static inline __Float64x2_t _sine_kernel(__Float64x2_t x, double *tbl, int sym)
{
  __Float64x2_t result;
  __Float64x2_t *vmod_ptr, vmod;
  __Float64x2_t *modulus;
  __Float64x2_t m0,m1,m2;
  __Float64x2_t t0, t1, a, a0, a1, a2, x0, x1, x02, c0, r1, k;
  __Float64x2_t hiref, loref, tbl2, tbl3;
  __Float64x2_t sinpoly_0, sinpoly_1, sinpoly_2, sinpoly_3, sinpoly_4;
  __Float64x2_t cospoly_0, cospoly_1, cospoly_2, cospoly_3, cospoly_4;
  double *sinpoly, *cospoly;
  double hiref_0, hiref_1, loref_0, loref_1, tbl2_0, tbl2_1, tbl3_0, tbl3_1;
  int sign_0, sign_1;
  unsigned int tablebase_0, tablebase_1;

  sinpoly = (double *)_sin_poly;
  cospoly = (double *)_cos_poly;
  sinpoly_0 = (__Float64x2_t) { sinpoly[0], sinpoly[0] };
  sinpoly_1 = (__Float64x2_t) { sinpoly[1], sinpoly[1] };
  sinpoly_2 = (__Float64x2_t) { sinpoly[2], sinpoly[2] };
  sinpoly_3 = (__Float64x2_t) { sinpoly[3], sinpoly[3] };
  sinpoly_4 = (__Float64x2_t) { sinpoly[4], sinpoly[4] };
  cospoly_0 = (__Float64x2_t) { cospoly[0], cospoly[0] };
  cospoly_1 = (__Float64x2_t) { cospoly[1], cospoly[1] };
  cospoly_2 = (__Float64x2_t) { cospoly[2], cospoly[2] };
  cospoly_3 = (__Float64x2_t) { cospoly[3], cospoly[3] };
  cospoly_4 = (__Float64x2_t) { cospoly[4], cospoly[4] };
  vmod_ptr = (__Float64x2_t *)_vmod;
  vmod = *vmod_ptr;
  modulus = (__Float64x2_t *)_modulus; 

  m0 = modulus[0];
  m1 = modulus[1];
  m2 = modulus[2];

  sign_0 = x[0] > 0 ? 0 : 1;
  sign_1 = x[1] > 0 ? 0 : 1;

  a = __builtin_aarch64_absv2df (x);
  a1 = a * vmod;  /* vmod is 1/(2*PI/32) */
  k = __builtin_aarch64_roundv2df (a1);
  tablebase_0 = (unsigned int) k[0];
  tablebase_1 = (unsigned int) k[1];
  t0 = k * m0;
  /* We need to explicitly use fma for the extra precision. */
  t1 = __builtin_aarch64_fmav2df (k, m0, -t0);
  a1 = k * m1;
  /* We need to explicitly use fma for the extra precision. */
  a2 =  __builtin_aarch64_fmav2df (k, m1, -a1);
  a0 = a - t0;
  a1 = a1 + t1; /* add remainder from 1st term to 2nd term */
  x0 = a0 - a1;
  x1 = a0 - x0;
  x1 = x1 - a1;
  r1 = x1 - a2;
  r1 = r1 - k * m2;
  
  if (sign_0)
    tablebase_0 += sym;
  if (sign_1)
    tablebase_1 += sym;

  tablebase_0 = (tablebase_0 << 2) & 0x7c;
  tablebase_1 = (tablebase_1 << 2) & 0x7c;
  hiref_0 = tbl[tablebase_0];
  loref_0 = tbl[tablebase_0 + 1];
  hiref_1 = tbl[tablebase_1];
  loref_1 = tbl[tablebase_1 + 1];
  hiref = (__Float64x2_t) { hiref_0, hiref_1 };
  loref = (__Float64x2_t) { loref_0, loref_1 };
  c0 = x0 - hiref;
  x1 = x0 - c0;
  x1 = x1 - hiref;
  x0 = c0;
  x1 = x1 + r1;
  x1 = x1 - loref;
  x02 = x0 * x0; /* reduced x^2 for poly */
  t0 = x02 * sinpoly_4 + sinpoly_3;
  t0 = x02 * t0 + sinpoly_2;
  t0 = x02 * t0 + sinpoly_1;
  t0 = x02 * t0 + sinpoly_0;
  t0 = t0 * x02;
  t0 = x0 * t0 + x1;
  t0 = t0 + x0;
  tbl3_0 = tbl[tablebase_0+3];
  tbl3_1 = tbl[tablebase_1+3];
  tbl3 = (__Float64x2_t) { tbl3_0, tbl3_1 };
  t0 = t0 * tbl3;                         /* cos of ref, gives cos(a)*sin(b) */
  t1 = x02 * cospoly_4 + cospoly_3;
  t1 = x02 * t1 + cospoly_2;
  t1 = x02 * t1 + cospoly_1;
  t1 = x02 * t1 + cospoly_0;
  t1 = x02 * t1;
  tbl2_0 = tbl[tablebase_0+2];
  tbl2_1 = tbl[tablebase_1+2];
  tbl2 = (__Float64x2_t) { tbl2_0, tbl2_1 };
  t1 = t1 * tbl2;                         /* sin of ref */
  result = t1 + t0;
  result = result + tbl2;                 /* add sine */
  return(result);
}

#define CUTOFF 1000.00

//
// sine entry point
//

__AARCH64_VECTOR_PCS_ATTR
__Float64x2_t _ZGVnN2v_sin(__Float64x2_t x)
{
__Float64x2_t result, c;
double d,e;
double *ptr;
int sym;

  c = __builtin_aarch64_absv2df (x);
  d = __builtin_aarch64_reduc_smax_scal_v2df (c);
  e = __builtin_aarch64_reduc_smin_scal_v2df (c);
  
  /* This algorithm is inexact for large numbers.  */
  if (d > CUTOFF)
    return (__Float64x2_t) { sin(x[0]), sin(x[1]) };

  /* _sine_kernel returns +0 for sin(-0) which is wrong.  */
  if (e == 0)
    return (__Float64x2_t) { sin(x[0]), sin(x[1]) };

  ptr = (double *)_sin_table;
  sym = 1<<4;
  result = _sine_kernel(x,ptr,sym);
  return (result);
}
weak_alias (_ZGVnN2v_sin, _ZGVnN2v___sin_finite)

//
// cosine entry point
//

__AARCH64_VECTOR_PCS_ATTR
__Float64x2_t _ZGVnN2v_cos(__Float64x2_t x)
{
__Float64x2_t result;
double c;
double *ptr;
int sym;

  c = __builtin_aarch64_reduc_smax_scal_v2df (__builtin_aarch64_absv2df (x));
  if (c > CUTOFF)
    return (__Float64x2_t) { cos(x[0]), cos(x[1]) };

  ptr = (double *)_cos_table;
  sym = 0;
  result = _sine_kernel(x,ptr,sym);
  return (result);
}
weak_alias (_ZGVnN2v_cos, _ZGVnN2v___cos_finite)
