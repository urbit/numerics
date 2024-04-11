/// @file

#include "jets/q.h"
#include "jets/w.h"

#include "noun.h"
#include "softfloat.h"
#include "softblas.h"

#include <assert.h>
#include <math.h>  // for pow()
#include <stdio.h>

  union half {
    float16_t h;
    c3_w c;
  };

  union sing {
    float32_t s;
    c3_w c;
  };

  union doub {
    float64_t d;
    c3_d c;
  };

  union quad {
    float128_t q;
    c3_d c[2];
  };

  //  $?(%n %u %d %z %a)
  static inline void
  _set_rounding(c3_w a)
  {
    // We could use SoftBLAS set_rounding() to set the SoftFloat
    // mode as well, but it's more explicit to do it here since
    // we may use SoftFloat in any given Lagoon jet and we want
    // you, dear developer, to see it set here.
    switch ( a )
    {
    default:
      u3m_bail(c3__fail);
      break;
    // %n - near
    case c3__n:
      softfloat_roundingMode = softfloat_round_near_even;
      softblas_roundingMode = 'n';
      break;
    // %z - zero
    case c3__z:
      softfloat_roundingMode = softfloat_round_minMag;
      softblas_roundingMode = 'z';
      break;
    // %u - up
    case c3__u:
      softfloat_roundingMode = softfloat_round_max;
      softblas_roundingMode = 'u';
      break;
    // %d - down
    case c3__d:
      softfloat_roundingMode = softfloat_round_min;
      softblas_roundingMode = 'd';
      break;
    // %a - away
    case c3__a:
      softfloat_roundingMode = softfloat_round_near_maxMag;
      softblas_roundingMode = 'a';
      break;
    }
  }

/* length of shape = x * y * z * w * ...
*/
  static inline c3_d _get_length(u3_noun shape)
  {
    c3_d len = 1;
    while (u3_nul != shape) {
      len = len * u3x_atom(u3h(shape));
      shape = u3t(shape);
    }
    return len;
  }

/* get dims from shape as array [x y z w ...]
*/
  static inline c3_d* _get_dims(u3_noun shape)
  {
    u3_atom len = u3qb_lent(shape);
    c3_d len_d = u3r_chub(0, len);
    c3_d* dims = (c3_d*)u3a_malloc(len_d*sizeof(c3_d));
    for (c3_d i = 0; i < len_d; i++) {
      dims[i] = u3r_chub(0, u3x_atom(u3h(shape)));
      shape = u3t(shape);
    }
    u3z(len);
    return dims;
  }

/* 
*/
  static inline size_t _get_array_length(c3_d* array)
  {
    size_t n = sizeof(array)/sizeof(array[0]);
    for (size_t i = 0; i < n; i++) {
      fprintf(stderr, "%x ", array[i]);
    }
    fprintf(stderr, " => %x \n", n);
    return n;
  }

/* add - axpy = 1*x+y
*/
  u3_noun
  u3qf_la_add_real(u3_noun x_data,
                   u3_noun y_data,
                   u3_noun shape,
                   u3_noun bloq,
                   u3_noun kind,
                   u3_noun fxp
                   )
  {
    //  Fence on valid bloq size.
    if (bloq < 4 || bloq > 7) {
      return u3_none;
    }

    //  Unpack the data as a byte array.  We assume total length < 2**64.
    // len_x is length in base units
    c3_d len_x = _get_length(shape);

    // siz_x is length in bytes
    c3_d siz_x = len_x * pow(2, bloq-3);
    // fprintf(stderr, "len_x: %d 0x%x units\r\n", len_x, len_x);
    // fprintf(stderr, "siz_x: %d 0x%x bytes\r\n", siz_x, siz_x);

    // x_bytes is the data array (w/o leading 0x1)
    c3_y* x_bytes = (c3_y*)u3a_malloc(siz_x*sizeof(c3_y));
    u3r_bytes(0, siz_x, x_bytes, x_data);

    // y_bytes is the data array (w/ leading 0x1, skipped by ?axpy)
    c3_y* y_bytes = (c3_y*)u3a_malloc((siz_x+1)*sizeof(c3_y));
    u3r_bytes(0, siz_x+1, y_bytes, y_data);

    // r_data is the result noun of [data]
    u3_noun r_data;

    //  Switch on the block size.
    switch (bloq) {
      case 4:
        haxpy(len_x, (float16_t){SB_REAL16_ONE}, (float16_t*)x_bytes, 1, (float16_t*)y_bytes, 1);
        break;

      case 5:
        saxpy(len_x, (float32_t){SB_REAL32_ONE}, (float32_t*)x_bytes, 1, (float32_t*)y_bytes, 1);
        break;

      case 6:
        daxpy(len_x, (float64_t){SB_REAL64_ONE}, (float64_t*)x_bytes, 1, (float64_t*)y_bytes, 1);
        break;

      case 7:
        qaxpy(len_x, (float128_t){SB_REAL128L_ONE,SB_REAL128U_ONE}, (float128_t*)x_bytes, 1, (float128_t*)y_bytes, 1);
        break;
    }

    r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), y_bytes);

    // fprintf(stderr, "y = [\r\n");
    // // fflush(stdout);
    // fflush(stderr);
    // for (c3_w i = 0; i <= siz_x; i++) {
    //   fprintf(stderr, ">>>>%d: ", i);
    //   fprintf(stderr, "%0x \r\n", y_bytes[siz_x-i]);
    // }
    // fprintf(stderr, "]\r\n");

    //  Clean up and return.
    u3a_free(x_bytes);
    u3a_free(y_bytes);

    fprintf(stderr, "here-end\r\n");
    return u3nc(u3nq(shape, bloq, kind, fxp), r_data);
  }

/* sub - axpy = -1*y+x
*/
  u3_noun
  u3qf_la_sub_real(u3_noun x_data,
                   u3_noun y_data,
                   u3_noun shape,
                   u3_noun bloq)
  {
    //  Unpack the data as a byte array.  We assume total length < 2**64.
    c3_d len_x = _get_length(shape);
    c3_d siz_x = len_x * pow(2, bloq - 3);
    c3_y* x_bytes = (c3_y*)u3a_malloc(siz_x*sizeof(c3_y));
    u3r_bytes(0, siz_x, x_bytes, y_data);  // XXX
    c3_y* y_bytes = (c3_y*)u3a_malloc((siz_x+1)*sizeof(c3_y));
    u3r_bytes(0, siz_x+1, y_bytes, x_data);  // XXX

    u3_noun r_data;

    //  Switch on the block size.
    switch (bloq) {
      case 4:
        haxpy(len_x, (float16_t){SB_REAL16_NEGONE}, (float16_t*)x_bytes, 1, (float16_t*)y_bytes, 1);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), y_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return r_data;

      case 5:
        saxpy(len_x, (float32_t){SB_REAL32_NEGONE}, (float32_t*)x_bytes, 1, (float32_t*)y_bytes, 1);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), y_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return r_data;

      case 6:
        daxpy(len_x, (float64_t){SB_REAL64_NEGONE}, (float64_t*)x_bytes, 1, (float64_t*)y_bytes, 1);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), y_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return r_data;

      case 7:
        qaxpy(len_x, (float128_t){SB_REAL128L_NEGONE,SB_REAL128U_NEGONE}, (float128_t*)x_bytes, 1, (float128_t*)y_bytes, 1);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), y_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return r_data;

      default:
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return u3_none;
    }
  }

/* mul - x.*y
   elementwise multiplication
*/
  u3_noun
  u3qf_la_mul_real(u3_noun x_data,
                   u3_noun y_data,
                   u3_noun shape,
                   u3_noun bloq)
  {
    //  Unpack the data as a byte array.  We assume total length < 2**64.
    c3_d len_x = _get_length(shape);
    c3_d siz_x = len_x * pow(2, bloq - 3);
    c3_y* x_bytes = (c3_y*)u3a_malloc(siz_x*sizeof(c3_y));
    u3r_bytes(0, siz_x, x_bytes, x_data);
    c3_y* y_bytes = (c3_y*)u3a_malloc((siz_x+1)*sizeof(c3_y));
    u3r_bytes(0, siz_x+1, y_bytes, y_data);

    u3_noun r_data;

    //  Switch on the block size.
    switch (bloq) {
      case 4:
        for (c3_d i = 0; i < len_x; i++) {
          ((float16_t*)y_bytes)[i] = f16_mul(((float16_t*)x_bytes)[i], ((float16_t*)y_bytes)[i]);
        }

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), y_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return r_data;

      case 5:
        for (c3_d i = 0; i < len_x; i++) {
          ((float32_t*)y_bytes)[i] = f32_mul(((float32_t*)x_bytes)[i], ((float32_t*)y_bytes)[i]);
        }

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), y_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return r_data;

      case 6:
        for (c3_d i = 0; i < len_x; i++) {
          ((float64_t*)y_bytes)[i] = f64_mul(((float64_t*)x_bytes)[i], ((float64_t*)y_bytes)[i]);
        }

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), y_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return r_data;

      case 7:
        for (c3_d i = 0; i < len_x; i++) {
          f128M_mul(&(((float128_t*)y_bytes)[i]), &(((float128_t*)x_bytes)[i]), &(((float128_t*)y_bytes)[i]));
        }

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), y_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return r_data;

      default:
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return u3_none;
    }
  }

/* div - x/y
   elementwise division
*/
  u3_noun
  u3qf_la_div_real(u3_noun x_data,
                   u3_noun y_data,
                   u3_noun shape,
                   u3_noun bloq)
  {
    //  Unpack the data as a byte array.  We assume total length < 2**64.
    c3_d len_x = _get_length(shape);
    c3_d siz_x = len_x * pow(2, bloq - 3);
    c3_y* x_bytes = (c3_y*)u3a_malloc(siz_x*sizeof(c3_y));
    u3r_bytes(0, siz_x, x_bytes, x_data);
    c3_y* y_bytes = (c3_y*)u3a_malloc((siz_x+1)*sizeof(c3_y));
    u3r_bytes(0, siz_x+1, y_bytes, y_data);

    u3_noun r_data;

    //  Switch on the block size.
    switch (bloq) {
      case 4:
        for (c3_d i = 0; i < len_x; i++) {
          ((float16_t*)y_bytes)[i] = f16_div(((float16_t*)x_bytes)[i], ((float16_t*)y_bytes)[i]);
        }

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), y_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return r_data;

      case 5:
        for (c3_d i = 0; i < len_x; i++) {
          ((float32_t*)y_bytes)[i] = f32_div(((float32_t*)x_bytes)[i], ((float32_t*)y_bytes)[i]);
        }

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), y_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return r_data;

      case 6:
        for (c3_d i = 0; i < len_x; i++) {
          ((float64_t*)y_bytes)[i] = f64_div(((float64_t*)x_bytes)[i], ((float64_t*)y_bytes)[i]);
        }

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), y_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return r_data;

      case 7:
        for (c3_d i = 0; i < len_x; i++) {
          f128M_div(&(((float128_t*)y_bytes)[i]), &(((float128_t*)x_bytes)[i]), &(((float128_t*)y_bytes)[i]));
        }

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), y_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return r_data;

      default:
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return u3_none;
    }
  }

/* adds - axpy = 1*x+n
*/
  u3_noun
  u3qf_la_adds_real(u3_noun x_data,
                    u3_noun shape,
                    u3_noun bloq,
                    u3_noun n)
  {
    //  Unpack the data as a byte array.  We assume total length < 2**64.
    c3_d len_x = _get_length(shape);
    fprintf(stderr, "len_x: %d 0x%x units\r\n", len_x, len_x);
    c3_d siz_x = len_x * pow(2, bloq - 3);
    fprintf(stderr, "siz_x: %d 0x%x bytes\r\n", siz_x, siz_x);
    c3_y* x_bytes = (c3_y*)u3a_malloc(siz_x*sizeof(c3_y));
    u3r_bytes(0, siz_x, x_bytes, x_data);
    c3_y* y_bytes = (c3_y*)u3a_malloc((siz_x+1)*sizeof(c3_y));

    float16_t n16;
    float32_t n32;
    float64_t n64;
    float128_t n128;

    u3_noun r_data;

    //  Switch on the block size.  We assume that n fits in the target block size; Hoon typecheck should prevent.
    switch (bloq) {
      case 4:
        u3r_bytes(0, 2, (c3_y*)&n16, n);
        // set y to [n]
        for (c3_d i = 0; i < len_x; i++) {
          ((float16_t*)y_bytes)[i] = n16;
        }
        y_bytes[siz_x] = 1;  // pin head
        haxpy(len_x, (float16_t){SB_REAL16_ONE}, (float16_t*)x_bytes, 1, (float16_t*)y_bytes, 1);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), y_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return r_data;

      case 5:
        u3r_bytes(0, 4, (c3_y*)&n32, n);
        // set y to [n]
        for (c3_d i = 0; i < len_x; i++) {
          ((float32_t*)y_bytes)[i] = n32;
        }
        y_bytes[siz_x] = 1;  // pin head
        saxpy(len_x, (float32_t){SB_REAL32_ONE}, (float32_t*)x_bytes, 1, (float32_t*)y_bytes, 1);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), y_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return r_data;

      case 6:
        u3r_bytes(0, 8, (c3_y*)&n64, n);
        // set y to [n]
        for (c3_d i = 0; i < len_x; i++) {
          ((float64_t*)y_bytes)[i] = n64;
        }
        y_bytes[siz_x] = 1;  // pin head
        daxpy(len_x, (float64_t){SB_REAL64_ONE}, (float64_t*)x_bytes, 1, (float64_t*)y_bytes, 1);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), y_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return r_data;

      case 7:
        u3r_bytes(0, 16, (c3_y*)&n128, n);
        // set y to [n]
        for (c3_d i = 0; i < len_x; i++) {
          ((float128_t*)y_bytes)[i] = (float128_t){n128.v[0], n128.v[1]};
        }
        y_bytes[siz_x] = 1;  // pin head
        qaxpy(len_x, (float128_t){SB_REAL128L_ONE,SB_REAL128U_ONE}, (float128_t*)x_bytes, 1, (float128_t*)y_bytes, 1);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), y_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return r_data;

      default:
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return u3_none;
    }
  }

/* subs - axpy = -1*n+x
*/
  u3_noun
  u3qf_la_subs_real(u3_noun x_data,
                    u3_noun shape,
                    u3_noun bloq,
                    u3_noun n)
  {
    //  Unpack the data as a byte array.  We assume total length < 2**64.
    c3_d len_x = _get_length(shape);
    c3_d siz_x = len_x * pow(2, bloq - 3);
    c3_y* x_bytes = (c3_y*)u3a_malloc((siz_x+1)*sizeof(c3_y));
    u3r_bytes(0, siz_x+1, x_bytes, x_data);
    c3_y* y_bytes = (c3_y*)u3a_malloc(siz_x*sizeof(c3_y));

    float16_t n16;
    float32_t n32;
    float64_t n64;
    float128_t n128;

    u3_noun r_data;

    //  Switch on the block size.  We assume that n fits in the target block size; Hoon typecheck should prevent.
    switch (bloq) {
      case 4:
        u3r_bytes(0, 2, (c3_y*)&n16, n);
        // set y to [n]
        for (c3_d i = 0; i < len_x; i++) {
          ((float16_t*)y_bytes)[i] = n16;
        }
        haxpy(len_x, (float16_t){SB_REAL16_NEGONE}, (float16_t*)y_bytes, 1, (float16_t*)x_bytes, 1);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), x_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return r_data;

      case 5:
        u3r_bytes(0, 4, (c3_y*)&n32, n);
        // set y to [n]
        for (c3_d i = 0; i < len_x; i++) {
          ((float32_t*)y_bytes)[i] = n32;
        }
        saxpy(len_x, (float32_t){SB_REAL32_NEGONE}, (float32_t*)y_bytes, 1, (float32_t*)x_bytes, 1);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), x_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return r_data;

      case 6:
        u3r_bytes(0, 8, (c3_y*)&n64, n);
        // set y to [n]
        for (c3_d i = 0; i < len_x; i++) {
          ((float64_t*)y_bytes)[i] = n64;
        }
        daxpy(len_x, (float64_t){SB_REAL64_NEGONE}, (float64_t*)y_bytes, 1, (float64_t*)x_bytes, 1);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), x_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return r_data;

      case 7:
        u3r_bytes(0, 16, (c3_y*)&n128, n);
        // set y to [n]
        for (c3_d i = 0; i < len_x; i++) {
          ((float128_t*)y_bytes)[i] = (float128_t){n128.v[0], n128.v[1]};
        }
        qaxpy(len_x, (float128_t){SB_REAL128L_NEGONE,SB_REAL128U_NEGONE}, (float128_t*)y_bytes, 1, (float128_t*)x_bytes, 1);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), x_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return r_data;

      default:
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return u3_none;
    }
  }

/* muls - ?scal n * x
   elementwise multiplication
*/
  u3_noun
  u3qf_la_muls_real(u3_noun x_data,
                    u3_noun shape,
                    u3_noun bloq,
                    u3_noun n)
  {
    //  Unpack the data as a byte array.  We assume total length < 2**64.
    c3_d len_x = _get_length(shape);
    c3_d siz_x = len_x * pow(2, bloq - 3);
    c3_y* x_bytes = (c3_y*)u3a_malloc((siz_x+1)*sizeof(c3_y));
    u3r_bytes(0, siz_x, x_bytes, x_data);
    x_bytes[siz_x] = 1;  // pin head

    float16_t n16;
    float32_t n32;
    float64_t n64;
    float128_t n128;

    u3_noun r_data;

    //  Switch on the block size.
    switch (bloq) {
      case 4:
        u3r_bytes(0, 2, (c3_y*)&n16, n);
        hscal(len_x, n16, (float16_t*)x_bytes, 1);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), x_bytes);

        //  Clean up.
        u3a_free(x_bytes);

        return r_data;

      case 5:
        u3r_bytes(0, 4, (c3_y*)&n32, n);
        sscal(len_x, n32, (float32_t*)x_bytes, 1);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), x_bytes);

        //  Clean up.
        u3a_free(x_bytes);

        return r_data;

      case 6:
        u3r_bytes(0, 8, (c3_y*)&n64, n);
        dscal(len_x, n64, (float64_t*)x_bytes, 1);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), x_bytes);

        //  Clean up.
        u3a_free(x_bytes);

        return r_data;

      case 7:
        u3r_bytes(0, 16, (c3_y*)&(n128.v[0]), n);
        qscal(len_x, n128, (float128_t*)x_bytes, 1);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), x_bytes);

        //  Clean up.
        u3a_free(x_bytes);

        return r_data;

      default:
        u3a_free(x_bytes);

        return u3_none;
    }
  }

/* divs - ?scal 1/n * x
   elementwise division
*/
  u3_noun
  u3qf_la_divs_real(u3_noun x_data,
                    u3_noun shape,
                    u3_noun bloq,
                    u3_noun n)
  {
    //  Unpack the data as a byte array.  We assume total length < 2**64.
    c3_d len_x = _get_length(shape);
    c3_d siz_x = len_x * pow(2, bloq - 3);
    c3_y* x_bytes = (c3_y*)u3a_malloc((siz_x+1)*sizeof(c3_y));
    u3r_bytes(0, siz_x, x_bytes, x_data);
    x_bytes[siz_x] = 1;  // pin head

    float16_t n16;
    float32_t n32;
    float64_t n64;
    float128_t n128;

    u3_noun r_data;

    //  Switch on the block size.
    switch (bloq) {
      case 4:
        u3r_bytes(0, 2, (c3_y*)&(n16.v), n);
        n16 = f16_div((float16_t){SB_REAL16_ONE}, n16);
        hscal(len_x, n16, (float16_t*)x_bytes, 1);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), x_bytes);

        //  Clean up.
        u3a_free(x_bytes);

        return r_data;

      case 5:
        u3r_bytes(0, 4, (c3_y*)&(n32.v), n);
        n32 = f32_div((float32_t){SB_REAL32_ONE}, n32);
        sscal(len_x, n32, (float32_t*)x_bytes, 1);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), x_bytes);

        //  Clean up.
        u3a_free(x_bytes);

        return r_data;

      case 6:
        u3r_bytes(0, 8, (c3_y*)&(n64.v), n);
        n64 = f64_div((float64_t){SB_REAL64_ONE}, n64);
        dscal(len_x, n64, (float64_t*)x_bytes, 1);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), x_bytes);

        //  Clean up.
        u3a_free(x_bytes);

        return r_data;

      case 7:
        u3r_bytes(0, 16, (c3_y*)&(n128.v[0]), n);
        f128M_div(&((float128_t){SB_REAL128L_ONE,SB_REAL128U_ONE}), &n128, &n128);
        qscal(len_x, n128, (float128_t*)x_bytes, 1);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), x_bytes);

        //  Clean up.
        u3a_free(x_bytes);

        return r_data;

      default:
        u3a_free(x_bytes);

        return u3_none;
    }
  }

/* dot - ?dot = x · y
*/
  u3_noun
  u3qf_la_dot_real(u3_noun x_data,
                   u3_noun y_data,
                   u3_noun shape,
                   u3_noun bloq)
  {
    //  Unpack the data as a byte array.  We assume total length < 2**64.
    c3_d len_x = _get_length(shape);
    c3_d siz_x = len_x * pow(2, bloq - 3);
    c3_y* x_bytes = (c3_y*)u3a_malloc(siz_x*sizeof(c3_y));
    u3r_bytes(0, siz_x, x_bytes, x_data);
    c3_y* y_bytes = (c3_y*)u3a_malloc((siz_x+1)*sizeof(c3_y));
    u3r_bytes(0, siz_x+1, y_bytes, y_data);

    u3_noun r_data;

    //  Switch on the block size.
    switch (bloq) {
      case 4:
        hdot(len_x, (float16_t*)x_bytes, 1, (float16_t*)y_bytes, 1);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), y_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return r_data;

      case 5:
        sdot(len_x, (float32_t*)x_bytes, 1, (float32_t*)y_bytes, 1);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), y_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return r_data;

      case 6:
        ddot(len_x, (float64_t*)x_bytes, 1, (float64_t*)y_bytes, 1);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), y_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return r_data;

      case 7:
        qdot(len_x, (float128_t*)x_bytes, 1, (float128_t*)y_bytes, 1);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes((siz_x+1)*sizeof(c3_y), y_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return r_data;

      default:
        u3a_free(x_bytes);
        u3a_free(y_bytes);

        return u3_none;
    }
  }

/* diag - diag(x)
*/
  u3_noun
  u3qf_la_diag(u3_noun x_data,
               u3_noun shape,
               u3_noun bloq)
  {
    //  Assert length of dims is 2.
    assert(u3qb_lent(shape) == 2);
    //  Unpack shape into an array of dimensions.
    c3_d *dims = _get_dims(shape);
    assert(dims[0] == dims[1]);

    //  Unpack the data as a byte array.  We assume total length < 2**64.
    c3_d len_x = _get_length(shape);
    c3_d siz_x = len_x * pow(2, bloq - 3);
    c3_d stride = dims[0] * pow(2, bloq - 3);
    c3_y* x_bytes = (c3_y*)u3a_malloc((siz_x+1)*sizeof(c3_y));
    u3r_bytes(0, siz_x+1, x_bytes, x_data);
    c3_d siz_y = stride * dims[1];
    c3_y* y_bytes = (c3_y*)u3a_malloc((siz_y+1)*sizeof(c3_y));

    u3_noun r_data;

    for (c3_d i = 0; i < dims[1]; i++) {
      for (c3_d j = 0; j < stride; j++) {
        fprintf(stderr, "i*s+j = %d*%d+%d = %d // x_bytes[i]: %lx\r\n", i, stride, j, i*stride+j, x_bytes[i*stride+j + i]);
        y_bytes[i*stride+j] = x_bytes[i*stride+j + i];
      }
    }
    y_bytes[siz_y] = 1;  // pin head

    //  Unpack the result back into a noun.
    r_data = u3i_bytes((siz_y+1)*sizeof(c3_y), y_bytes);
    
    u3a_free(x_bytes);
    u3a_free(y_bytes);
    u3a_free(dims);

    return r_data;
  }

/* trace - tr(x)
*/
  u3_noun
  u3qf_la_trace_real(u3_noun x_data,
                     u3_noun shape,
                     u3_noun bloq)
  {
    u3_noun diag_data = u3qf_la_diag(x_data, shape, bloq);
    c3_d len_x0 = _get_dims(shape)[0];
    return u3qf_la_dot_real(diag_data, diag_data, u3nt(len_x0, 0x1, u3_nul), bloq);
  }

/* mmul
*/
  u3_noun
  u3qf_la_mmul_real(u3_noun x_data,
                    u3_noun y_data,
                    u3_noun x_shape,
                    u3_noun y_shape,
                    u3_noun bloq)
  {
    //  Unpack the data as a byte array.  We assume total length < 2**64.
    c3_d M = u3h(x_shape);
    c3_d Na = u3h(u3t(x_shape));
    c3_d Nb = u3h(y_shape);
    c3_d P = u3h(u3t(y_shape));

    assert(u3_nul == u3t(u3t(x_shape)));
    assert(Na == Nb);
    c3_d N = Na;
    assert(u3_nul == u3t(u3t(y_shape)));

    c3_y* x_bytes = (c3_y*)u3a_malloc((M*N)*sizeof(c3_y));
    u3r_bytes(0, M*N, x_bytes, x_data);
    c3_y* y_bytes = (c3_y*)u3a_malloc((N*P)*sizeof(c3_y));
    u3r_bytes(0, N*P, y_bytes, y_data);
    c3_y* c_bytes = (c3_y*)u3a_malloc((M*P)*sizeof(c3_y));

    u3_noun r_data;

    //  Switch on the block size.
    switch (bloq) {
      case 4:
        hgemm('N', 'N', M, N, P, (float16_t){SB_REAL16_ONE}, (float16_t*)x_bytes, N, (float16_t*)y_bytes, N, (float16_t){SB_REAL16_ZERO}, (float16_t*)c_bytes, P);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes(M*P, c_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);
        u3a_free(c_bytes);

        return u3nc(u3nq(u3nl(M, P, u3_none), bloq, c3__real, u3_nul), r_data);

      case 5:
        sgemm('N', 'N', M, N, P, (float32_t){SB_REAL32_ONE}, (float32_t*)x_bytes, N, (float32_t*)y_bytes, N, (float32_t){SB_REAL32_ZERO}, (float32_t*)c_bytes, P);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes(M*P, c_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);
        u3a_free(c_bytes);

        return u3nc(u3nq(u3nl(M, P, u3_none), bloq, c3__real, u3_nul), r_data);

      case 6:
        dgemm('N', 'N', M, N, P, (float64_t){SB_REAL64_ONE}, (float64_t*)x_bytes, N, (float64_t*)y_bytes, N, (float64_t){SB_REAL64_ZERO}, (float64_t*)c_bytes, P);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes(M*P, c_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);
        u3a_free(c_bytes);

        return u3nc(u3nq(u3nl(M, P, u3_none), bloq, c3__real, u3_nul), r_data);

      case 7:
        qgemm('N', 'N', M, N, P, (float128_t){SB_REAL128L_ONE,SB_REAL128U_ONE}, (float128_t*)x_bytes, N, (float128_t*)y_bytes, N, (float128_t){SB_REAL128L_ZERO,SB_REAL128U_ZERO}, (float128_t*)c_bytes, P);

        //  Unpack the result back into a noun.
        r_data = u3i_bytes(M*P, c_bytes);

        //  Clean up.
        u3a_free(x_bytes);
        u3a_free(y_bytes);
        u3a_free(c_bytes);

        return u3nc(u3nq(u3nl(M, P, u3_none), bloq, c3__real, u3_nul), r_data);

      default:
        u3a_free(x_bytes);
        u3a_free(y_bytes);
        u3a_free(c_bytes);
        
        return u3_none;
    }
  }

  u3_noun
  u3wf_la_add(u3_noun cor)
  {
    // Each argument is a ray, [=meta data=@ux]
    u3_noun x_meta, x_data,
            y_meta, y_data;

    if ( c3n == u3r_mean(cor,
                         u3x_sam_4, &x_meta,
                         u3x_sam_5, &x_data,
                         u3x_sam_6, &y_meta,
                         u3x_sam_7, &y_data,
                         0) ||
         c3n == u3ud(x_data) ||
         c3n == u3ud(y_data) )
    {
      return u3m_bail(c3__exit);
    } else {
      u3_noun x_shape, x_bloq, x_kind, x_fxp,
              y_shape, y_bloq, y_kind, y_fxp,
              rnd;
      if ( c3n == u3r_mean(x_meta,
                            2, &x_shape,
                            6, &x_bloq,
                           14, &x_kind,
                           15, &x_fxp,
                            0) ||
           c3n == u3r_mean(y_meta,
                            2, &y_shape,
                            6, &y_bloq,
                           14, &y_kind,
                           15, &y_fxp,
                            0) ||
           c3n == u3ud(x_bloq) ||
           c3n == u3ud(y_bloq) ||
           c3n == u3ud(x_kind) ||
           c3n == u3ud(y_kind) ||
           c3n == u3r_sing(x_shape, y_shape) ||
           c3n == u3r_sing(x_bloq, y_bloq) ||
           c3n == u3r_sing(x_kind, y_kind) ||
           //  fxp does not need to match here so no check
           c3n == u3r_mean(cor, u3x_con_sam, &rnd, 0)
         )
      {
        return u3m_bail(c3__exit);
      } else {
        switch (x_kind) {
          case c3__real:
            fprintf(stderr, "\r\n>>  add_real (rnd:  0x%lx)\r\n", rnd);
            fprintf(stderr, ">>  length: %lx\r\n", _get_length(x_shape));
            _set_rounding(rnd);
            // u3_noun r_data = u3qf_la_add_real(x_data, y_data, x_shape, x_bloq);
            fprintf(stderr, ">>  x_shape: [%lx %lx %lx]\r\n", u3h(x_shape), u3h(u3t(x_shape)), u3t(u3t(x_shape)));
            fprintf(stderr, ">>  x_bloq: 0x%lx\r\n", x_bloq);
            fprintf(stderr, ">>  x_kind: 0x%lx\r\n", x_kind);
            // return u3nc(u3nq(x_shape, x_bloq, x_kind, x_fxp), y_data);
            return u3qf_la_add_real(x_data, y_data, x_shape, x_bloq, x_kind, x_fxp);
            break;

          default:
            return u3_none;
        }
      }
    }
  }

  u3_noun
  u3wf_la_sub(u3_noun cor)
  {
    // Each argument is a ray, [=meta data=@ux]
    u3_noun x_meta, x_data,
            y_meta, y_data;

    if ( c3n == u3r_mean(cor,
                         u3x_sam_4, &x_meta,
                         u3x_sam_5, &x_data,
                         u3x_sam_6, &y_meta,
                         u3x_sam_7, &y_data,
                         0) ||
         c3n == u3ud(x_data) ||
         c3n == u3ud(y_data) )
    {
      return u3m_bail(c3__exit);
    } else {
      u3_noun x_shape, x_bloq, x_kind, x_fxp,
              y_shape, y_bloq, y_kind, y_fxp,
              rnd;
      if ( c3n == u3r_mean(x_meta,
                            2, &x_shape,
                            6, &x_bloq,
                           14, &x_kind,
                           15, &x_fxp,
                            0) ||
           c3n == u3r_mean(y_meta,
                            2, &y_shape,
                            6, &y_bloq,
                           14, &y_kind,
                           15, &y_fxp,
                            0) ||
           c3n == u3r_sing(x_shape, y_shape) ||
           c3n == u3r_sing(x_bloq, y_bloq) ||
           c3n == u3r_sing(x_kind, y_kind) ||
           //  fxp does not need to match so no check
           c3n == u3r_mean(cor, u3x_con_sam, &rnd, 0)
         )
      {
        return u3m_bail(c3__exit);
      } else {
        switch (x_kind) {
          case c3__real:
            _set_rounding(rnd);
            u3_noun r_data = u3qf_la_sub_real(x_data, y_data, x_shape, x_bloq);
            return u3nc(u3nq(x_shape, x_bloq, x_kind, x_fxp), r_data);
            break;

          default:
            return u3_none;
        }
      }
    }
  }

  u3_noun
  u3wf_la_mul(u3_noun cor)
  {
    // Each argument is a ray, [=meta data=@ux]
    u3_noun x_meta, x_data,
            y_meta, y_data;

    if ( c3n == u3r_mean(cor,
                         u3x_sam_4, &x_meta,
                         u3x_sam_5, &x_data,
                         u3x_sam_6, &y_meta,
                         u3x_sam_7, &y_data,
                         0) ||
         c3n == u3ud(x_data) ||
         c3n == u3ud(y_data) )
    {
      return u3m_bail(c3__exit);
    } else {
      u3_noun x_shape, x_bloq, x_kind, x_fxp,
              y_shape, y_bloq, y_kind, y_fxp,
              rnd;
      if ( c3n == u3r_mean(x_meta,
                            2, &x_shape,
                            6, &x_bloq,
                           14, &x_kind,
                           15, &x_fxp,
                            0) ||
           c3n == u3r_mean(y_meta,
                            2, &y_shape,
                            6, &y_bloq,
                           14, &y_kind,
                           15, &y_fxp,
                            0) ||
           c3n == u3r_sing(x_shape, y_shape) ||
           c3n == u3r_sing(x_bloq, y_bloq) ||
           c3n == u3r_sing(x_kind, y_kind) ||
           //  fxp does not need to match so no check
           c3n == u3r_mean(cor, u3x_con_sam, &rnd, 0)
         )
      {
        return u3m_bail(c3__exit);
      } else {
        switch (x_kind) {
          case c3__real:
            _set_rounding(rnd);
            u3_noun r_data = u3qf_la_mul_real(x_data, y_data, x_shape, x_bloq);
            return u3nc(u3nq(x_shape, x_bloq, x_kind, x_fxp), r_data);
            break;

          default:
            return u3_none;
        }
      }
    }
  }

  u3_noun
  u3wf_la_div(u3_noun cor)
  {
    // Each argument is a ray, [=meta data=@ux]
    u3_noun x_meta, x_data,
            y_meta, y_data;

    if ( c3n == u3r_mean(cor,
                         u3x_sam_4, &x_meta,
                         u3x_sam_5, &x_data,
                         u3x_sam_6, &y_meta,
                         u3x_sam_7, &y_data,
                         0) ||
         c3n == u3ud(x_data) ||
         c3n == u3ud(y_data) )
    {
      return u3m_bail(c3__exit);
    } else {
      u3_noun x_shape, x_bloq, x_kind, x_fxp,
              y_shape, y_bloq, y_kind, y_fxp,
              rnd;
      if ( c3n == u3r_mean(x_meta,
                            2, &x_shape,
                            6, &x_bloq,
                           14, &x_kind,
                           15, &x_fxp,
                            0) ||
           c3n == u3r_mean(y_meta,
                            2, &y_shape,
                            6, &y_bloq,
                           14, &y_kind,
                           15, &y_fxp,
                            0) ||
           c3n == u3r_sing(x_shape, y_shape) ||
           c3n == u3r_sing(x_bloq, y_bloq) ||
           c3n == u3r_sing(x_kind, y_kind) ||
           //  fxp does not need to match so no check
           c3n == u3r_mean(cor, u3x_con_sam, &rnd, 0)
         )
      {
        return u3m_bail(c3__exit);
      } else {
        switch (x_kind) {
          case c3__real:
            _set_rounding(rnd);
            u3_noun r_data = u3qf_la_div_real(x_data, y_data, x_shape, x_bloq);
            // fprintf(stderr, "x_shape[0]: %x\r\n", u3h(x_shape));
            // fprintf(stderr, "x_shape[1]: %x\r\n", u3h(u3t(x_shape)));
            // fprintf(stderr, "x_shape[2]: %x\r\n", u3t(u3t(x_shape)));
            return u3nc(u3nq(x_shape, x_bloq, x_kind, x_fxp), r_data);
            break;

          default:
            return u3_none;
        }
      }
    }
  }

  u3_noun
  u3wf_la_adds(u3_noun cor)
  {
    // Each argument is a ray, [=meta data=@ux]
    u3_noun x_meta, x_data, n;

    if ( c3n == u3r_mean(cor,
                         u3x_sam_4, &x_meta,
                         u3x_sam_5, &x_data,
                         u3x_sam_3, &n,
                         0) ||
         c3n == u3ud(x_data) ||
         c3n == u3ud(n) )
    {
      return u3m_bail(c3__exit);
    } else {
      u3_noun x_shape, x_bloq, x_kind, x_fxp,
              rnd;
      if ( c3n == u3r_mean(x_meta,
                            2, &x_shape,
                            6, &x_bloq,
                           14, &x_kind,
                           15, &x_fxp,
                            0) ||
           //  shape does not matter so no check
           //  bloq does not matter so no check
           //  kind does not matter so no check
           //  fxp does not need to match so no check
           c3n == u3r_mean(cor, u3x_con_sam, &rnd, 0)
         )
      {
        return u3m_bail(c3__exit);
      } else {
        switch (x_kind) {
          case c3__real:
            _set_rounding(rnd);
            u3_noun r_data = u3qf_la_adds_real(x_data, x_shape, x_bloq, n);
            return u3nc(u3nq(x_shape, x_bloq, x_kind, x_fxp), r_data);
            break;

          default:
            return u3_none;
        }
      }
    }
  }

  u3_noun
  u3wf_la_subs(u3_noun cor)
  {
    // Each argument is a ray, [=meta data=@ux]
    u3_noun x_meta, x_data, n;

    if ( c3n == u3r_mean(cor,
                         u3x_sam_4, &x_meta,
                         u3x_sam_5, &x_data,
                         u3x_sam_3, &n,
                         0) ||
         c3n == u3ud(x_data) ||
         c3n == u3ud(n) )
    {
      return u3m_bail(c3__exit);
    } else {
      u3_noun x_shape, x_bloq, x_kind, x_fxp,
              rnd;
      if ( c3n == u3r_mean(x_meta,
                            2, &x_shape,
                            6, &x_bloq,
                           14, &x_kind,
                           15, &x_fxp,
                            0) ||
           //  shape does not matter so no check
           //  bloq does not matter so no check
           //  kind does not matter so no check
           //  fxp does not need to match so no check
           c3n == u3r_mean(cor, u3x_con_sam, &rnd, 0)
         )
      {
        return u3m_bail(c3__exit);
      } else {
        switch (x_kind) {
          case c3__real:
            _set_rounding(rnd);
            u3_noun r_data = u3qf_la_subs_real(x_data, x_shape, x_bloq, n);
            return u3nc(u3nq(x_shape, x_bloq, x_kind, x_fxp), r_data);
            break;

          default:
            return u3_none;
        }
      }
    }
  }

  u3_noun
  u3wf_la_muls(u3_noun cor)
  {
    // Each argument is a ray, [=meta data=@ux]
    u3_noun x_meta, x_data, n;

    if ( c3n == u3r_mean(cor,
                         u3x_sam_4, &x_meta,
                         u3x_sam_5, &x_data,
                         u3x_sam_3, &n,
                         0) ||
         c3n == u3ud(x_data) ||
         c3n == u3ud(n) )
    {
      return u3m_bail(c3__exit);
    } else {
      u3_noun x_shape, x_bloq, x_kind, x_fxp,
              rnd;
      if ( c3n == u3r_mean(x_meta,
                            2, &x_shape,
                            6, &x_bloq,
                           14, &x_kind,
                           15, &x_fxp,
                            0) ||
           //  shape does not matter so no check
           //  bloq does not matter so no check
           //  kind does not matter so no check
           //  fxp does not need to match so no check
           c3n == u3r_mean(cor, u3x_con_sam, &rnd, 0)
         )
      {
        return u3m_bail(c3__exit);
      } else {
        switch (x_kind) {
          case c3__real:
            _set_rounding(rnd);
            u3_noun r_data = u3qf_la_muls_real(x_data, x_shape, x_bloq, n);
            return u3nc(u3nq(x_shape, x_bloq, x_kind, x_fxp), r_data);
            break;

          default:
            return u3_none;
        }
      }
    }
  }

  u3_noun
  u3wf_la_divs(u3_noun cor)
  {
    // Each argument is a ray, [=meta data=@ux]
    u3_noun x_meta, x_data, n;

    if ( c3n == u3r_mean(cor,
                         u3x_sam_4, &x_meta,
                         u3x_sam_5, &x_data,
                         u3x_sam_3, &n,
                         0) ||
         c3n == u3ud(x_data) ||
         c3n == u3ud(n) )
    {
      return u3m_bail(c3__exit);
    } else {
      u3_noun x_shape, x_bloq, x_kind, x_fxp,
              rnd;
      if ( c3n == u3r_mean(x_meta,
                            2, &x_shape,
                            6, &x_bloq,
                           14, &x_kind,
                           15, &x_fxp,
                            0) ||
           //  shape does not matter so no check
           //  bloq does not matter so no check
           //  kind does not matter so no check
           //  fxp does not need to match so no check
           c3n == u3r_mean(cor, u3x_con_sam, &rnd, 0)
         )
      {
        return u3m_bail(c3__exit);
      } else {
        switch (x_kind) {
          case c3__real:
            _set_rounding(rnd);
            u3_noun r_data = u3qf_la_divs_real(x_data, x_shape, x_bloq, n);
            return u3nc(u3nq(x_shape, x_bloq, x_kind, x_fxp), r_data);
            break;

          default:
            return u3_none;
        }
      }
    }
  }

  u3_noun
  u3wf_la_dot(u3_noun cor)
  {
    // Each argument is a ray, [=meta data=@ux]
    u3_noun x_meta, x_data,
            y_meta, y_data;

    if ( c3n == u3r_mean(cor,
                         u3x_sam_4, &x_meta,
                         u3x_sam_5, &x_data,
                         u3x_sam_6, &y_meta,
                         u3x_sam_7, &y_data,
                         0) ||
         c3n == u3ud(x_data) ||
         c3n == u3ud(y_data) )
    {
      return u3m_bail(c3__exit);
    } else {
      u3_noun x_shape, x_bloq, x_kind, x_fxp,
              y_shape, y_bloq, y_kind, y_fxp,
              rnd;
      if ( c3n == u3r_mean(x_meta,
                            2, &x_shape,
                            6, &x_bloq,
                           14, &x_kind,
                           15, &x_fxp,
                            0) ||
           c3n == u3r_mean(y_meta,
                            2, &y_shape,
                            6, &y_bloq,
                           14, &y_kind,
                           15, &y_fxp,
                            0) ||
           c3n == u3r_sing(x_shape, y_shape) ||
           c3n == u3r_sing(x_bloq, y_bloq) ||
           c3n == u3r_sing(x_kind, y_kind) ||
           c3n == u3r_sing(x_fxp, y_fxp) ||
           c3n == u3r_mean(cor, u3x_con_sam, &rnd, 0)
         )
      {
        return u3m_bail(c3__exit);
      } else {
        switch (x_kind) {
          case c3__real:
            _set_rounding(rnd);
            u3_noun r_data = u3qf_la_dot_real(x_data, y_data, x_shape, x_bloq);
            return u3nc(u3nq(x_shape, x_bloq, x_kind, x_fxp), r_data);
            break;

          default:
            return u3_none;
        }
      }
    }
  }

  u3_noun
  u3wf_la_diag(u3_noun cor)
  {
    // Each argument is a ray, [=meta data=@ux]
    u3_noun x_meta, x_data;

    if ( c3n == u3r_mean(cor,
                         u3x_sam_2, &x_meta,
                         u3x_sam_3, &x_data,
                         0) ||
         c3n == u3ud(x_data) )
    {
      return u3m_bail(c3__exit);
    } else {
      u3_noun x_shape, x_bloq, x_kind, x_fxp,
              rnd;
      if ( c3n == u3r_mean(x_meta,
                            2, &x_shape,
                            6, &x_bloq,
                           14, &x_kind,
                           15, &x_fxp,
                            0)
          //  c3n == u3r_sing(x_shape, y_shape) ||
          //  c3n == u3r_sing(x_bloq, y_bloq) ||
          //  c3n == u3r_sing(x_kind, y_kind) ||
          //  c3n == u3r_sing(x_fxp, y_fxp) ||
          //  c3n == u3r_mean(cor, u3x_con_sam, &rnd, 0)
         )
      {
        return u3m_bail(c3__exit);
      } else {
        u3_noun r_data = u3qf_la_diag(x_data, x_shape, x_bloq);
        c3_d len_x0 = _get_dims(x_shape)[0];
        return u3nc(u3nq(u3nt(len_x0, 0x1, u3_nul), x_bloq, x_kind, x_fxp), r_data);
      }
    }
  }

  u3_noun
  u3wf_la_trace(u3_noun cor)
  {
    // Each argument is a ray, [=meta data=@ux]
    u3_noun x_meta, x_data;

    if ( c3n == u3r_mean(cor,
                         u3x_sam_2, &x_meta,
                         u3x_sam_3, &x_data,
                         0) ||
         c3n == u3ud(x_data) )
    {
      return u3m_bail(c3__exit);
    } else {
      u3_noun x_shape, x_bloq, x_kind, x_fxp,
              rnd;
      if ( c3n == u3r_mean(x_meta,
                            2, &x_shape,
                            6, &x_bloq,
                           14, &x_kind,
                           15, &x_fxp,
                            0)
          //  c3n == u3r_sing(x_shape, y_shape) ||
          //  c3n == u3r_sing(x_bloq, y_bloq) ||
          //  c3n == u3r_sing(x_kind, y_kind) ||
          //  c3n == u3r_sing(x_fxp, y_fxp) ||
          //  c3n == u3r_mean(cor, u3x_con_sam, &rnd, 0)
         )
      {
        return u3m_bail(c3__exit);
      } else {
        switch (x_kind) {
          case c3__real:
            _set_rounding(rnd);
            u3_noun r_data = u3qf_la_trace_real(x_data, x_shape, x_bloq);
            return u3nc(u3nq(u3nt(0x1, 0x1, u3_nul), x_bloq, x_kind, x_fxp), r_data);
            break;

          default:
            return u3_none;
        }
      }
    }
  }

  u3_noun
  u3wf_la_mmul(u3_noun cor)
  {
    // Each argument is a ray, [=meta data=@ux]
    u3_noun x_meta, x_data,
            y_meta, y_data;

    if ( c3n == u3r_mean(cor,
                         u3x_sam_4, &x_meta,
                         u3x_sam_5, &x_data,
                         u3x_sam_6, &y_meta,
                         u3x_sam_7, &y_data,
                         0) ||
         c3n == u3ud(x_data) ||
         c3n == u3ud(y_data) )
    {
      return u3m_bail(c3__exit);
    } else {
      u3_noun x_shape, x_bloq, x_kind, x_fxp,
              y_shape, y_bloq, y_kind, y_fxp,
              rnd;
      if ( c3n == u3r_mean(x_meta,
                            2, &x_shape,
                            6, &x_bloq,
                           14, &x_kind,
                           15, &x_fxp,
                            0) ||
           c3n == u3r_mean(y_meta,
                            2, &y_shape,
                            6, &y_bloq,
                           14, &y_kind,
                           15, &y_fxp,
                            0) ||
           c3n == u3r_sing(x_bloq, y_bloq) ||
           c3n == u3r_sing(x_kind, y_kind) ||
           //  fxp does not need to match so no check
           c3n == u3r_mean(cor, u3x_con_sam, &rnd, 0)
         )
      {
        return u3m_bail(c3__exit);
      } else {
        switch (x_kind) {
          case c3__real:
            _set_rounding(rnd);
            return u3qf_la_mmul_real(x_data, y_data, x_shape, y_shape, x_bloq);
            break;

          default:
            return u3_none;
        }
      }
    }
  }
