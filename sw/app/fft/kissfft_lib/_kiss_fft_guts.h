/*
 *  Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
 *  This file is part of KISS FFT - https://github.com/mborgerding/kissfft
 *
 *  SPDX-License-Identifier: BSD-3-Clause
 *  See COPYING file for more information.
 */

/* kiss_fft.h
   defines kiss_fft_scalar as either short or a float type
   and defines
   typedef struct { kiss_fft_scalar r; kiss_fft_scalar i; }kiss_fft_cpx; */

#ifndef _kiss_fft_guts_h
#define _kiss_fft_guts_h

#include "kiss_fft.h"
#include "kiss_fft_log.h"
#include <limits.h>

#define MAXFACTORS 32
/* e.g. an fft of length 128 has 4 factors
 as far as kissfft is concerned
 4*4*4*2
 */

struct kiss_fft_state{
    int nfft;
    int inverse;

    // THE C_EXP OPERATIONS ARE PRECALCULATED
    int *factors;
    kiss_fft_cpx *twiddles;

//     int factors[2*MAXFACTORS];
//     kiss_fft_cpx twiddles[1];
};

/*
  Explanation of macros dealing with complex math:

   C_MUL(m,a,b)         : m = a*b
   C_FIXDIV( c , div )  : if a fixed point impl., c /= div. noop otherwise
   C_SUB( res, a,b)     : res = a - b
   C_SUBFROM( res , a)  : res -= a
   C_ADDTO( res , a)    : res += a
 * */
#ifdef FIXED_POINT
#include <stdint.h>
#if (FIXED_POINT==32)
# define FRACBITS 31
# define SAMPPROD int64_t
#define SAMP_MAX INT32_MAX
#define SAMP_MIN INT32_MIN
#else
# define FRACBITS 15
# define SAMPPROD int32_t
#define SAMP_MAX INT16_MAX
#define SAMP_MIN INT16_MIN
#endif

#if defined(CHECK_OVERFLOW)
#  define CHECK_OVERFLOW_OP(a,op,b)  \
    if ( (SAMPPROD)(a) op (SAMPPROD)(b) > SAMP_MAX || (SAMPPROD)(a) op (SAMPPROD)(b) < SAMP_MIN ) { \
        KISS_FFT_WARNING("overflow (%d " #op" %d) = %ld", (a),(b),(SAMPPROD)(a) op (SAMPPROD)(b)); }
#endif


#   define smul(a,b) ( (SAMPPROD)(a)*(b) )
#   define sround( x )  (kiss_fft_scalar)( ( (x) + (1<<(FRACBITS-1)) ) >> FRACBITS )

#   define S_MUL(a,b) sround( smul(a,b) )
#if defined(KISSFFT_USE_IM_CUSTOM) && (FIXED_POINT==16) && defined(__riscv)
static inline uint32_t kissfft_pack_cpx(kiss_fft_cpx a) {
    return ((uint16_t)a.r) | ((uint32_t)(uint16_t)a.i << 16);
}

static inline kiss_fft_cpx kissfft_unpack_cpx(uint32_t v) {
    kiss_fft_cpx r;
    r.r = (int16_t)(v & 0xFFFFu);
    r.i = (int16_t)(v >> 16);
    return r;
}

typedef uint32_t kissfft_u32_alias_t __attribute__((__may_alias__));

static inline uint32_t kissfft_load_u32(const kiss_fft_cpx *p) {
    const kissfft_u32_alias_t *pa = (const kissfft_u32_alias_t *)__builtin_assume_aligned(p, 4);
    return *pa;
}

static inline void kissfft_store_u32(kiss_fft_cpx *p, uint32_t v) {
    kissfft_u32_alias_t *pa = (kissfft_u32_alias_t *)__builtin_assume_aligned(p, 4);
    *pa = v;
}

static inline uint32_t kissfft_cadd_u32(uint32_t a, uint32_t b) {
    uint32_t r;
    __asm__ volatile (".insn r 0x0B, 0, 0x00, %0, %1, %2" : "=r"(r) : "r"(a), "r"(b));
    return r;
}

static inline uint32_t kissfft_csub_u32(uint32_t a, uint32_t b) {
    uint32_t r;
    __asm__ volatile (".insn r 0x0B, 1, 0x00, %0, %1, %2" : "=r"(r) : "r"(a), "r"(b));
    return r;
}

static inline uint32_t kissfft_cmul_u32(uint32_t a, uint32_t b) {
    uint32_t r;
    __asm__ volatile (".insn r 0x0B, 2, 0x00, %0, %1, %2" : "=r"(r) : "r"(a), "r"(b));
    return r;
}

static inline uint32_t kissfft_cmadd_u32(uint32_t acc, uint32_t a, uint32_t b) {
    __asm__ volatile (".insn r 0x0B, 3, 0x00, %0, %1, %2" : "+r"(acc) : "r"(a), "r"(b));
    return acc;
}

static inline uint32_t kissfft_cmsub_u32(uint32_t acc, uint32_t a, uint32_t b) {
    __asm__ volatile (".insn r 0x0B, 4, 0x00, %0, %1, %2" : "+r"(acc) : "r"(a), "r"(b));
    return acc;
}

static inline void kissfft_twld_u32(uint32_t idx, uint32_t val) {
    __asm__ volatile (".insn r 0x0B, 5, 0x00, x0, %0, %1" :: "r"(idx), "r"(val));
}

static inline void kissfft_twcfg_u32(uint32_t mode) {
    uint32_t z = 0;
    __asm__ volatile (".insn r 0x0B, 6, 0x00, x0, %0, %1" :: "r"(mode), "r"(z));
}

static inline uint32_t kissfft_bfy2_lo_u32(uint32_t a, uint32_t b, uint32_t tw) {
    uint32_t r;
    __asm__ volatile (".insn r4 0x2B, 0, 0, %0, %1, %2, %3" : "=r"(r) : "r"(a), "r"(b), "r"(tw));
    return r;
}

static inline uint32_t kissfft_bfy2_hi_u32(uint32_t dep) {
    uint32_t r;
    __asm__ volatile (".insn r4 0x2B, 1, 0, %0, %1, %1, %1" : "=r"(r) : "r"(dep));
    return r;
}

static inline void kissfft_bfy4_s0(uint32_t a, uint32_t b, uint32_t t1) {
    uint32_t tmp;
    __asm__ volatile (".insn r4 0x5B, 0, 0, %0, %1, %2, %3" : "=r"(tmp) : "r"(a), "r"(b), "r"(t1));
    (void)tmp;
}

static inline void kissfft_bfy4_s1(uint32_t c, uint32_t d, uint32_t t2) {
    uint32_t tmp;
    __asm__ volatile (".insn r4 0x5B, 1, 0, %0, %1, %2, %3" : "=r"(tmp) : "r"(c), "r"(d), "r"(t2));
    (void)tmp;
}

static inline void kissfft_bfy4_s2(uint32_t t3) {
    uint32_t z = 0;
    uint32_t tmp;
    __asm__ volatile (".insn r4 0x5B, 2, 0, %0, %1, %2, %3" : "=r"(tmp) : "r"(z), "r"(z), "r"(t3));
    (void)tmp;
}

static inline uint32_t kissfft_bfy4_o0(uint32_t dep) {
    uint32_t r;
    __asm__ volatile (".insn r4 0x5B, 3, 0, %0, %1, %1, %1" : "=r"(r) : "r"(dep));
    return r;
}

static inline uint32_t kissfft_bfy4_o1(uint32_t dep) {
    uint32_t r;
    __asm__ volatile (".insn r4 0x5B, 4, 0, %0, %1, %1, %1" : "=r"(r) : "r"(dep));
    return r;
}

static inline uint32_t kissfft_bfy4_o2(uint32_t dep) {
    uint32_t r;
    __asm__ volatile (".insn r4 0x5B, 5, 0, %0, %1, %1, %1" : "=r"(r) : "r"(dep));
    return r;
}

static inline uint32_t kissfft_bfy4_o3(uint32_t dep) {
    uint32_t r;
    __asm__ volatile (".insn r4 0x5B, 6, 0, %0, %1, %1, %1" : "=r"(r) : "r"(dep));
    return r;
}

#   define C_MUL(m,a,b) \
      do{ \
          uint32_t _a = kissfft_pack_cpx(a); \
          uint32_t _b = kissfft_pack_cpx(b); \
          (m) = kissfft_unpack_cpx(kissfft_cmul_u32(_a, _b)); \
      }while(0)

#   define C_ADD(res,a,b) \
      do { \
          uint32_t _a = kissfft_pack_cpx(a); \
          uint32_t _b = kissfft_pack_cpx(b); \
          (res) = kissfft_unpack_cpx(kissfft_cadd_u32(_a, _b)); \
      }while(0)

#   define C_SUB(res,a,b) \
      do { \
          uint32_t _a = kissfft_pack_cpx(a); \
          uint32_t _b = kissfft_pack_cpx(b); \
          (res) = kissfft_unpack_cpx(kissfft_csub_u32(_a, _b)); \
      }while(0)

#   define C_ADDTO(res,a) \
      do { \
          uint32_t _r = kissfft_pack_cpx(res); \
          uint32_t _a = kissfft_pack_cpx(a); \
          (res) = kissfft_unpack_cpx(kissfft_cadd_u32(_r, _a)); \
      }while(0)

#   define C_SUBFROM(res,a) \
      do { \
          uint32_t _r = kissfft_pack_cpx(res); \
          uint32_t _a = kissfft_pack_cpx(a); \
          (res) = kissfft_unpack_cpx(kissfft_csub_u32(_r, _a)); \
      }while(0)

#   define C_MULACC(res,a,b) \
      do{ \
          uint32_t _acc = kissfft_pack_cpx(res); \
          uint32_t _a = kissfft_pack_cpx(a); \
          uint32_t _b = kissfft_pack_cpx(b); \
          (res) = kissfft_unpack_cpx(kissfft_cmadd_u32(_acc, _a, _b)); \
      }while(0)

#   define C_MULSUB(res,a,b) \
      do{ \
          uint32_t _acc = kissfft_pack_cpx(res); \
          uint32_t _a = kissfft_pack_cpx(a); \
          uint32_t _b = kissfft_pack_cpx(b); \
          (res) = kissfft_unpack_cpx(kissfft_cmsub_u32(_acc, _a, _b)); \
      }while(0)

#   define C_BFY2_LO(a,b,tw) kissfft_bfy2_lo_u32((a),(b),(tw))
#   define C_BFY2_HI(dep) kissfft_bfy2_hi_u32((dep))
#   define C_BFY4_S0(a,b,t1) kissfft_bfy4_s0((a),(b),(t1))
#   define C_BFY4_S1(c,d,t2) kissfft_bfy4_s1((c),(d),(t2))
#   define C_BFY4_S2(t3) kissfft_bfy4_s2((t3))
#   define C_BFY4_O0(dep) kissfft_bfy4_o0((dep))
#   define C_BFY4_O1(dep) kissfft_bfy4_o1((dep))
#   define C_BFY4_O2(dep) kissfft_bfy4_o2((dep))
#   define C_BFY4_O3(dep) kissfft_bfy4_o3((dep))
#else
#   define C_MUL(m,a,b) \
      do{ (m).r = sround( smul((a).r,(b).r) - smul((a).i,(b).i) ); \
          (m).i = sround( smul((a).r,(b).i) + smul((a).i,(b).r) ); }while(0)
#   define C_MULACC(res,a,b) \
    do{ kiss_fft_cpx _t; C_MUL(_t, a, b); C_ADDTO(res, _t); }while(0)
#   define C_MULSUB(res,a,b) \
    do{ kiss_fft_cpx _t; C_MUL(_t, a, b); C_SUBFROM(res, _t); }while(0)
#endif

#   define DIVSCALAR(x,k) \
    (x) = sround( smul(  x, SAMP_MAX/k ) )

#   define C_FIXDIV(c,div) \
    do {    DIVSCALAR( (c).r , div);  \
        DIVSCALAR( (c).i  , div); }while (0)

#   define C_MULBYSCALAR( c, s ) \
    do{ (c).r =  sround( smul( (c).r , s ) ) ;\
        (c).i =  sround( smul( (c).i , s ) ) ; }while(0)

#else  /* not FIXED_POINT*/

#   define S_MUL(a,b) ( (a)*(b) )
#define C_MUL(m,a,b) \
    do{ (m).r = (a).r*(b).r - (a).i*(b).i;\
        (m).i = (a).r*(b).i + (a).i*(b).r; }while(0)
#   define C_FIXDIV(c,div) /* NOOP */
#   define C_MULBYSCALAR( c, s ) \
    do{ (c).r *= (s);\
        (c).i *= (s); }while(0)
#endif

#ifndef CHECK_OVERFLOW_OP
#  define CHECK_OVERFLOW_OP(a,op,b) /* noop */
#endif

#ifndef C_ADD
#define  C_ADD( res, a,b)\
    do { \
        CHECK_OVERFLOW_OP((a).r,+,(b).r)\
        CHECK_OVERFLOW_OP((a).i,+,(b).i)\
        (res).r=(a).r+(b).r;  (res).i=(a).i+(b).i; \
    }while(0)
#endif
#ifndef C_SUB
#define  C_SUB( res, a,b)\
    do { \
        CHECK_OVERFLOW_OP((a).r,-,(b).r)\
        CHECK_OVERFLOW_OP((a).i,-,(b).i)\
        (res).r=(a).r-(b).r;  (res).i=(a).i-(b).i; \
    }while(0)
#endif
#ifndef C_ADDTO
#define C_ADDTO( res , a)\
    do { \
        CHECK_OVERFLOW_OP((res).r,+,(a).r)\
        CHECK_OVERFLOW_OP((res).i,+,(a).i)\
        (res).r += (a).r;  (res).i += (a).i;\
    }while(0)
#endif

#ifndef C_SUBFROM
#define C_SUBFROM( res , a)\
    do {\
        CHECK_OVERFLOW_OP((res).r,-,(a).r)\
        CHECK_OVERFLOW_OP((res).i,-,(a).i)\
        (res).r -= (a).r;  (res).i -= (a).i; \
    }while(0)
#endif


#ifdef FIXED_POINT
#  define KISS_FFT_COS(phase)  floor(.5+SAMP_MAX * cos (phase))
#  define KISS_FFT_SIN(phase)  floor(.5+SAMP_MAX * sin (phase))
#  define HALF_OF(x) ((x)>>1)
#elif defined(USE_SIMD)
#  define KISS_FFT_COS(phase) _mm_set1_ps( cos(phase) )
#  define KISS_FFT_SIN(phase) _mm_set1_ps( sin(phase) )
#  define HALF_OF(x) ((x)*_mm_set1_ps(.5))
#else
#  define KISS_FFT_COS(phase) (kiss_fft_scalar) cos(phase)
#  define KISS_FFT_SIN(phase) (kiss_fft_scalar) sin(phase)
#  define HALF_OF(x) ((x)*((kiss_fft_scalar).5))
#endif

#define  kf_cexp(x,phase) \
    do{ \
        (x)->r = KISS_FFT_COS(phase);\
        (x)->i = KISS_FFT_SIN(phase);\
    }while(0)


/* a debugging function */
#define pcpx(c)\
    KISS_FFT_DEBUG("%g + %gi\n",(double)((c)->r),(double)((c)->i))


#ifdef KISS_FFT_USE_ALLOCA
// define this to allow use of alloca instead of malloc for temporary buffers
// Temporary buffers are used in two case:
// 1. FFT sizes that have "bad" factors. i.e. not 2,3 and 5
// 2. "in-place" FFTs.  Notice the quotes, since kissfft does not really do an in-place transform.
#include <alloca.h>
#define  KISS_FFT_TMP_ALLOC(nbytes) alloca(nbytes)
#define  KISS_FFT_TMP_FREE(ptr)
#else
#define  KISS_FFT_TMP_ALLOC(nbytes) KISS_FFT_MALLOC(nbytes)
#define  KISS_FFT_TMP_FREE(ptr) KISS_FFT_FREE(ptr)
#endif

#endif /* _kiss_fft_guts_h */

