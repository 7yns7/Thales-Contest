/*
 *  Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
 *  This file is part of KISS FFT - https://github.com/mborgerding/kissfft
 *
 *  SPDX-License-Identifier: BSD-3-Clause
 *  See COPYING file for more information.
 */

#include "_kiss_fft_guts.h"
/* The guts header contains all the multiplication and addition macros that are defined for
 fixed or floating point complex numbers.  It also delares the kf_ internal functions.
 */

static void kf_bfly2(
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const kiss_fft_cfg st,
        int m
        )
{
    kiss_fft_cpx * Fout2;
    kiss_fft_cpx * tw1 = st->twiddles;
#ifndef KISSFFT_USE_IM_CUSTOM
    kiss_fft_cpx t;
#endif
    Fout2 = Fout + m;
    do{
        C_FIXDIV(*Fout,2); C_FIXDIV(*Fout2,2);

#ifdef KISSFFT_USE_BFY2
    {
        uint32_t a = kissfft_pack_cpx(*Fout);
        uint32_t b = kissfft_pack_cpx(*Fout2);
        uint32_t tw = kissfft_pack_cpx(*tw1);
        uint32_t lo = C_BFY2_LO(a, b, tw);
            uint32_t hi = C_BFY2_HI(lo);
        *Fout = kissfft_unpack_cpx(lo);
        *Fout2 = kissfft_unpack_cpx(hi);
    }
#elif defined(KISSFFT_USE_IM_CUSTOM)
    {
        kiss_fft_cpx f0 = *Fout;
        kiss_fft_cpx f1 = *Fout;
        C_MULSUB(f1, *Fout2, *tw1);
        C_MULACC(f0, *Fout2, *tw1);
        *Fout2 = f1;
        *Fout  = f0;
    }
#else
        C_MUL (t,  *Fout2 , *tw1);
        C_SUB( *Fout2 ,  *Fout , t );
        C_ADDTO( *Fout ,  t );
#endif
        tw1 += fstride;
        ++Fout2;
        ++Fout;
    }while (--m);
}

static void kf_bfly4(
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const kiss_fft_cfg st,
        const size_t m
        )
{
    kiss_fft_cpx *tw1,*tw2,*tw3;
    kiss_fft_cpx scratch[6];
    size_t k=m;
    const size_t m2=2*m;
    const size_t m3=3*m;
#ifdef KISSFFT_USE_TWIDDLE_RF
    size_t tw1_idx = 0;
    size_t tw2_idx = 0;
    size_t tw3_idx = 0;
    const size_t tw1_step = fstride;
    const size_t tw2_step = fstride*2;
    const size_t tw3_step = fstride*3;
#endif


    tw3 = tw2 = tw1 = st->twiddles;

#ifdef KISSFFT_MANUAL_UNROLL4
    while (k >= 2) {
    #ifndef KISSFFT_USE_BFY4_SCALE
        C_FIXDIV(*Fout,4); C_FIXDIV(Fout[m],4); C_FIXDIV(Fout[m2],4); C_FIXDIV(Fout[m3],4);
    #endif

    #ifndef KISSFFT_USE_BFY4
        C_MUL(scratch[0],Fout[m] , *tw1 );
        C_MUL(scratch[1],Fout[m2] , *tw2 );
        C_MUL(scratch[2],Fout[m3] , *tw3 );
    #endif

    #ifdef KISSFFT_USE_BFY4
        {
            uint32_t a = kissfft_pack_cpx(*Fout);
            uint32_t b = kissfft_pack_cpx(Fout[m]);
            uint32_t c = kissfft_pack_cpx(Fout[m2]);
            uint32_t d = kissfft_pack_cpx(Fout[m3]);
#ifdef KISSFFT_USE_TWIDDLE_RF
            uint32_t t1 = (uint32_t)tw1_idx;
            uint32_t t2 = (uint32_t)tw2_idx;
            uint32_t t3 = (uint32_t)tw3_idx;
#else
            uint32_t t1 = kissfft_pack_cpx(*tw1);
            uint32_t t2 = kissfft_pack_cpx(*tw2);
            uint32_t t3 = kissfft_pack_cpx(*tw3);
#endif
            C_BFY4_S0(a, b, t1);
            C_BFY4_S1(c, d, t2);
            C_BFY4_S2(t3);
            uint32_t o0 = C_BFY4_O0(t3);
            uint32_t o1 = C_BFY4_O1(t3);
            uint32_t o2 = C_BFY4_O2(t3);
            uint32_t o3 = C_BFY4_O3(t3);
            *Fout = kissfft_unpack_cpx(o0);
            Fout[m] = kissfft_unpack_cpx(o1);
            Fout[m2] = kissfft_unpack_cpx(o2);
            Fout[m3] = kissfft_unpack_cpx(o3);
        }
        tw1 += fstride;
        tw2 += fstride*2;
        tw3 += fstride*3;
    #ifdef KISSFFT_USE_TWIDDLE_RF
        tw1_idx += tw1_step;
        tw2_idx += tw2_step;
        tw3_idx += tw3_step;
    #endif
    #else
    #ifdef KISSFFT_USE_IM_CUSTOM
        {
            kiss_fft_cpx f0 = *Fout;
            kiss_fft_cpx f5 = *Fout;
            C_MULACC(f0, Fout[m2], *tw2);
            C_MULSUB(f5, Fout[m2], *tw2);
            *Fout = f0;
            scratch[5] = f5;
        }
    #else
        C_SUB( scratch[5] , *Fout, scratch[1] );
        C_ADDTO(*Fout, scratch[1]);
    #endif
        C_ADD( scratch[3] , scratch[0] , scratch[2] );
        C_SUB( scratch[4] , scratch[0] , scratch[2] );
        C_SUB( Fout[m2], *Fout, scratch[3] );
        tw1 += fstride;
        tw2 += fstride*2;
        tw3 += fstride*3;
        C_ADDTO( *Fout , scratch[3] );

        if(st->inverse) {
            Fout[m].r = scratch[5].r - scratch[4].i;
            Fout[m].i = scratch[5].i + scratch[4].r;
            Fout[m3].r = scratch[5].r + scratch[4].i;
            Fout[m3].i = scratch[5].i - scratch[4].r;
        }else{
            Fout[m].r = scratch[5].r + scratch[4].i;
            Fout[m].i = scratch[5].i - scratch[4].r;
            Fout[m3].r = scratch[5].r - scratch[4].i;
            Fout[m3].i = scratch[5].i + scratch[4].r;
        }
#endif
        ++Fout;

    #ifndef KISSFFT_USE_BFY4_SCALE
        C_FIXDIV(*Fout,4); C_FIXDIV(Fout[m],4); C_FIXDIV(Fout[m2],4); C_FIXDIV(Fout[m3],4);
    #endif

    #ifndef KISSFFT_USE_BFY4
        C_MUL(scratch[0],Fout[m] , *tw1 );
        C_MUL(scratch[1],Fout[m2] , *tw2 );
        C_MUL(scratch[2],Fout[m3] , *tw3 );
    #endif

    #ifdef KISSFFT_USE_BFY4
        {
            uint32_t a = kissfft_pack_cpx(*Fout);
            uint32_t b = kissfft_pack_cpx(Fout[m]);
            uint32_t c = kissfft_pack_cpx(Fout[m2]);
            uint32_t d = kissfft_pack_cpx(Fout[m3]);
#ifdef KISSFFT_USE_TWIDDLE_RF
            uint32_t t1 = (uint32_t)tw1_idx;
            uint32_t t2 = (uint32_t)tw2_idx;
            uint32_t t3 = (uint32_t)tw3_idx;
#else
            uint32_t t1 = kissfft_pack_cpx(*tw1);
            uint32_t t2 = kissfft_pack_cpx(*tw2);
            uint32_t t3 = kissfft_pack_cpx(*tw3);
#endif
            C_BFY4_S0(a, b, t1);
            C_BFY4_S1(c, d, t2);
            C_BFY4_S2(t3);
            uint32_t o0 = C_BFY4_O0(t3);
            uint32_t o1 = C_BFY4_O1(t3);
            uint32_t o2 = C_BFY4_O2(t3);
            uint32_t o3 = C_BFY4_O3(t3);
            *Fout = kissfft_unpack_cpx(o0);
            Fout[m] = kissfft_unpack_cpx(o1);
            Fout[m2] = kissfft_unpack_cpx(o2);
            Fout[m3] = kissfft_unpack_cpx(o3);
        }
        tw1 += fstride;
        tw2 += fstride*2;
        tw3 += fstride*3;
    #ifdef KISSFFT_USE_TWIDDLE_RF
        tw1_idx += tw1_step;
        tw2_idx += tw2_step;
        tw3_idx += tw3_step;
    #endif
    #else
    #ifdef KISSFFT_USE_IM_CUSTOM
        {
            kiss_fft_cpx f0 = *Fout;
            kiss_fft_cpx f5 = *Fout;
            C_MULACC(f0, Fout[m2], *tw2);
            C_MULSUB(f5, Fout[m2], *tw2);
            *Fout = f0;
            scratch[5] = f5;
        }
    #else
        C_SUB( scratch[5] , *Fout, scratch[1] );
        C_ADDTO(*Fout, scratch[1]);
    #endif
        C_ADD( scratch[3] , scratch[0] , scratch[2] );
        C_SUB( scratch[4] , scratch[0] , scratch[2] );
        C_SUB( Fout[m2], *Fout, scratch[3] );
        tw1 += fstride;
        tw2 += fstride*2;
        tw3 += fstride*3;
        C_ADDTO( *Fout , scratch[3] );

        if(st->inverse) {
            Fout[m].r = scratch[5].r - scratch[4].i;
            Fout[m].i = scratch[5].i + scratch[4].r;
            Fout[m3].r = scratch[5].r + scratch[4].i;
            Fout[m3].i = scratch[5].i - scratch[4].r;
        }else{
            Fout[m].r = scratch[5].r + scratch[4].i;
            Fout[m].i = scratch[5].i - scratch[4].r;
            Fout[m3].r = scratch[5].r - scratch[4].i;
            Fout[m3].i = scratch[5].i + scratch[4].r;
        }
#endif
        ++Fout;
        k -= 2;
    }

    if (k) {
    #ifndef KISSFFT_USE_BFY4_SCALE
        C_FIXDIV(*Fout,4); C_FIXDIV(Fout[m],4); C_FIXDIV(Fout[m2],4); C_FIXDIV(Fout[m3],4);
    #endif

    #ifndef KISSFFT_USE_BFY4
        C_MUL(scratch[0],Fout[m] , *tw1 );
        C_MUL(scratch[1],Fout[m2] , *tw2 );
        C_MUL(scratch[2],Fout[m3] , *tw3 );
    #endif

    #ifdef KISSFFT_USE_BFY4
        {
            uint32_t a = kissfft_pack_cpx(*Fout);
            uint32_t b = kissfft_pack_cpx(Fout[m]);
            uint32_t c = kissfft_pack_cpx(Fout[m2]);
            uint32_t d = kissfft_pack_cpx(Fout[m3]);
#ifdef KISSFFT_USE_TWIDDLE_RF
            uint32_t t1 = (uint32_t)tw1_idx;
            uint32_t t2 = (uint32_t)tw2_idx;
            uint32_t t3 = (uint32_t)tw3_idx;
#else
            uint32_t t1 = kissfft_pack_cpx(*tw1);
            uint32_t t2 = kissfft_pack_cpx(*tw2);
            uint32_t t3 = kissfft_pack_cpx(*tw3);
#endif
            C_BFY4_S0(a, b, t1);
            C_BFY4_S1(c, d, t2);
            C_BFY4_S2(t3);
            uint32_t o0 = C_BFY4_O0(t3);
            uint32_t o1 = C_BFY4_O1(t3);
            uint32_t o2 = C_BFY4_O2(t3);
            uint32_t o3 = C_BFY4_O3(t3);
            *Fout = kissfft_unpack_cpx(o0);
            Fout[m] = kissfft_unpack_cpx(o1);
            Fout[m2] = kissfft_unpack_cpx(o2);
            Fout[m3] = kissfft_unpack_cpx(o3);
        }
        tw1 += fstride;
        tw2 += fstride*2;
        tw3 += fstride*3;
    #ifdef KISSFFT_USE_TWIDDLE_RF
        tw1_idx += tw1_step;
        tw2_idx += tw2_step;
        tw3_idx += tw3_step;
    #endif
    #else
    #ifdef KISSFFT_USE_IM_CUSTOM
        {
            kiss_fft_cpx f0 = *Fout;
            kiss_fft_cpx f5 = *Fout;
            C_MULACC(f0, Fout[m2], *tw2);
            C_MULSUB(f5, Fout[m2], *tw2);
            *Fout = f0;
            scratch[5] = f5;
        }
    #else
        C_SUB( scratch[5] , *Fout, scratch[1] );
        C_ADDTO(*Fout, scratch[1]);
    #endif
        C_ADD( scratch[3] , scratch[0] , scratch[2] );
        C_SUB( scratch[4] , scratch[0] , scratch[2] );
        C_SUB( Fout[m2], *Fout, scratch[3] );
        tw1 += fstride;
        tw2 += fstride*2;
        tw3 += fstride*3;
        C_ADDTO( *Fout , scratch[3] );

        if(st->inverse) {
            Fout[m].r = scratch[5].r - scratch[4].i;
            Fout[m].i = scratch[5].i + scratch[4].r;
            Fout[m3].r = scratch[5].r + scratch[4].i;
            Fout[m3].i = scratch[5].i - scratch[4].r;
        }else{
            Fout[m].r = scratch[5].r + scratch[4].i;
            Fout[m].i = scratch[5].i - scratch[4].r;
            Fout[m3].r = scratch[5].r - scratch[4].i;
            Fout[m3].i = scratch[5].i + scratch[4].r;
        }
#endif
        ++Fout;
    }
#else
    do {
    #ifndef KISSFFT_USE_BFY4_SCALE
        C_FIXDIV(*Fout,4); C_FIXDIV(Fout[m],4); C_FIXDIV(Fout[m2],4); C_FIXDIV(Fout[m3],4);
    #endif

    #ifndef KISSFFT_USE_BFY4
        C_MUL(scratch[0],Fout[m] , *tw1 );
        C_MUL(scratch[1],Fout[m2] , *tw2 );
        C_MUL(scratch[2],Fout[m3] , *tw3 );
    #endif

    #ifdef KISSFFT_USE_BFY4
        {
            uint32_t a = kissfft_pack_cpx(*Fout);
            uint32_t b = kissfft_pack_cpx(Fout[m]);
            uint32_t c = kissfft_pack_cpx(Fout[m2]);
            uint32_t d = kissfft_pack_cpx(Fout[m3]);
#ifdef KISSFFT_USE_TWIDDLE_RF
            uint32_t t1 = (uint32_t)tw1_idx;
            uint32_t t2 = (uint32_t)tw2_idx;
            uint32_t t3 = (uint32_t)tw3_idx;
#else
            uint32_t t1 = kissfft_pack_cpx(*tw1);
            uint32_t t2 = kissfft_pack_cpx(*tw2);
            uint32_t t3 = kissfft_pack_cpx(*tw3);
#endif
            C_BFY4_S0(a, b, t1);
            C_BFY4_S1(c, d, t2);
            C_BFY4_S2(t3);
            uint32_t o0 = C_BFY4_O0(t3);
            uint32_t o1 = C_BFY4_O1(t3);
            uint32_t o2 = C_BFY4_O2(t3);
            uint32_t o3 = C_BFY4_O3(t3);
            *Fout = kissfft_unpack_cpx(o0);
            Fout[m] = kissfft_unpack_cpx(o1);
            Fout[m2] = kissfft_unpack_cpx(o2);
            Fout[m3] = kissfft_unpack_cpx(o3);
        }
        tw1 += fstride;
        tw2 += fstride*2;
        tw3 += fstride*3;
    #ifdef KISSFFT_USE_TWIDDLE_RF
        tw1_idx += tw1_step;
        tw2_idx += tw2_step;
        tw3_idx += tw3_step;
    #endif
    #else
    #ifdef KISSFFT_USE_IM_CUSTOM
        {
            kiss_fft_cpx f0 = *Fout;
            kiss_fft_cpx f5 = *Fout;
            C_MULACC(f0, Fout[m2], *tw2);
            C_MULSUB(f5, Fout[m2], *tw2);
            *Fout = f0;
            scratch[5] = f5;
        }
    #else
        C_SUB( scratch[5] , *Fout, scratch[1] );
        C_ADDTO(*Fout, scratch[1]);
    #endif
        C_ADD( scratch[3] , scratch[0] , scratch[2] );
        C_SUB( scratch[4] , scratch[0] , scratch[2] );
        C_SUB( Fout[m2], *Fout, scratch[3] );
        tw1 += fstride;
        tw2 += fstride*2;
        tw3 += fstride*3;
        C_ADDTO( *Fout , scratch[3] );

        if(st->inverse) {
            Fout[m].r = scratch[5].r - scratch[4].i;
            Fout[m].i = scratch[5].i + scratch[4].r;
            Fout[m3].r = scratch[5].r + scratch[4].i;
            Fout[m3].i = scratch[5].i - scratch[4].r;
        }else{
            Fout[m].r = scratch[5].r + scratch[4].i;
            Fout[m].i = scratch[5].i - scratch[4].r;
            Fout[m3].r = scratch[5].r - scratch[4].i;
            Fout[m3].i = scratch[5].i + scratch[4].r;
        }
#endif
        ++Fout;
    }while(--k);
#endif
}

static void kf_bfly3(
         kiss_fft_cpx * Fout,
         const size_t fstride,
         const kiss_fft_cfg st,
         size_t m
         )
{
     size_t k=m;
     const size_t m2 = 2*m;
     kiss_fft_cpx *tw1,*tw2;
     kiss_fft_cpx scratch[5];
     kiss_fft_cpx epi3;
     epi3 = st->twiddles[fstride*m];

     tw1=tw2=st->twiddles;

     do{
         C_FIXDIV(*Fout,3); C_FIXDIV(Fout[m],3); C_FIXDIV(Fout[m2],3);

         C_MUL(scratch[1],Fout[m] , *tw1);
         C_MUL(scratch[2],Fout[m2] , *tw2);

         C_ADD(scratch[3],scratch[1],scratch[2]);
         C_SUB(scratch[0],scratch[1],scratch[2]);
         tw1 += fstride;
         tw2 += fstride*2;

         Fout[m].r = Fout->r - HALF_OF(scratch[3].r);
         Fout[m].i = Fout->i - HALF_OF(scratch[3].i);

         C_MULBYSCALAR( scratch[0] , epi3.i );

         C_ADDTO(*Fout,scratch[3]);

         Fout[m2].r = Fout[m].r + scratch[0].i;
         Fout[m2].i = Fout[m].i - scratch[0].r;

         Fout[m].r -= scratch[0].i;
         Fout[m].i += scratch[0].r;

         ++Fout;
     }while(--k);
}

static void kf_bfly5(
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const kiss_fft_cfg st,
        int m
        )
{
    kiss_fft_cpx *Fout0,*Fout1,*Fout2,*Fout3,*Fout4;
    int u;
    kiss_fft_cpx scratch[13];
    kiss_fft_cpx * twiddles = st->twiddles;
    kiss_fft_cpx *tw;
    kiss_fft_cpx ya,yb;
    ya = twiddles[fstride*m];
    yb = twiddles[fstride*2*m];

    Fout0=Fout;
    Fout1=Fout0+m;
    Fout2=Fout0+2*m;
    Fout3=Fout0+3*m;
    Fout4=Fout0+4*m;

    tw=st->twiddles;
    for ( u=0; u<m; ++u ) {
        C_FIXDIV( *Fout0,5); C_FIXDIV( *Fout1,5); C_FIXDIV( *Fout2,5); C_FIXDIV( *Fout3,5); C_FIXDIV( *Fout4,5);
        scratch[0] = *Fout0;

        C_MUL(scratch[1] ,*Fout1, tw[u*fstride]);
        C_MUL(scratch[2] ,*Fout2, tw[2*u*fstride]);
        C_MUL(scratch[3] ,*Fout3, tw[3*u*fstride]);
        C_MUL(scratch[4] ,*Fout4, tw[4*u*fstride]);

        C_ADD( scratch[7],scratch[1],scratch[4]);
        C_SUB( scratch[10],scratch[1],scratch[4]);
        C_ADD( scratch[8],scratch[2],scratch[3]);
        C_SUB( scratch[9],scratch[2],scratch[3]);

    #ifdef KISSFFT_USE_IM_CUSTOM
        C_ADDTO(*Fout0, scratch[7]);
        C_ADDTO(*Fout0, scratch[8]);
    #else
        Fout0->r += scratch[7].r + scratch[8].r;
        Fout0->i += scratch[7].i + scratch[8].i;
    #endif

        scratch[5].r = scratch[0].r + S_MUL(scratch[7].r,ya.r) + S_MUL(scratch[8].r,yb.r);
        scratch[5].i = scratch[0].i + S_MUL(scratch[7].i,ya.r) + S_MUL(scratch[8].i,yb.r);

        scratch[6].r =  S_MUL(scratch[10].i,ya.i) + S_MUL(scratch[9].i,yb.i);
        scratch[6].i = -S_MUL(scratch[10].r,ya.i) - S_MUL(scratch[9].r,yb.i);

        C_SUB(*Fout1,scratch[5],scratch[6]);
        C_ADD(*Fout4,scratch[5],scratch[6]);

        scratch[11].r = scratch[0].r + S_MUL(scratch[7].r,yb.r) + S_MUL(scratch[8].r,ya.r);
        scratch[11].i = scratch[0].i + S_MUL(scratch[7].i,yb.r) + S_MUL(scratch[8].i,ya.r);
        scratch[12].r = - S_MUL(scratch[10].i,yb.i) + S_MUL(scratch[9].i,ya.i);
        scratch[12].i = S_MUL(scratch[10].r,yb.i) - S_MUL(scratch[9].r,ya.i);

        C_ADD(*Fout2,scratch[11],scratch[12]);
        C_SUB(*Fout3,scratch[11],scratch[12]);

        ++Fout0;++Fout1;++Fout2;++Fout3;++Fout4;
    }
}

/* perform the butterfly for one stage of a mixed radix FFT */
static void kf_bfly_generic(
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const kiss_fft_cfg st,
        int m,
        int p
        )
{
    int u,k,q1,q;
    kiss_fft_cpx * twiddles = st->twiddles;
    kiss_fft_cpx t;
    int Norig = st->nfft;

    kiss_fft_cpx * scratch = (kiss_fft_cpx*)KISS_FFT_TMP_ALLOC(sizeof(kiss_fft_cpx)*p);
    if (scratch == NULL){
        KISS_FFT_ERROR("Memory allocation failed.");
        return;
    }

    for ( u=0; u<m; ++u ) {
        k=u;
        for ( q1=0 ; q1<p ; ++q1 ) {
            scratch[q1] = Fout[ k  ];
            C_FIXDIV(scratch[q1],p);
            k += m;
        }

        k=u;
        for ( q1=0 ; q1<p ; ++q1 ) {
            int twidx=0;
            Fout[ k ] = scratch[0];
            for (q=1;q<p;++q ) {
                twidx += fstride * k;
                if (twidx>=Norig) twidx-=Norig;
#ifdef KISSFFT_USE_IM_CUSTOM
                C_MULACC( Fout[ k ] , scratch[q] , twiddles[twidx] );
#else
                C_MUL(t,scratch[q] , twiddles[twidx] );
                C_ADDTO( Fout[ k ] ,t);
#endif
            }
            k += m;
        }
    }
    KISS_FFT_TMP_FREE(scratch);
}

static
void kf_work(
        kiss_fft_cpx * Fout,
        const kiss_fft_cpx * f,
        const size_t fstride,
        int in_stride,
        int * factors,
        const kiss_fft_cfg st
        )
{
    kiss_fft_cpx * Fout_beg=Fout;
    const int p=*factors++; /* the radix  */
    const int m=*factors++; /* stage's fft length/p */
    const kiss_fft_cpx * Fout_end = Fout + p*m;

#ifdef _OPENMP
    // use openmp extensions at the
    // top-level (not recursive)
    if (fstride==1 && p<=5 && m!=1)
    {
        int k;

        // execute the p different work units in different threads
#       pragma omp parallel for
        for (k=0;k<p;++k)
            kf_work( Fout +k*m, f+ fstride*in_stride*k,fstride*p,in_stride,factors,st);
        // all threads have joined by this point

        switch (p) {
            case 2: kf_bfly2(Fout,fstride,st,m); break;
            case 3: kf_bfly3(Fout,fstride,st,m); break;
            case 4: kf_bfly4(Fout,fstride,st,m); break;
            case 5: kf_bfly5(Fout,fstride,st,m); break;
            default: kf_bfly_generic(Fout,fstride,st,m,p); break;
        }
        return;
    }
#endif

    if (m==1) {
        do{
            *Fout = *f;
            f += fstride*in_stride;
        }while(++Fout != Fout_end );
    }else{
        do{
            // recursive call:
            // DFT of size m*p performed by doing
            // p instances of smaller DFTs of size m,
            // each one takes a decimated version of the input
            kf_work( Fout , f, fstride*p, in_stride, factors,st);
            f += fstride*in_stride;
        }while( (Fout += m) != Fout_end );
    }

    Fout=Fout_beg;

    // recombine the p smaller DFTs
    switch (p) {
        case 2: kf_bfly2(Fout,fstride,st,m); break;
        case 3: kf_bfly3(Fout,fstride,st,m); break;
        case 4: kf_bfly4(Fout,fstride,st,m); break;
        case 5: kf_bfly5(Fout,fstride,st,m); break;
        default: kf_bfly_generic(Fout,fstride,st,m,p); break;
    }
}

/*  facbuf is populated by p1,m1,p2,m2, ...
    where
    p[i] * m[i] = m[i-1]
    m0 = n                  */
static
void kf_factor(int n,int * facbuf)
{
    int p=4;
    double floor_sqrt;
    floor_sqrt = floor( sqrt((double)n) );

    /*factor out powers of 4, powers of 2, then any remaining primes */
    do {
        while (n % p) {
            switch (p) {
                case 4: p = 2; break;
                case 2: p = 3; break;
                default: p += 2; break;
            }
            if (p > floor_sqrt)
                p = n;          /* no more factors, skip to end */
        }
        n /= p;
        *facbuf++ = p;
        *facbuf++ = n;
    } while (n > 1);
}

/*
 *
 * User-callable function to allocate all necessary storage space for the fft.
 *
 * The return value is a contiguous block of memory, allocated with malloc.  As such,
 * It can be freed with free(), rather than a kiss_fft-specific function.
 * */
extern kiss_fft_cpx g_twiddles[512];
extern int g_factors[10];
kiss_fft_cfg kiss_fft_alloc(int nfft,int inverse_fft,void * mem,size_t * lenmem )
{

    KISS_FFT_ALIGN_CHECK(mem)

    kiss_fft_cfg st=NULL;
    size_t memneeded = KISS_FFT_ALIGN_SIZE_UP(sizeof(struct kiss_fft_state)
        + sizeof(kiss_fft_cpx)*(nfft-1)); /* twiddle factors*/

    if ( lenmem==NULL ) {
        st = ( kiss_fft_cfg)KISS_FFT_MALLOC( memneeded );
    }else{
        if (mem != NULL && *lenmem >= memneeded)
            st = (kiss_fft_cfg)mem;
        *lenmem = memneeded;
    }
    if (st) {
        int i;
        st->nfft=nfft;
        st->inverse = inverse_fft;

	// THE C_EXP OPERATIONS ARE PRECALCULATED

	st->twiddles = (kiss_fft_cpx *)&g_twiddles;
	st->factors = (int *)&g_factors;

        // for (i=0;i<nfft;++i) {
        //     const double pi=3.141592653589793238462643383279502884197169399375105820974944;
        //     double phase = -2*pi*i / nfft;
        //     if (st->inverse)
        //         phase *= -1;
        //     kf_cexp(st->twiddles+i, phase );
        // }

        // kf_factor(nfft,st->factors);
    }

    return st;
}

#ifdef KISSFFT_USE_TWIDDLE_RF
void kiss_fft_hw_load_twiddles(kiss_fft_cfg st)
{
#if defined(KISSFFT_USE_IM_CUSTOM) && (FIXED_POINT==16) && defined(__riscv)
    int nfft = st->nfft;
    for (int i = 0; i < nfft; ++i) {
        uint32_t v = kissfft_pack_cpx(st->twiddles[i]);
        kissfft_twld_u32((uint32_t)i, v);
    }
    uint32_t mode = 1;
#ifdef KISSFFT_USE_BFY4_SCALE
    mode |= 2;
#endif
    kissfft_twcfg_u32(mode);
#else
    (void)st;
#endif
}
#endif


void kiss_fft_stride(kiss_fft_cfg st,const kiss_fft_cpx *fin,kiss_fft_cpx *fout,int in_stride)
{
    if (fin == fout) {
        //NOTE: this is not really an in-place FFT algorithm.
        //It just performs an out-of-place FFT into a temp buffer
        if (fout == NULL){
            KISS_FFT_ERROR("fout buffer NULL.");
        return;
        }

        kiss_fft_cpx * tmpbuf = (kiss_fft_cpx*)KISS_FFT_TMP_ALLOC( sizeof(kiss_fft_cpx)*st->nfft);
        if (tmpbuf == NULL){
            KISS_FFT_ERROR("Memory allocation error.");
        return;
        }

        kf_work(tmpbuf,fin,1,in_stride, st->factors,st);
        memcpy(fout,tmpbuf,sizeof(kiss_fft_cpx)*st->nfft);
        KISS_FFT_TMP_FREE(tmpbuf);
    }else{
        kf_work( fout, fin, 1,in_stride, st->factors,st );
    }
}

void kiss_fft(kiss_fft_cfg cfg,const kiss_fft_cpx *fin,kiss_fft_cpx *fout)
{
    kiss_fft_stride(cfg,fin,fout,1);
}


void kiss_fft_cleanup(void)
{
    // nothing needed any more
}

int kiss_fft_next_fast_size(int n)
{
    while(1) {
        int m=n;
        while ( (m%2) == 0 ) m/=2;
        while ( (m%3) == 0 ) m/=3;
        while ( (m%5) == 0 ) m/=5;
        if (m<=1)
            break; /* n is completely factorable by twos, threes, and fives */
        n++;
    }
    return n;
}
