# FFT Software Integration Guide

## Overview

This document describes how the kissfft library has been modified to use the
custom RISC-V ISA extensions added to the CV32A6 core for FFT acceleration.

For hardware details, see:
- [`docs/FFT_ACCELERATION.md`](../../../docs/FFT_ACCELERATION.md) — Full architecture doc
- [`core/FFT_ISA_EXTENSIONS.md`](../../../core/FFT_ISA_EXTENSIONS.md) — Hardware reference

---

## Modified Files

| File | Purpose |
|------|---------|
| `_kiss_fft_guts.h` | Inline assembly intrinsics + macro overrides |
| `kiss_fft.c` | Butterfly functions with HW instruction protocol |
| `fft_int16_main.c` | Twiddle preloading + FFT invocation |
| `Makefile` (parent) | Compile flags for feature tiers |

---

## Packed Complex Format

Each complex number `{re, im}` is packed into a single `uint32_t`:

```
 bits:  31         16  15          0
       ┌──────────────┬──────────────┐
       │   imag (i)   │   real (r)   │
       └──────────────┴──────────────┘
```

Helper functions in `_kiss_fft_guts.h`:

```c
// Pack kiss_fft_cpx → uint32_t
static inline uint32_t kissfft_pack_cpx(kiss_fft_cpx a) {
    return ((uint16_t)a.r) | ((uint32_t)(uint16_t)a.i << 16);
}

// Unpack uint32_t → kiss_fft_cpx
static inline kiss_fft_cpx kissfft_unpack_cpx(uint32_t v) {
    kiss_fft_cpx r;
    r.r = (int16_t)(v & 0xFFFFu);
    r.i = (int16_t)(v >> 16);
    return r;
}
```

---

## Intrinsics Reference

### Tier 1 — Packed Complex ALU (OpcodeCustom0 = 0x0B)

| Function | Encoding | Description |
|----------|----------|-------------|
| `kissfft_cadd_u32(a, b)` | `.insn r 0x0B, 0, 0x00` | Complex add |
| `kissfft_csub_u32(a, b)` | `.insn r 0x0B, 1, 0x00` | Complex sub |
| `kissfft_cmul_u32(a, b)` | `.insn r 0x0B, 2, 0x00` | Complex mul (Q15) |
| `kissfft_cmadd_u32(acc, a, b)` | `.insn r 0x0B, 3, 0x00` | Complex fused mul-add |
| `kissfft_cmsub_u32(acc, a, b)` | `.insn r 0x0B, 4, 0x00` | Complex fused mul-sub |

### Tier 2 — Butterfly-2 (OpcodeCustom1 = 0x2B)

| Function | Encoding | Description |
|----------|----------|-------------|
| `kissfft_bfy2_lo_u32(a, b, tw)` | `.insn r4 0x2B, 0, 0` | BFY2 low output (a + tw*b) |
| `kissfft_bfy2_hi_u32(dep)` | `.insn r4 0x2B, 1, 0` | BFY2 high output (a - tw*b) |

### Tier 3 — Butterfly-4 (OpcodeCustom2 = 0x5B)

| Function | Encoding | Description |
|----------|----------|-------------|
| `kissfft_bfy4_s0(a, b, t1)` | `.insn r4 0x5B, 0, 0` | Setup: load a, b, tw1 |
| `kissfft_bfy4_s1(c, d, t2)` | `.insn r4 0x5B, 1, 0` | Setup: load c, d, tw2 |
| `kissfft_bfy4_s2(t3)` | `.insn r4 0x5B, 2, 0` | Compute (tw3 via rs3) |
| `kissfft_bfy4_o0(dep)` | `.insn r4 0x5B, 3, 0` | Read output 0 |
| `kissfft_bfy4_o1(dep)` | `.insn r4 0x5B, 4, 0` | Read output 1 |
| `kissfft_bfy4_o2(dep)` | `.insn r4 0x5B, 5, 0` | Read output 2 |
| `kissfft_bfy4_o3(dep)` | `.insn r4 0x5B, 6, 0` | Read output 3 |

### Tier 4 — Twiddle Register File (OpcodeCustom0 = 0x0B)

| Function | Encoding | Description |
|----------|----------|-------------|
| `kissfft_twld_u32(idx, val)` | `.insn r 0x0B, 5, 0x00` | Write twiddle_mem[idx] = val |
| `kissfft_twcfg_u32(mode)` | `.insn r 0x0B, 6, 0x00` | Configure cache (bit0=en, bit1=scale) |

---

## BFY4 Usage Protocol

The BFY4 instructions **must** be called in this exact sequence:

```c
// 1. Load 4 complex values from memory
uint32_t a = kissfft_pack_cpx(*Fout);
uint32_t b = kissfft_pack_cpx(Fout[m]);
uint32_t c = kissfft_pack_cpx(Fout[m2]);
uint32_t d = kissfft_pack_cpx(Fout[m3]);

// 2. Get twiddle factors (indices or packed values)
uint32_t t1 = (uint32_t)tw1_idx;   // with twiddle register file
uint32_t t2 = (uint32_t)tw2_idx;
uint32_t t3 = (uint32_t)tw3_idx;

// 3. Setup phase (order matters: S0 before S1 before S2)
C_BFY4_S0(a, b, t1);   // Latches a, b, tw1
C_BFY4_S1(c, d, t2);   // Latches c, d, tw2
C_BFY4_S2(t3);          // Computes all 4 outputs using tw3

// 4. Read results (any order, all depend on t3)
uint32_t o0 = C_BFY4_O0(t3);   // Output 0
uint32_t o1 = C_BFY4_O1(t3);   // Output 1
uint32_t o2 = C_BFY4_O2(t3);   // Output 2
uint32_t o3 = C_BFY4_O3(t3);   // Output 3

// 5. Store results back
*Fout      = kissfft_unpack_cpx(o0);
Fout[m]    = kissfft_unpack_cpx(o1);
Fout[m2]   = kissfft_unpack_cpx(o2);
Fout[m3]   = kissfft_unpack_cpx(o3);
```

### Dependency Chain

The `dep` argument in `BFY4_O*` is passed to create a **data dependency** that
prevents the compiler from reordering the reads before the compute. All O0-O3
use `t3` as their dependency (not each other), so the hardware can service them
without inter-instruction stalls.

---

## Twiddle Register File Preloading

Before calling `kiss_fft()`, the twiddle factors must be loaded into hardware:

```c
// In fft_int16_main.c:
#ifdef KISSFFT_USE_TWIDDLE_RF
    kiss_fft_hw_load_twiddles(cfg);
#endif
```

The implementation in `kiss_fft.c`:

```c
void kiss_fft_hw_load_twiddles(kiss_fft_cfg st) {
    int nfft = st->nfft;
    for (int i = 0; i < nfft; ++i) {
        uint32_t v = kissfft_pack_cpx(st->twiddles[i]);
        kissfft_twld_u32((uint32_t)i, v);      // TWLD: twiddle_mem[i] = v
    }
    uint32_t mode = 1;                          // bit0 = cache enable
#ifdef KISSFFT_USE_BFY4_SCALE
    mode |= 2;                                  // bit1 = ÷4 scale enable
#endif
    kissfft_twcfg_u32(mode);                    // TWCFG
}
```

---

## Fused ÷4 Scaling

When `KISSFFT_USE_BFY4_SCALE` is defined:

1. The `TWCFG` instruction is called with bit1=1, enabling hardware ÷4 scaling
2. The `C_FIXDIV(*Fout, 4)` calls in `kf_bfly4()` are **removed** (guarded by
   `#ifndef KISSFFT_USE_BFY4_SCALE`)
3. The BFY4 hardware automatically applies `q15_div4()` to all inputs during
   `BFY4_S0` and `BFY4_S1`

This is semantically equivalent but saves ~8 instructions per butterfly iteration.

---

## Compile Flags Summary

| Flag | Tier | Effect |
|------|------|--------|
| `KISSFFT_USE_IM_CUSTOM` | 1 | Enable packed complex intrinsics |
| `KISSFFT_USE_BFY2` | 2 | Use BFY2/BFY2H in radix-2 butterfly |
| `KISSFFT_USE_BFY4` | 3 | Use BFY4 protocol in radix-4 butterfly |
| `KISSFFT_MANUAL_UNROLL4` | — | Manually unroll radix-4 loop ×2 |
| `KISSFFT_USE_TWIDDLE_RF` | 4 | Use hardware twiddle register file (indices) |
| `KISSFFT_USE_BFY4_SCALE` | 5 | Fuse ÷4 scaling into BFY4 hardware |
| `FIXED_POINT` | — | Use int16 fixed-point (set to 16) |

### Build Command

```bash
# From sw/app/
make fft
```

The Makefile automatically passes all flags to both the kissfft library and the
main application.

---

## Performance Notes

### Instruction Count Evolution (approximate, 512-point FFT)

| Configuration | Instructions |
|--------------|-------------|
| Baseline RV32IM (no extensions) | ~150,000+ |
| + Packed complex ops (Tier 1) | ~100,000 |
| + BFY2 (Tier 2) | ~90,000 |
| + BFY4 (Tier 3) | ~70,000-80,000 |
| + Twiddle register file (Tier 4) | Testing needed |
| + Fused scaling (Tier 5) | Testing needed |

> **Note**: Tiers 4 and 5 are implemented but performance needs to be validated
> by running `make fft` and checking the instruction/cycle counts.
