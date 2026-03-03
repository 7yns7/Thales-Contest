# FFT Acceleration — Architecture & Design Document

## Table of Contents

1. [Overview](#overview)
2. [Contest Constraints](#contest-constraints)
3. [High-Level Architecture](#high-level-architecture)
4. [Custom ISA Extensions Summary](#custom-isa-extensions-summary)
5. [Tier 1 — Packed Complex ALU Operations](#tier-1--packed-complex-alu-operations)
6. [Tier 2 — Radix-2 Butterfly Unit (BFY2)](#tier-2--radix-2-butterfly-unit-bfy2)
7. [Tier 3 — Radix-4 Butterfly Unit (BFY4)](#tier-3--radix-4-butterfly-unit-bfy4)
8. [Tier 4 — Twiddle Factor Cache](#tier-4--twiddle-factor-cache)
9. [Tier 5 — Fused Scaling (÷4 Fusion)](#tier-5--fused-scaling-4-fusion)
10. [Data Flow Diagrams](#data-flow-diagrams)
11. [Instruction Encoding Reference](#instruction-encoding-reference)
12. [Software Integration](#software-integration)
13. [Build Flags](#build-flags)
14. [Files Modified](#files-modified)

---

## Overview

This document describes all hardware (RTL) and software modifications made to the
**CV32A6** (CVA6) RISC-V core to accelerate a **512-point fixed-point (int16, Q15)
radix-4 FFT** using the [kissfft](https://github.com/mborgerding/kissfft) library.

The approach is **incremental**: five tiers of ISA extensions are layered on top
of the base RV32IM core, each one building on the previous:

```
 ┌─────────────────────────────────────────────────────────┐
 │          Tier 5: Fused ÷4 Scaling in BFY4               │
 ├─────────────────────────────────────────────────────────┤
 │          Tier 4: Twiddle Factor Cache (512 entries)      │
 ├─────────────────────────────────────────────────────────┤
 │          Tier 3: Radix-4 Butterfly (BFY4_S0..S2, O0..O3)│
 ├─────────────────────────────────────────────────────────┤
 │          Tier 2: Radix-2 Butterfly (BFY2 / BFY2H)       │
 ├─────────────────────────────────────────────────────────┤
 │          Tier 1: Packed Complex ALU (CADD/CSUB/CMUL/...) │
 ├─────────────────────────────────────────────────────────┤
 │                  Base CV32A6 (RV32IM_Zicsr)              │
 └─────────────────────────────────────────────────────────┘
```

---

## Contest Constraints

| Rule | Constraint |
|------|-----------|
| Core | CV32A6 (single-issue, in-order) or CVXIF coprocessor |
| ISA base | RV32IM_Zicsr |
| Extensions | Custom ISA extensions allowed (fine-grained ops) |
| Forbidden | Hardwired accelerator for full or partial FFT ≥ 8 points |
| Cache | No increase in cache size allowed |
| Frequency | Must not decrease by more than 20% |
| Correctness | Bit-exact match with reference output |
| Other apps | CoreMark and MNIST must still run correctly |

---

## High-Level Architecture

```
                             CV32A6 Core (Modified)
  ┌──────────────────────────────────────────────────────────────────┐
  │                                                                  │
  │  ┌─────────┐    ┌───────────┐    ┌──────────────────────────┐   │
  │  │ FRONTEND │───▶│  DECODER  │───▶│      ISSUE STAGE         │   │
  │  └─────────┘    │ (Custom0  │    │  (reads rs1,rs2,rs3)     │   │
  │                 │  Custom1  │    └────┬───────────┬──────────┘   │
  │                 │  Custom2) │         │           │              │
  │                 └───────────┘         ▼           ▼              │
  │                               ┌────────────┐ ┌────────────────┐ │
  │                               │    ALU     │ │      MULT      │ │
  │                               │            │ │                │ │
  │                               │ ● CADD16   │ │ ● Multiplier   │ │
  │                               │ ● CSUB16   │ │ ● Divider      │ │
  │                               │ ● CMUL16   │ │ ● BFY2 unit    │ │
  │                               │ ● CMADD16  │ │ ● BFY4 unit    │ │
  │                               │ ● CMSUB16  │ │ ● Twiddle register file│ │
  │                               │            │ │ ● TWLD/TWCFG   │ │
  │                               └────────────┘ └────────────────┘ │
  │                                       │           │              │
  │                                       ▼           ▼              │
  │                               ┌──────────────────────────┐      │
  │                               │     COMMIT / WRITEBACK    │      │
  │                               └──────────────────────────┘      │
  └──────────────────────────────────────────────────────────────────┘
```

**Key insight**: All custom instructions use existing functional units (ALU or MULT)
and the standard register file. No new pipeline stages or external accelerators are
added — this keeps the design contest-compliant.

---

## Custom ISA Extensions Summary

| Instruction | Opcode | FU | Latency | Description |
|------------|--------|-----|---------|-------------|
| `CADD16` | Custom-0 (0x0B), f3=0 | ALU | 1 cyc | Packed complex add |
| `CSUB16` | Custom-0 (0x0B), f3=1 | ALU | 1 cyc | Packed complex sub |
| `CMUL16` | Custom-0 (0x0B), f3=2 | ALU | 1 cyc | Packed complex mul (Q15) |
| `CMADD16` | Custom-0 (0x0B), f3=3 | ALU | 1 cyc | Packed complex mul-add (Q15) |
| `CMSUB16` | Custom-0 (0x0B), f3=4 | ALU | 1 cyc | Packed complex mul-sub (Q15) |
| `TWLD` | Custom-0 (0x0B), f3=5 | MULT | 1 cyc | Twiddle register file load |
| `TWCFG` | Custom-0 (0x0B), f3=6 | MULT | 1 cyc | Twiddle register file config |
| `BFY2` | Custom-1 (0x2B), f3=0 | MULT | 2 cyc | Butterfly-2 low output |
| `BFY2H` | Custom-1 (0x2B), f3=1 | MULT | 1 cyc | Butterfly-2 high output |
| `BFY4_S0` | Custom-2 (0x5B), f3=0 | MULT | 1 cyc | BFY4 setup: load a,b,tw1 |
| `BFY4_S1` | Custom-2 (0x5B), f3=1 | MULT | 1 cyc | BFY4 setup: load c,d,tw2 |
| `BFY4_S2` | Custom-2 (0x5B), f3=2 | MULT | 1 cyc | BFY4 compute (tw3 via rs3) |
| `BFY4_O0` | Custom-2 (0x5B), f3=3 | MULT | 1 cyc | Read BFY4 output 0 |
| `BFY4_O1` | Custom-2 (0x5B), f3=4 | MULT | 1 cyc | Read BFY4 output 1 |
| `BFY4_O2` | Custom-2 (0x5B), f3=5 | MULT | 1 cyc | Read BFY4 output 2 |
| `BFY4_O3` | Custom-2 (0x5B), f3=6 | MULT | 1 cyc | Read BFY4 output 3 |

---

## Tier 1 — Packed Complex ALU Operations

### Concept

Each complex number `z = (re, im)` is packed into a single 32-bit register:

```
  31         16 15          0
 ┌─────────────┬─────────────┐
 │   imag (i)  │   real (r)  │   ← int16 + int16 = 32 bits
 └─────────────┴─────────────┘
```

This halves the register pressure: one complex value = one GPR.

### Instructions

All Tier 1 ops use **R-type** encoding on `OpcodeCustom0` (0x0B):

| Mnemonic | funct3 | Semantics |
|----------|--------|-----------|
| `CADD16` | 000 | `rd = (rs1.i + rs2.i, rs1.r + rs2.r)` |
| `CSUB16` | 001 | `rd = (rs1.i - rs2.i, rs1.r - rs2.r)` |
| `CMUL16` | 010 | `rd = Q15_round(rs1 × rs2)` — complex multiply |
| `CMADD16` | 011 | `rd = Q15_round(rs1 × rs2) + rd_old` — fused mul-add (rs3=rd) |
| `CMSUB16` | 100 | `rd = rd_old - Q15_round(rs1 × rs2)` — fused mul-sub (rs3=rd) |

**Q15 rounding**: `q15_round(x) = (x + 16384) >> 15`

### Hardware location: `core/alu.sv`

The combinational datapath computes all five results in parallel and the result
MUX selects the correct one based on `fu_data_i.operation`. For `CMADD16`/`CMSUB16`,
the accumulator comes from the `imm` field (which carries the rs3 value via
`MUX_RD_RS3` in the decoder).

```
          rs1 packed          rs2 packed
             │                    │
     ┌───────┴───────┐   ┌───────┴───────┐
     │  split 16|16  │   │  split 16|16  │
     └───┬───────┬───┘   └───┬───────┬───┘
       a_r     a_i         b_r     b_i
         │       │           │       │
         ├───────┼───────────┼───────┤
         │       │   4 muls  │       │
         │       ▼           ▼       │
         │   ┌──────────────────┐    │
         │   │ a_r*b_r  a_i*b_i │    │
         │   │ a_r*b_i  a_i*b_r │    │
         │   └────────┬─────────┘    │
         │      ┌─────┴──────┐       │
         │      │ sub / add  │       │
         │      │ Q15 round  │       │
         │      └─────┬──────┘       │
         │         mul_packed        │
         ├─────► add_packed          │
         ├─────► sub_packed          │
         │                           │
         ▼         Result MUX        ▼
```

---

## Tier 2 — Radix-2 Butterfly Unit (BFY2)

### Concept

A radix-2 butterfly computes:
```
out0 = a + twiddle * b
out1 = a - twiddle * b
```

In software this requires: 1 complex multiply + 2 complex add/sub = ~10 scalar ops.
The BFY2 instruction fuses this into a **2-cycle pipelined** operation in the MULT unit.

### Instructions

Uses **R4-type** encoding on `OpcodeCustom1` (0x2B):

| Mnemonic | funct3 | Operands | Result |
|----------|--------|----------|--------|
| `BFY2` | 000 | rs1=a, rs2=b, rs3=twiddle | rd = out0 (a + tw*b) |
| `BFY2H` | 001 | rs1=dependency | rd = out1 (a - tw*b) |

### Pipeline

```
 Cycle 0 (BFY2 issued):
   ┌─────────────────────────────┐
   │ Capture a, b, twiddle       │
   │ Compute 4 partial products: │
   │   b_r*t_r, b_i*t_i,        │
   │   b_r*t_i, b_i*t_r         │
   │ Store in pipeline regs      │
   └─────────────────────────────┘
                 │
                 ▼
 Cycle 1 (result available for BFY2):
   ┌─────────────────────────────┐
   │ real = b_rr - b_ii          │
   │ imag = b_ri + b_ir          │
   │ Q15 round both              │
   │ out0 = a + rounded          │
   │ out1 = a - rounded (stored) │
   │ → rd = out0                 │
   └─────────────────────────────┘
                 │
                 ▼
 Cycle 2+ (BFY2H issued):
   ┌─────────────────────────────┐
   │ → rd = out1 (from register) │
   └─────────────────────────────┘
```

### Hardware location: `core/mult.sv` (BFY2 section)

---

## Tier 3 — Radix-4 Butterfly Unit (BFY4)

### Concept

A complete radix-4 butterfly computes 4 outputs from 4 inputs and 3 twiddle factors.
In the standard algorithm:

```
scratch[0] = tw1 * b
scratch[1] = tw2 * c
scratch[2] = tw3 * d

p0 = a + scratch[1]    p1 = a - scratch[1]
p2 = scratch[0] + scratch[2]    p3 = scratch[0] - scratch[2]

out0 = p0 + p2
out1 = p1 + j*p3     (forward FFT)
out2 = p0 - p2
out3 = p1 - j*p3     (forward FFT)
```

This is done via a **3-phase stateful protocol** (setup → compute → readback):

### Instructions

Uses **R4-type** encoding on `OpcodeCustom2` (0x5B):

| Phase | Mnemonic | funct3 | Operands | Action |
|-------|----------|--------|----------|--------|
| Setup | `BFY4_S0` | 000 | rs1=a, rs2=b, rs3=tw1 | Capture a, b, twiddle1 |
| Setup | `BFY4_S1` | 001 | rs1=c, rs2=d, rs3=tw2 | Capture c, d, twiddle2 |
| Compute | `BFY4_S2` | 010 | rs3=tw3 | Compute all 4 outputs |
| Readback | `BFY4_O0` | 011 | — | rd = output[0] |
| Readback | `BFY4_O1` | 100 | — | rd = output[1] |
| Readback | `BFY4_O2` | 101 | — | rd = output[2] |
| Readback | `BFY4_O3` | 110 | — | rd = output[3] |

### Dataflow

```
 ┌─────── BFY4_S0 ───────┐     ┌─────── BFY4_S1 ───────┐
 │  rs1 → a (packed)      │     │  rs1 → c (packed)      │
 │  rs2 → b (packed)      │     │  rs2 → d (packed)      │
 │  rs3 → tw1 or tw_idx   │     │  rs3 → tw2 or tw_idx   │
 │                        │     │                        │
 │  Latch: a_r, a_i       │     │  Latch: c_r, c_i       │
 │         b_r, b_i       │     │         d_r, d_i       │
 │         t1_r, t1_i     │     │         t2_r, t2_i     │
 └────────────────────────┘     └────────────────────────┘
               │                            │
               └────────────┬───────────────┘
                            ▼
                 ┌─────── BFY4_S2 ───────┐
                 │  rs3 → tw3 (or idx)   │
                 │                        │
                 │  Compute:              │
                 │   b' = Q15(tw1 × b)   │
                 │   c' = Q15(tw2 × c)   │
                 │   d' = Q15(tw3 × d)   │
                 │                        │
                 │   p0 = a + c'          │
                 │   p1 = a - c'          │
                 │   p2 = b' + d'         │
                 │   p3 = b' - d'         │
                 │                        │
                 │   o0 = p0 + p2         │
                 │   o1 = p1 + j·p3       │
                 │   o2 = p0 - p2         │
                 │   o3 = p1 - j·p3       │
                 │                        │
                 │  Latch: o0,o1,o2,o3    │
                 └────────────────────────┘
                            │
              ┌─────────────┼─────────────┐
              ▼             ▼             ▼             ▼
        ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐
        │ BFY4_O0  │ │ BFY4_O1  │ │ BFY4_O2  │ │ BFY4_O3  │
        │ rd = o0  │ │ rd = o1  │ │ rd = o2  │ │ rd = o3  │
        └──────────┘ └──────────┘ └──────────┘ └──────────┘
```

### Dependency Chain Optimization

The four `BFY4_O*` instructions all depend on the same register (the last `t3`
value from `BFY4_S2`), **not** on each other. This avoids a sequential read
dependency chain and allows the hardware to service them back-to-back without
pipeline stalls.

### Hardware location: `core/mult.sv` (BFY4 section)

The `BFY4_S2` computation is a **large combinational block** containing:
- 3× complex multiplications (12 multipliers 16×16→32)
- 3× Q15 rounding
- Butterfly additions/subtractions
- All in one cycle (combinational)

---

## Tier 4 — Twiddle Factor Cache

### Concept

The FFT's inner loop repeatedly loads twiddle factors from memory. For N=512,
there are up to 512 unique twiddle factors. A **hardware twiddle register file** stores
all 512 packed complex values (32 bits each) in flip-flops inside the MULT unit,
eliminating memory loads entirely.

### Instructions

| Mnemonic | funct3 | Action |
|----------|--------|--------|
| `TWLD` | 101 | `twiddle_mem[rs1] = rs2` — write one entry |
| `TWCFG` | 110 | `twiddle_rf_en = rs1[0]; bfy4_scale_en = rs1[1]` — enable/configure |

### Hardware

```
  ┌──────────────────────────────────────────────┐
  │            Twiddle Register File (in mult.sv)         │
  │                                              │
  │   twiddle_mem[0..511]  (512 × 32-bit regs)  │
  │                                              │
  │   TWLD: twiddle_mem[rs1] ← rs2              │
  │   TWCFG: cache_en ← rs1[0]                  │
  │          scale_en ← rs1[1]                  │
  │                                              │
  │   MUX: twiddle_word =                        │
  │     cache_en ? twiddle_mem[imm[8:0]]         │
  │              : imm (direct register value)   │
  │                                              │
  └──────────────────────────────────────────────┘
              │
              ▼
    Used by BFY4_S0 (tw1), BFY4_S1 (tw2), BFY4_S2 (tw3)
    and by BFY2 (twiddle operand)
```

### Software Usage

At startup (before FFT):
```c
// Load all twiddles into hardware cache
for (int i = 0; i < 512; i++) {
    uint32_t v = pack(twiddles[i]);
    TWLD(i, v);    // kissfft_twld_u32(i, v)
}
TWCFG(0b11);       // enable cache + enable scale
```

In the FFT loop, twiddle factors are passed as **indices** instead of values:
```c
uint32_t t1 = (uint32_t)tw1_idx;   // just an integer index
BFY4_S0(a, b, t1);                 // hardware reads twiddle_mem[t1]
```

---

## Tier 5 — Fused Scaling (÷4 Fusion)

### Concept

In fixed-point FFT, each radix-4 stage applies `C_FIXDIV(x, 4)` to prevent
overflow. This is normally:

```c
x.r = Q15_round(x.r * (SAMP_MAX/4));
x.i = Q15_round(x.i * (SAMP_MAX/4));
```

This costs **8 multiply-round operations** per butterfly (4 complex values × 2 components).

The fusion moves this scaling **inside** the BFY4 hardware: when `bfy4_scale_en`
is set (via `TWCFG`), the `BFY4_S0` and `BFY4_S1` setup phases automatically
divide each input by 4 before latching.

### Hardware function

```systemverilog
function automatic logic signed [15:0] q15_div4(input logic signed [15:0] x);
    logic signed [31:0] prod;
    prod = x * 32'sd8191;        // 8191 ≈ 32767/4
    q15_div4 = q15_round(prod);  // Q15 round
endfunction
```

### Dataflow

```
              BFY4_S0 (with scale_en = 1)
  ┌─────────────────────────────────────────┐
  │  rs1 → a_packed                         │
  │  rs2 → b_packed                         │
  │                                         │
  │  a_r_scaled = q15_div4(a_packed[15:0])  │
  │  a_i_scaled = q15_div4(a_packed[31:16]) │
  │  b_r_scaled = q15_div4(b_packed[15:0])  │
  │  b_i_scaled = q15_div4(b_packed[31:16]) │
  │                                         │
  │  Latch scaled values                    │
  └─────────────────────────────────────────┘
```

This eliminates 8 software instructions per butterfly iteration.

---

## Data Flow Diagrams

### Complete BFY4 with Twiddle Register File + Scale — One Butterfly Iteration

```
 Software:                          Hardware (mult.sv):
 ─────────                          ───────────────────
 
 lw a, [Fout]                       
 lw b, [Fout+m]                     
 lw c, [Fout+2m]                    
 lw d, [Fout+3m]                    
                                    
 BFY4_S0(a, b, tw1_idx)  ───────▶  Latch ÷4(a), ÷4(b)
                                    tw1 = twiddle_mem[tw1_idx]
                                    
 BFY4_S1(c, d, tw2_idx)  ───────▶  Latch ÷4(c), ÷4(d)
                                    tw2 = twiddle_mem[tw2_idx]
                                    
 BFY4_S2(tw3_idx)        ───────▶  tw3 = twiddle_mem[tw3_idx]
                                    Compute all 4 outputs
                                    ┌─────────────────────────┐
                                    │ 3× complex multiply     │
                                    │ 3× Q15 round            │
                                    │ Butterfly add/sub       │
                                    │ j-rotation for fwd FFT  │
                                    └─────────────────────────┘
                                    Latch o0, o1, o2, o3
                                    
 o0 = BFY4_O0(dep)       ◀───────  Read o0
 o1 = BFY4_O1(dep)       ◀───────  Read o1
 o2 = BFY4_O2(dep)       ◀───────  Read o2
 o3 = BFY4_O3(dep)       ◀───────  Read o3
                                    
 sw o0, [Fout]                      
 sw o1, [Fout+m]                    
 sw o2, [Fout+2m]                   
 sw o3, [Fout+3m]                   
```

### Instruction Count Comparison (per radix-4 butterfly)

| Version | Approx. instructions | Notes |
|---------|---------------------|-------|
| Baseline (scalar) | ~50-60 | 3 C_MUL + adds + stores + fixdiv |
| + Packed ops (Tier 1) | ~30-35 | Pack/unpack overhead |
| + BFY4 (Tier 3) | ~15-18 | 7 custom insns + 8 load/store |
| + Twiddle register file (Tier 4) | ~13-16 | Remove twiddle loads |
| + Scale fusion (Tier 5) | ~11-14 | Remove C_FIXDIV |

---

## Instruction Encoding Reference

### Custom-0 (OpcodeCustom0 = 0x0B) — R-type

```
 31    25  24  20  19  15  14 12  11   7  6    0
┌────────┬───────┬───────┬──────┬───────┬───────┐
│ funct7 │  rs2  │  rs1  │ f3   │  rd   │0001011│
│0000000 │       │       │      │       │ 0x0B  │
└────────┴───────┴───────┴──────┴───────┴───────┘

f3=000: CADD16    rd = rs1 + rs2  (packed complex)
f3=001: CSUB16    rd = rs1 - rs2  (packed complex)
f3=010: CMUL16    rd = Q15(rs1 × rs2) (packed complex)
f3=011: CMADD16   rd = Q15(rs1 × rs2) + rd_old (rs3=rd)
f3=100: CMSUB16   rd = rd_old - Q15(rs1 × rs2) (rs3=rd)
f3=101: TWLD      twiddle_mem[rs1] = rs2 (rd=x0)
f3=110: TWCFG     config ← rs1 (rd=x0)
```

### Custom-1 (OpcodeCustom1 = 0x2B) — R4-type

```
 31  27 26 25  24  20  19  15  14 12  11   7  6    0
┌──────┬─────┬───────┬───────┬──────┬───────┬───────┐
│ rs3  │ f2  │  rs2  │  rs1  │  f3  │  rd   │0101011│
│      │ 00  │       │       │      │       │ 0x2B  │
└──────┴─────┴───────┴───────┴──────┴───────┴───────┘

f3=000: BFY2     rd = bfly2_lo(rs1=a, rs2=b, rs3=tw) [2-cycle]
f3=001: BFY2H    rd = bfly2_hi(rs1=dep)               [1-cycle]
```

### Custom-2 (OpcodeCustom2 = 0x5B) — R4-type

```
 31  27 26 25  24  20  19  15  14 12  11   7  6    0
┌──────┬─────┬───────┬───────┬──────┬───────┬───────┐
│ rs3  │ f2  │  rs2  │  rs1  │  f3  │  rd   │1011011│
│      │ 00  │       │       │      │       │ 0x5B  │
└──────┴─────┴───────┴───────┴──────┴───────┴───────┘

f3=000: BFY4_S0  setup(rs1=a, rs2=b, rs3=tw1)  rd=x0
f3=001: BFY4_S1  setup(rs1=c, rs2=d, rs3=tw2)  rd=x0
f3=010: BFY4_S2  compute(rs3=tw3)              rd=x0
f3=011: BFY4_O0  rd = output[0]
f3=100: BFY4_O1  rd = output[1]
f3=101: BFY4_O2  rd = output[2]
f3=110: BFY4_O3  rd = output[3]
```

---

## Software Integration

### Inline Assembly Intrinsics (`_kiss_fft_guts.h`)

All custom instructions are exposed via `static inline` functions using GCC
inline assembly with the `.insn` directive:

```c
// Example: packed complex multiply
static inline uint32_t kissfft_cmul_u32(uint32_t a, uint32_t b) {
    uint32_t r;
    __asm__ volatile (
        ".insn r 0x0B, 2, 0x00, %0, %1, %2"
        : "=r"(r) : "r"(a), "r"(b)
    );
    return r;
}

// Example: BFY4 setup stage 0
static inline void kissfft_bfy4_s0(uint32_t a, uint32_t b, uint32_t t1) {
    uint32_t tmp;
    __asm__ volatile (
        ".insn r4 0x5B, 0, 0, %0, %1, %2, %3"
        : "=r"(tmp) : "r"(a), "r"(b), "r"(t1)
    );
}
```

### Compile-time Feature Flags

The software uses preprocessor guards to enable/disable each tier:

```c
#ifdef KISSFFT_USE_IM_CUSTOM     // Tier 1: packed complex ops
#ifdef KISSFFT_USE_BFY2          // Tier 2: radix-2 butterfly
#ifdef KISSFFT_USE_BFY4          // Tier 3: radix-4 butterfly
#ifdef KISSFFT_USE_TWIDDLE_RF // Tier 4: twiddle register file
#ifdef KISSFFT_USE_BFY4_SCALE    // Tier 5: fused ÷4 scaling
#ifdef KISSFFT_MANUAL_UNROLL4    // Optimization: loop unrolling
```

---

## Build Flags

In `sw/app/Makefile`, the kissfft library is compiled with:

```makefile
CFLAGS += -DKISSFFT_USE_IM_CUSTOM    \
          -DKISSFFT_USE_BFY2         \
          -DKISSFFT_USE_BFY4         \
          -DKISSFFT_MANUAL_UNROLL4   \
          -DKISSFFT_USE_TWIDDLE_RF \
          -DKISSFFT_USE_BFY4_SCALE
```

And the main application additionally gets:
```makefile
RISCV_CFLAGS += -DKISSFFT_USE_TWIDDLE_RF \
                -DKISSFFT_USE_BFY4_SCALE
```

---

## Files Modified

### Hardware (RTL)

| File | Changes |
|------|---------|
| `core/include/ariane_pkg.sv` | Added `CADD16`, `CSUB16`, `CMUL16`, `CMADD16`, `CMSUB16`, `TWLD`, `TWCFG`, `BFY2`, `BFY2H`, `BFY4_S0..S2`, `BFY4_O0..O3` to `fu_op` enum. Updated `is_rs3_gpr()` for `CMADD16`, `CMSUB16`, `BFY2`, `BFY4_S0..S2`. |
| `core/decoder.sv` | Added decode logic for `OpcodeCustom0` (ALU ops + TWLD/TWCFG), `OpcodeCustom1` (BFY2/BFY2H), `OpcodeCustom2` (BFY4_S0..O3). |
| `core/alu.sv` | Added packed complex arithmetic datapath (`CADD16`, `CSUB16`, `CMUL16`, `CMADD16`, `CMSUB16`) with Q15 rounding. |
| `core/mult.sv` | Added BFY2 2-cycle pipeline, BFY4 stateful unit (setup/compute/readback), 512-entry twiddle register file, `q15_div4` scaling function, output arbitration. |

### Software

| File | Changes |
|------|---------|
| `sw/app/fft/kissfft_lib/_kiss_fft_guts.h` | Added inline asm intrinsics for all 16 custom instructions. Added pack/unpack helpers. Redefined C_MUL, C_ADD, C_SUB, etc. to use hardware ops. |
| `sw/app/fft/kissfft_lib/kiss_fft.c` | Modified `kf_bfly4()` to use BFY4 protocol with twiddle register file indices. Added `kiss_fft_hw_load_twiddles()`. Conditional C_FIXDIV removal under `KISSFFT_USE_BFY4_SCALE`. |
| `sw/app/fft/fft_int16_main.c` | Added `kiss_fft_hw_load_twiddles(cfg)` call before FFT execution. |
| `sw/app/Makefile` | Added compile flags for all tiers. Added FORCE dependency for lib rebuild. |

---

## Q15 Fixed-Point Arithmetic Reference

All complex multiplications use Q15 (1.15) fixed-point format:

```
Value range: [-1.0, +0.999969] mapped to [-32768, +32767]

Multiply:   result = (a * b + 16384) >> 15
Divide by 4: result = (x * 8191 + 16384) >> 15    (≈ x * 0.25)
```

The `+16384` bias provides correct rounding (round-to-nearest).
