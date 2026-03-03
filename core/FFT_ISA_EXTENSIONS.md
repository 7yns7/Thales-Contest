# CVA6 FFT ISA Extensions — Hardware Reference

> Detailed hardware implementation guide for all custom instructions added to the
> CV32A6 core for FFT acceleration. For the high-level architecture overview, see
> [`docs/FFT_ACCELERATION.md`](../docs/FFT_ACCELERATION.md).

---

## Table of Contents

1. [Modified RTL Files](#modified-rtl-files)
2. [Instruction Decode Flow](#instruction-decode-flow)
3. [ALU Datapath — Packed Complex Ops](#alu-datapath--packed-complex-ops)
4. [MULT Unit — Butterfly & Twiddle Register File](#mult-unit--butterfly--twiddle-register-file)
5. [BFY2 Pipeline Detail](#bfy2-pipeline-detail)
6. [BFY4 Stateful Machine](#bfy4-stateful-machine)
7. [Twiddle Register File Memory Map](#twiddle-register-file-memory-map)
8. [Output Arbitration](#output-arbitration)
9. [Timing & Synthesis Considerations](#timing--synthesis-considerations)

---

## Modified RTL Files

```
core/
├── include/
│   └── ariane_pkg.sv        ← fu_op enum + is_rs3_gpr()
├── decoder.sv                ← OpcodeCustom0/1/2 decode
├── alu.sv                    ← CADD16, CSUB16, CMUL16, CMADD16, CMSUB16
└── mult.sv                   ← BFY2, BFY2H, BFY4_S0..O3, TWLD, TWCFG
```

---

## Instruction Decode Flow

### decoder.sv — Custom Opcode Routing

```
                    instruction_i[6:0]
                          │
              ┌───────────┼───────────┐
              ▼           ▼           ▼
         OpcodeCustom0  OpcodeCustom1  OpcodeCustom2
           (0x0B)        (0x2B)        (0x5B)
              │           │           │
         ┌────┴────┐   ┌──┴──┐    ┌──┴──────────┐
         │ funct3  │   │ f3  │    │    funct3    │
         ├─────────┤   ├─────┤    ├──────────────┤
         │0: CADD16│   │0:BFY2│   │0: BFY4_S0   │
         │1: CSUB16│   │1:BFY2H│  │1: BFY4_S1   │
         │2: CMUL16│   └─────┘    │2: BFY4_S2   │
         │3:CMADD16│              │3: BFY4_O0   │
         │4:CMSUB16│              │4: BFY4_O1   │
         │5: TWLD  │              │5: BFY4_O2   │
         │6: TWCFG │              │6: BFY4_O3   │
         └─────────┘              └──────────────┘
              │                         │
         FU routing:                FU routing:
         f3<5 → ALU                 ALL → MULT
         f3≥5 → MULT
```

### rs3 Handling (3-operand instructions)

Some instructions need a third source register:

| Instruction | rs3 source | imm_select |
|------------|-----------|------------|
| `CMADD16` | rd field (accumulator) | `MUX_RD_RS3` |
| `CMSUB16` | rd field (accumulator) | `MUX_RD_RS3` |
| `BFY2` | rs3 field (twiddle) | `RS3` |
| `BFY4_S0` | rs3 field (tw1 / idx) | `RS3` |
| `BFY4_S1` | rs3 field (tw2 / idx) | `RS3` |
| `BFY4_S2` | rs3 field (tw3 / idx) | `RS3` |

The `is_rs3_gpr()` function in `ariane_pkg.sv` tells the issue stage to read
the third operand from the register file:

```systemverilog
function automatic logic is_rs3_gpr(input fu_op op);
    unique case (op)
        CMADD16, CMSUB16, BFY2, BFY4_S0, BFY4_S1, BFY4_S2: return 1'b1;
        default: return 1'b0;
    endcase
endfunction
```

---

## ALU Datapath — Packed Complex Ops

### Architecture in `alu.sv`

The five packed complex operations are implemented as **purely combinational**
logic alongside the existing ALU operations:

```
                  fu_data_i.operand_a           fu_data_i.operand_b
                       │                              │
              ┌────────┴────────┐            ┌────────┴────────┐
              │ [31:16] [15:0]  │            │ [31:16] [15:0]  │
              │  a_i     a_r    │            │  b_i     b_r    │
              └──┬────────┬─────┘            └──┬────────┬─────┘
                 │        │                     │        │
                 │   ┌────┴─────────────────────┴────┐   │
                 │   │       4× 16-bit multipliers   │   │
                 │   │   rr = a_r × b_r              │   │
                 │   │   ii = a_i × b_i              │   │
                 │   │   ri = a_r × b_i              │   │
                 │   │   ir = a_i × b_r              │   │
                 │   └────────────┬──────────────────┘   │
                 │                │                       │
                 │    ┌───────────┴───────────┐           │
                 │    │ real32 = rr - ii      │           │
                 │    │ imag32 = ri + ir      │           │
                 │    │ real_q15 = Q15(real32) │           │
                 │    │ imag_q15 = Q15(imag32) │           │
                 │    └───────────┬───────────┘           │
                 │                │                       │
          ┌──────┼────────────────┼───────────────────────┤
          │      ▼                ▼                       ▼
          │  ┌───────┐    ┌──────────────┐         ┌──────────┐
          │  │ ADD16  │    │   MUL Q15    │         │  SUB16   │
          │  │a_r+b_r │    │{imag_q15,    │         │ a_r-b_r  │
          │  │a_i+b_i │    │ real_q15}    │         │ a_i-b_i  │
          │  └───┬────┘    └──────┬───────┘         └────┬─────┘
          │      │                │                      │
          │      ▼                ▼                      ▼
          │  cpx_add_packed  cpx_mul_packed       cpx_sub_packed
          │
          │               ┌─── fu_data_i.imm (rs3/accumulator)
          │               │
          │     ┌─────────┴──────────┐
          │     │   MADD / MSUB      │
          │     │ madd = mul + acc   │
          │     │ msub = acc - mul   │
          │     └────────┬───────────┘
          │              │
          │     cpx_madd_packed / cpx_msub_packed
          │
          └─────────────────────────┐
                                    ▼
                           ┌────────────────┐
                           │  Result MUX    │
                           │  (unique case) │
                           └────────┬───────┘
                                    │
                                result_o
```

### Resource Cost

- **4 multipliers** (16×16→32): shared for CMUL16, CMADD16, CMSUB16
- **2 adders** (17-bit): for CADD16/CSUB16
- **2 rounding adders** (32-bit): for Q15 rounding
- Completely combinational, no added pipeline latency

---

## MULT Unit — Butterfly & Twiddle Register File

### Top-level Organization in `mult.sv`

```
mult.sv
├── Standard multiplier (i_multiplier)
├── Standard divider (i_div / serdiv)
├── BFY2 pipeline (2-cycle)
├── BFY4 stateful unit
├── Twiddle register file (512 × 32-bit)
├── TWLD / TWCFG control
└── Output arbitration
```

### Signal Summary

```systemverilog
// Twiddle register file
logic [31:0] twiddle_mem [0:511];     // 512-entry cache
logic twiddle_rf_en_q;             // register file enabled
logic bfy4_scale_en_q;                // ÷4 scaling enabled

// Twiddle word MUX
// When cache is ON: read from twiddle_mem indexed by imm[8:0]
// When cache is OFF: pass imm directly (register value)
assign twiddle_word = twiddle_rf_en_q
    ? twiddle_mem[fu_data_i.imm[8:0]]
    : fu_data_i.imm;

// BFY2 pipeline registers
logic bfy_s0_valid_q, bfy_s1_valid_q;  // 2-stage valid
logic signed [15:0] bfy_a_r_q, bfy_a_i_q;
logic signed [31:0] bfy_mul_rr_q, bfy_mul_ii_q, bfy_mul_ri_q, bfy_mul_ir_q;

// BFY4 state registers
logic signed [15:0] bfy4_a_r_q, bfy4_a_i_q;  // from S0
logic signed [15:0] bfy4_b_r_q, bfy4_b_i_q;  // from S0
logic signed [15:0] bfy4_c_r_q, bfy4_c_i_q;  // from S1
logic signed [15:0] bfy4_d_r_q, bfy4_d_i_q;  // from S1
logic signed [15:0] bfy4_t1_r_q, bfy4_t1_i_q;
logic signed [15:0] bfy4_t2_r_q, bfy4_t2_i_q;
logic [31:0] bfy4_o0_q, bfy4_o1_q, bfy4_o2_q, bfy4_o3_q;
```

---

## BFY2 Pipeline Detail

```
 Cycle N:  BFY2 issued
           ┌──────────────────────────────────────┐
           │ bfy_issue = 1                        │
           │ Capture:                             │
           │   b_r = operand_b[15:0]              │
           │   b_i = operand_b[31:16]             │
           │   t_r = imm[15:0]  (twiddle_word)    │
           │   t_i = imm[31:16] (twiddle_word)    │
           │ Compute (registered):                │
           │   bfy_mul_rr_q ← b_r × t_r          │
           │   bfy_mul_ii_q ← b_i × t_i          │
           │   bfy_mul_ri_q ← b_r × t_i          │
           │   bfy_mul_ir_q ← b_i × t_r          │
           │ Also capture:                        │
           │   bfy_a_r_q ← operand_a[15:0]       │
           │   bfy_a_i_q ← operand_a[31:16]      │
           └──────────────────────────────────────┘
                              │
                              ▼
 Cycle N+1: bfy_s0_valid_q = 1
           ┌──────────────────────────────────────┐
           │ real32 = bfy_mul_rr_q - bfy_mul_ii_q │
           │ imag32 = bfy_mul_ri_q + bfy_mul_ir_q │
           │ real_q15 = Q15_round(real32)         │
           │ imag_q15 = Q15_round(imag32)         │
           │                                      │
           │ out0_r = a_r + real_q15              │
           │ out0_i = a_i + imag_q15              │
           │ out1_r = a_r - real_q15              │
           │ out1_i = a_i - imag_q15              │
           │                                      │
           │ bfy_result ← pack(out0)    → rd      │
           │ bfy_hi_q   ← pack(out1)   → saved    │
           └──────────────────────────────────────┘
                              │
                              ▼
 Cycle N+2: BFY2H issued (whenever software needs out1)
           ┌──────────────────────────────────────┐
           │ bfyh_result = bfy_hi_q   → rd        │
           └──────────────────────────────────────┘
```

---

## BFY4 Stateful Machine

### State Transitions

```
        ┌──────────┐  BFY4_S0   ┌──────────┐  BFY4_S1   ┌──────────┐
        │  IDLE    │───────────▶│ A,B,TW1  │───────────▶│ C,D,TW2  │
        │          │            │ latched  │            │ latched  │
        └──────────┘            └──────────┘            └──────────┘
                                                              │
                                                         BFY4_S2
                                                              │
                                                              ▼
        ┌──────────┐  BFY4_O*   ┌──────────────────────────────────┐
        │  IDLE    │◀───────────│  COMPUTE & LATCH OUTPUTS         │
        │          │            │  o0, o1, o2, o3 available        │
        └──────────┘            └──────────────────────────────────┘
```

> Note: The BFY4 unit is **not** a true state machine with explicit states.
> It uses registered values that persist between instructions. The software
> **must** follow the S0→S1→S2→O0..O3 sequence for correct operation.

### BFY4_S2 Combinational Block

This is the largest combinational path. It performs in **one cycle**:

```
 ┌─────────────────────────────────────────────────────────────────┐
 │                      BFY4_S2 COMPUTE BLOCK                      │
 │                                                                 │
 │  Inputs (from registers):                                       │
 │    a_r, a_i, b_r, b_i, c_r, c_i, d_r, d_i                     │
 │    t1_r, t1_i, t2_r, t2_i                                      │
 │  Input (from twiddle_word): t3_r, t3_i                         │
 │                                                                 │
 │  ┌─── Complex Multiply b' = tw1 × b ──────────────────────┐    │
 │  │  b_rr = b_r × t1_r    b_ii = b_i × t1_i               │    │
 │  │  b_ri = b_r × t1_i    b_ir = b_i × t1_r               │    │
 │  │  b_r32 = b_rr - b_ii  b_i32 = b_ri + b_ir             │    │
 │  │  b_r16 = Q15(b_r32)   b_i16 = Q15(b_i32)              │    │
 │  └─────────────────────────────────────────────────────────┘    │
 │                                                                 │
 │  ┌─── Complex Multiply c' = tw2 × c ──────────────────────┐    │
 │  │  (same structure as above)                              │    │
 │  │  c_r16 = Q15(c_r32)   c_i16 = Q15(c_i32)              │    │
 │  └─────────────────────────────────────────────────────────┘    │
 │                                                                 │
 │  ┌─── Complex Multiply d' = tw3 × d ──────────────────────┐    │
 │  │  (same structure as above)                              │    │
 │  │  d_r16 = Q15(d_r32)   d_i16 = Q15(d_i32)              │    │
 │  └─────────────────────────────────────────────────────────┘    │
 │                                                                 │
 │  ┌─── Butterfly Sums ─────────────────────────────────────┐    │
 │  │  p0 = a + c'      p1 = a - c'                         │    │
 │  │  p2 = b' + d'     p3 = b' - d'                        │    │
 │  └────────────────────────────────────────────────────────┘    │
 │                                                                 │
 │  ┌─── Final Outputs ──────────────────────────────────────┐    │
 │  │  o0 = p0 + p2                                          │    │
 │  │  o2 = p0 - p2                                          │    │
 │  │  o1_r = p1_r + p3_i    o1_i = p1_i - p3_r  (j·p3)    │    │
 │  │  o3_r = p1_r - p3_i    o3_i = p1_i + p3_r  (-j·p3)   │    │
 │  └────────────────────────────────────────────────────────┘    │
 │                                                                 │
 │  Resources: 12 multipliers (16×16), 6 rounding adders,         │
 │             ~16 add/sub (17-bit), 4 pack operations             │
 └─────────────────────────────────────────────────────────────────┘
```

### Optional ÷4 Scaling (Tier 5)

When `bfy4_scale_en_q = 1` (set by TWCFG bit[1]):

```
 BFY4_S0:
   a_r_latched = q15_div4(operand_a[15:0])    instead of  operand_a[15:0]
   a_i_latched = q15_div4(operand_a[31:16])   instead of  operand_a[31:16]
   b_r_latched = q15_div4(operand_b[15:0])    instead of  operand_b[15:0]
   b_i_latched = q15_div4(operand_b[31:16])   instead of  operand_b[31:16]

 BFY4_S1:
   c_r_latched = q15_div4(operand_a[15:0])
   ...same pattern...

 Where q15_div4(x) = Q15_round(x × 8191)
   8191 ≈ 32767/4 = SAMP_MAX/4
```

This adds **8 extra multipliers** (16×32) on the S0/S1 input path,
but removes 8 multiply instructions from software per butterfly.

---

## Twiddle Register File Memory Map

```
 Address    Content
 ────────   ──────────────────────────────
 0x000      twiddle[0]   = {im[15:0], re[15:0]}
 0x001      twiddle[1]   = {im[15:0], re[15:0]}
 ...
 0x1FF      twiddle[511] = {im[15:0], re[15:0]}
 ────────   ──────────────────────────────
 Total: 512 entries × 32 bits = 16,384 bits = 2 KB
 Storage: Flip-flop based (in mult.sv)
```

### TWCFG Configuration Register

```
 Bit 0: twiddle_rf_en  (0=pass-through, 1=use register file)
 Bit 1: bfy4_scale_en     (0=no scaling, 1=÷4 on BFY4 inputs)
```

---

## Output Arbitration

The MULT unit has multiple result sources. Priority (highest first):

```
 1. BFY4   (bfy4_o_valid)  ← BFY4_S*, BFY4_O*, TWLD, TWCFG
 2. BFY2   (bfy_valid)     ← BFY2 result
 3. BFY2H  (bfyh_valid)    ← BFY2H result
 4. MUL    (mul_valid)      ← Standard multiplier
 5. DIV    (div_valid)      ← Serial divider
```

```systemverilog
assign mult_trans_id_o =
    (bfy4_o_valid) ? bfy4_trans_id :
    (bfy_valid)    ? bfy_trans_id  :
    (bfyh_valid)   ? bfyh_trans_id :
    (mul_valid)    ? mul_trans_id  :
                     div_trans_id;

assign result_o =
    (bfy4_o_valid) ? bfy4_result :
    (bfy_valid)    ? bfy_result  :
    (bfyh_valid)   ? bfyh_result :
    (mul_valid)    ? mul_result  :
                     div_result;

assign mult_valid_o = div_valid | mul_valid | bfy_valid | bfyh_valid | bfy4_o_valid;
```

---

## Timing & Synthesis Considerations

### Critical Path

The `BFY4_S2` combinational block is the widest single-cycle computation:
- 12× 16-bit multiplications
- 6× 32-bit additions (Q15 rounding)
- 16× 17-bit additions (butterfly)

**Target**: 40 MHz (25 ns period) on Zybo Z7-20 (xc7z020-clg400).
**Measured**: WNS = 4.446 ns → Fmax ≈ 48.6 MHz → **no timing violation**.

### FPGA Resource Estimates

| Resource | Estimated Usage | Notes |
|----------|----------------|-------|
| LUTs | ~2000-3000 | Multipliers, adders, MUX |
| FFs | ~2500 | Twiddle register file (512×32) + pipeline regs |
| DSPs | 0-12 | Depends on synthesis inference |
| BRAM | 0 | Twiddle register file uses FF (not BRAM) |

### Area Optimization Options

- **Twiddle register file → BRAM**: Move `twiddle_mem` to a BRAM block to save ~16K FFs.
  Would require 1-cycle read latency handling.
- **BFY4 pipelining**: Split `BFY4_S2` across 2 cycles to reduce critical path.
  Would add 1 cycle latency per butterfly.
