

module mult
  import ariane_pkg::*;
#(
    parameter config_pkg::cva6_cfg_t CVA6Cfg = config_pkg::cva6_cfg_empty,
    parameter type fu_data_t = logic
) (
    // Subsystem Clock - SUBSYSTEM
    input  logic                                 clk_i,
    // Asynchronous reset active low - SUBSYSTEM
    input  logic                                 rst_ni,
    // Flush - CONTROLLER
    input  logic                                 flush_i,
    // FU data needed to execute instruction - ISSUE_STAGE
    input  fu_data_t                             fu_data_i,
    // Mult instruction is valid - ISSUE_STAGE
    input  logic                                 mult_valid_i,
    // Mult result - ISSUE_STAGE
    output logic     [         CVA6Cfg.XLEN-1:0] result_o,
    // Mult result is valid - ISSUE_STAGE
    output logic                                 mult_valid_o,
    // Mutl is ready - ISSUE_STAGE
    output logic                                 mult_ready_o,
    // Mult transaction ID - ISSUE_STAGE
    output logic     [CVA6Cfg.TRANS_ID_BITS-1:0] mult_trans_id_o
);
  function automatic logic signed [15:0] q15_round(input logic signed [31:0] x);
    logic signed [31:0] tmp;
    begin
      tmp = x + 32'sd16384;
      q15_round = $signed(tmp) >>> 15;
    end
  endfunction

  function automatic logic signed [15:0] q15_div4(input logic signed [15:0] x);
    logic signed [31:0] prod;
    begin
      prod = x * 32'sd8191;
      q15_div4 = q15_round(prod);
    end
  endfunction

  localparam int unsigned TWIDDLE_DEPTH = 512;
  localparam int unsigned TWIDDLE_ADDR_BITS = $clog2(TWIDDLE_DEPTH);
  logic [31:0] twiddle_mem [0:TWIDDLE_DEPTH-1];
  logic twiddle_rf_en_q;
  logic bfy4_scale_en_q;
  logic [CVA6Cfg.XLEN-1:0] twiddle_word;
  assign twiddle_word = (twiddle_rf_en_q) ? twiddle_mem[fu_data_i.imm[TWIDDLE_ADDR_BITS-1:0]] : fu_data_i.imm;

  logic mul_valid;
  logic div_valid;
  logic bfy_valid;
  logic bfyh_valid;
  logic div_ready_i;  // receiver of division result is able to accept the result
  logic mult_ready_div;
  logic [CVA6Cfg.TRANS_ID_BITS-1:0] mul_trans_id;
  logic [CVA6Cfg.TRANS_ID_BITS-1:0] div_trans_id;
  logic [CVA6Cfg.TRANS_ID_BITS-1:0] bfy_trans_id;
  logic [CVA6Cfg.TRANS_ID_BITS-1:0] bfyh_trans_id;
  logic [CVA6Cfg.XLEN-1:0] mul_result;
  logic [CVA6Cfg.XLEN-1:0] div_result;
  logic [CVA6Cfg.XLEN-1:0] bfy_result;
  logic [CVA6Cfg.XLEN-1:0] bfyh_result;
  logic bfy4_o_valid;
  logic [CVA6Cfg.TRANS_ID_BITS-1:0] bfy4_trans_id;
  logic [CVA6Cfg.XLEN-1:0] bfy4_result;

  logic div_valid_op;
  logic mul_valid_op;
  // Input Arbitration

  assign mul_valid_op = ~flush_i && mult_valid_i && (fu_data_i.operation inside { MUL, MULH, MULHU, MULHSU, MULW, CLMUL, CLMULH, CLMULR });

  assign div_valid_op = ~flush_i && mult_valid_i && (fu_data_i.operation inside { DIV, DIVU, DIVW, DIVUW, REM, REMU, REMW, REMUW });

  logic bfy_issue;
  logic bfyh_issue;
  assign bfy_issue  = ~flush_i && mult_valid_i && (fu_data_i.operation == BFY2);
  assign bfyh_issue = ~flush_i && mult_valid_i && (fu_data_i.operation == BFY2H);

  logic twld_issue;
  logic twcfg_issue;
  assign twld_issue  = ~flush_i && mult_valid_i && (fu_data_i.operation == TWLD);
  assign twcfg_issue = ~flush_i && mult_valid_i && (fu_data_i.operation == TWCFG);

  // ---------------------
  // Output Arbitration
  // ---------------------
  // we give precedence to multiplication as the divider supports stalling and the multiplier is
  // just a dumb pipelined multiplier
  assign div_ready_i = (mul_valid | bfy_valid | bfyh_valid | bfy4_o_valid) ? 1'b0 : 1'b1;
  assign mult_trans_id_o = (bfy4_o_valid) ? bfy4_trans_id : (bfy_valid) ? bfy_trans_id : (bfyh_valid) ? bfyh_trans_id : (mul_valid) ? mul_trans_id : div_trans_id;
  assign result_o = (bfy4_o_valid) ? bfy4_result : (bfy_valid) ? bfy_result : (bfyh_valid) ? bfyh_result : (mul_valid) ? mul_result : div_result;
  assign mult_valid_o = div_valid | mul_valid | bfy_valid | bfyh_valid | bfy4_o_valid;
  // mult_ready_o = division as the multiplication will unconditionally be ready to accept new requests

  // ---------------------
  // Multiplication
  // ---------------------
  multiplier #(
      .CVA6Cfg(CVA6Cfg)
  ) i_multiplier (
      .clk_i,
      .rst_ni,
      .trans_id_i     (fu_data_i.trans_id),
      .operation_i    (fu_data_i.operation),
      .operand_a_i    (fu_data_i.operand_a),
      .operand_b_i    (fu_data_i.operand_b),
      .result_o       (mul_result),
      .mult_valid_i   (mul_valid_op),
      .mult_valid_o   (mul_valid),
      .mult_trans_id_o(mul_trans_id)
  );

  // ---------------------
  // Butterfly unit (BFY2)
  // ---------------------
  logic bfy_s0_valid_q, bfy_s1_valid_q;
  logic [CVA6Cfg.TRANS_ID_BITS-1:0] bfy_trans_id_q;
  logic signed [15:0] bfy_a_r_q, bfy_a_i_q;
  logic signed [31:0] bfy_mul_rr_q, bfy_mul_ii_q, bfy_mul_ri_q, bfy_mul_ir_q;
  logic [CVA6Cfg.XLEN-1:0] bfy_hi_q;
  logic bfy_hi_valid_q;

  logic bfyh_pending_q;
  logic [CVA6Cfg.TRANS_ID_BITS-1:0] bfyh_trans_id_q;

  always_ff @(posedge clk_i or negedge rst_ni) begin
    if (~rst_ni) begin
        twiddle_rf_en_q <= 1'b0;
        bfy4_scale_en_q <= 1'b0;
      bfy_s0_valid_q <= 1'b0;
      bfy_s1_valid_q <= 1'b0;
      bfy_trans_id_q <= '0;
      bfy_a_r_q <= '0;
      bfy_a_i_q <= '0;
      bfy_mul_rr_q <= '0;
      bfy_mul_ii_q <= '0;
      bfy_mul_ri_q <= '0;
      bfy_mul_ir_q <= '0;
      bfy_hi_q <= '0;
      bfy_hi_valid_q <= 1'b0;
      bfyh_pending_q <= 1'b0;
      bfyh_trans_id_q <= '0;
    end else if (flush_i) begin
      bfy_s0_valid_q <= 1'b0;
      bfy_s1_valid_q <= 1'b0;
      bfy_hi_valid_q <= 1'b0;
      bfyh_pending_q <= 1'b0;
    end else begin
      if (twcfg_issue) begin
        twiddle_rf_en_q <= fu_data_i.operand_a[0];
        bfy4_scale_en_q <= fu_data_i.operand_a[1];
      end
      if (twld_issue) begin
        twiddle_mem[fu_data_i.operand_a[TWIDDLE_ADDR_BITS-1:0]] <= fu_data_i.operand_b;
      end
      bfy_s0_valid_q <= bfy_issue;
      bfy_s1_valid_q <= bfy_s0_valid_q;

      if (bfy_issue) begin
        logic signed [15:0] b_r, b_i, t_r, t_i;
        b_r = fu_data_i.operand_b[15:0];
        b_i = fu_data_i.operand_b[31:16];
        t_r = fu_data_i.imm[15:0];
        t_i = fu_data_i.imm[31:16];
        bfy_mul_rr_q <= b_r * t_r;
        bfy_mul_ii_q <= b_i * t_i;
        bfy_mul_ri_q <= b_r * t_i;
        bfy_mul_ir_q <= b_i * t_r;
        bfy_a_r_q <= fu_data_i.operand_a[15:0];
        bfy_a_i_q <= fu_data_i.operand_a[31:16];
        bfy_trans_id_q <= fu_data_i.trans_id;
      end

      if (bfy_s0_valid_q) begin
        logic signed [31:0] real32;
        logic signed [31:0] imag32;
        logic signed [15:0] real_q15;
        logic signed [15:0] imag_q15;
        logic signed [16:0] out0_r;
        logic signed [16:0] out0_i;
        logic signed [16:0] out1_r;
        logic signed [16:0] out1_i;
        real32 = bfy_mul_rr_q - bfy_mul_ii_q;
        imag32 = bfy_mul_ri_q + bfy_mul_ir_q;
        real_q15 = q15_round(real32);
        imag_q15 = q15_round(imag32);
        out0_r = bfy_a_r_q + real_q15;
        out0_i = bfy_a_i_q + imag_q15;
        out1_r = bfy_a_r_q - real_q15;
        out1_i = bfy_a_i_q - imag_q15;
        bfy_result <= {out0_i[15:0], out0_r[15:0]};
        bfy_hi_q <= {out1_i[15:0], out1_r[15:0]};
        bfy_hi_valid_q <= 1'b1;
      end

      if (bfyh_issue) begin
        bfyh_pending_q <= 1'b1;
        bfyh_trans_id_q <= fu_data_i.trans_id;
      end else begin
        bfyh_pending_q <= 1'b0;
      end
    end
  end

  assign bfy_valid = bfy_s1_valid_q;
  assign bfy_trans_id = bfy_trans_id_q;

  assign bfyh_valid = bfyh_pending_q;
  assign bfyh_result = bfy_hi_q;
  assign bfyh_trans_id = bfyh_trans_id_q;

  // ---------------------
  // Butterfly-4 unit (stateful)
  // ---------------------
  logic bfy4_s0, bfy4_s1, bfy4_s2;
  logic bfy4_o0, bfy4_o1, bfy4_o2, bfy4_o3;
  logic bfy4_s_valid;

  assign bfy4_s0 = ~flush_i && mult_valid_i && (fu_data_i.operation == BFY4_S0);
  assign bfy4_s1 = ~flush_i && mult_valid_i && (fu_data_i.operation == BFY4_S1);
  assign bfy4_s2 = ~flush_i && mult_valid_i && (fu_data_i.operation == BFY4_S2);
  assign bfy4_o0 = ~flush_i && mult_valid_i && (fu_data_i.operation == BFY4_O0);
  assign bfy4_o1 = ~flush_i && mult_valid_i && (fu_data_i.operation == BFY4_O1);
  assign bfy4_o2 = ~flush_i && mult_valid_i && (fu_data_i.operation == BFY4_O2);
  assign bfy4_o3 = ~flush_i && mult_valid_i && (fu_data_i.operation == BFY4_O3);
  assign bfy4_s_valid = bfy4_s0 | bfy4_s1 | bfy4_s2;

  logic signed [15:0] bfy4_a_r_q, bfy4_a_i_q;
  logic signed [15:0] bfy4_b_r_q, bfy4_b_i_q;
  logic signed [15:0] bfy4_c_r_q, bfy4_c_i_q;
  logic signed [15:0] bfy4_d_r_q, bfy4_d_i_q;
  logic signed [15:0] bfy4_t1_r_q, bfy4_t1_i_q;
  logic signed [15:0] bfy4_t2_r_q, bfy4_t2_i_q;
  logic signed [15:0] bfy4_t3_r_q, bfy4_t3_i_q;

  logic [CVA6Cfg.XLEN-1:0] bfy4_o0_q, bfy4_o1_q, bfy4_o2_q, bfy4_o3_q;
  logic [CVA6Cfg.XLEN-1:0] bfy4_o0_d, bfy4_o1_d, bfy4_o2_d, bfy4_o3_d;
  logic signed [15:0] bfy4_a_r_s, bfy4_a_i_s;
  logic signed [15:0] bfy4_b_r_s, bfy4_b_i_s;
  logic signed [15:0] bfy4_c_r_s, bfy4_c_i_s;
  logic signed [15:0] bfy4_d_r_s, bfy4_d_i_s;

  always_comb begin
    bfy4_a_r_s = fu_data_i.operand_a[15:0];
    bfy4_a_i_s = fu_data_i.operand_a[31:16];
    bfy4_b_r_s = fu_data_i.operand_b[15:0];
    bfy4_b_i_s = fu_data_i.operand_b[31:16];
    bfy4_c_r_s = fu_data_i.operand_a[15:0];
    bfy4_c_i_s = fu_data_i.operand_a[31:16];
    bfy4_d_r_s = fu_data_i.operand_b[15:0];
    bfy4_d_i_s = fu_data_i.operand_b[31:16];

    if (bfy4_scale_en_q) begin
      bfy4_a_r_s = q15_div4(fu_data_i.operand_a[15:0]);
      bfy4_a_i_s = q15_div4(fu_data_i.operand_a[31:16]);
      bfy4_b_r_s = q15_div4(fu_data_i.operand_b[15:0]);
      bfy4_b_i_s = q15_div4(fu_data_i.operand_b[31:16]);
      bfy4_c_r_s = q15_div4(fu_data_i.operand_a[15:0]);
      bfy4_c_i_s = q15_div4(fu_data_i.operand_a[31:16]);
      bfy4_d_r_s = q15_div4(fu_data_i.operand_b[15:0]);
      bfy4_d_i_s = q15_div4(fu_data_i.operand_b[31:16]);
    end
  end

  always_comb begin : bfy4_compute
    bfy4_o0_d = bfy4_o0_q;
    bfy4_o1_d = bfy4_o1_q;
    bfy4_o2_d = bfy4_o2_q;
    bfy4_o3_d = bfy4_o3_q;

    if (bfy4_s2) begin
      logic signed [15:0] t3_r, t3_i;
      logic signed [31:0] b_rr, b_ii, b_ri, b_ir;
      logic signed [31:0] c_rr, c_ii, c_ri, c_ir;
      logic signed [31:0] d_rr, d_ii, d_ri, d_ir;
      logic signed [31:0] b_r32, b_i32;
      logic signed [31:0] c_r32, c_i32;
      logic signed [31:0] d_r32, d_i32;
      logic signed [15:0] b_r16, b_i16;
      logic signed [15:0] c_r16, c_i16;
      logic signed [15:0] d_r16, d_i16;
      logic signed [16:0] p0_r, p0_i;
      logic signed [16:0] p1_r, p1_i;
      logic signed [16:0] p2_r, p2_i;
      logic signed [16:0] p3_r, p3_i;
      logic signed [16:0] o0_r, o0_i;
      logic signed [16:0] o1_r, o1_i;
      logic signed [16:0] o2_r, o2_i;
      logic signed [16:0] o3_r, o3_i;

      t3_r = twiddle_word[15:0];
      t3_i = twiddle_word[31:16];

      b_rr = bfy4_b_r_q * bfy4_t1_r_q;
      b_ii = bfy4_b_i_q * bfy4_t1_i_q;
      b_ri = bfy4_b_r_q * bfy4_t1_i_q;
      b_ir = bfy4_b_i_q * bfy4_t1_r_q;
      b_r32 = b_rr - b_ii;
      b_i32 = b_ri + b_ir;
      b_r16 = q15_round(b_r32);
      b_i16 = q15_round(b_i32);

      c_rr = bfy4_c_r_q * bfy4_t2_r_q;
      c_ii = bfy4_c_i_q * bfy4_t2_i_q;
      c_ri = bfy4_c_r_q * bfy4_t2_i_q;
      c_ir = bfy4_c_i_q * bfy4_t2_r_q;
      c_r32 = c_rr - c_ii;
      c_i32 = c_ri + c_ir;
      c_r16 = q15_round(c_r32);
      c_i16 = q15_round(c_i32);

      d_rr = bfy4_d_r_q * t3_r;
      d_ii = bfy4_d_i_q * t3_i;
      d_ri = bfy4_d_r_q * t3_i;
      d_ir = bfy4_d_i_q * t3_r;
      d_r32 = d_rr - d_ii;
      d_i32 = d_ri + d_ir;
      d_r16 = q15_round(d_r32);
      d_i16 = q15_round(d_i32);

      p0_r = bfy4_a_r_q + c_r16;
      p0_i = bfy4_a_i_q + c_i16;
      p1_r = bfy4_a_r_q - c_r16;
      p1_i = bfy4_a_i_q - c_i16;
      p2_r = b_r16 + d_r16;
      p2_i = b_i16 + d_i16;
      p3_r = b_r16 - d_r16;
      p3_i = b_i16 - d_i16;

      o0_r = p0_r + p2_r;
      o0_i = p0_i + p2_i;
      o2_r = p0_r - p2_r;
      o2_i = p0_i - p2_i;

      o1_r = p1_r + p3_i;
      o1_i = p1_i - p3_r;
      o3_r = p1_r - p3_i;
      o3_i = p1_i + p3_r;

      bfy4_o0_d = {o0_i[15:0], o0_r[15:0]};
      bfy4_o1_d = {o1_i[15:0], o1_r[15:0]};
      bfy4_o2_d = {o2_i[15:0], o2_r[15:0]};
      bfy4_o3_d = {o3_i[15:0], o3_r[15:0]};
    end
  end

  always_ff @(posedge clk_i or negedge rst_ni) begin
    if (~rst_ni) begin
      bfy4_a_r_q <= '0; bfy4_a_i_q <= '0;
      bfy4_b_r_q <= '0; bfy4_b_i_q <= '0;
      bfy4_c_r_q <= '0; bfy4_c_i_q <= '0;
      bfy4_d_r_q <= '0; bfy4_d_i_q <= '0;
      bfy4_t1_r_q <= '0; bfy4_t1_i_q <= '0;
      bfy4_t2_r_q <= '0; bfy4_t2_i_q <= '0;
      bfy4_t3_r_q <= '0; bfy4_t3_i_q <= '0;
      bfy4_o0_q <= '0; bfy4_o1_q <= '0; bfy4_o2_q <= '0; bfy4_o3_q <= '0;
      bfy4_o_valid <= 1'b0;
      bfy4_trans_id <= '0;
      bfy4_result <= '0;
    end else if (flush_i) begin
      bfy4_o_valid <= 1'b0;
    end else begin
      bfy4_o_valid <= 1'b0;

      if (bfy4_s0) begin
        bfy4_a_r_q <= bfy4_a_r_s;
        bfy4_a_i_q <= bfy4_a_i_s;
        bfy4_b_r_q <= bfy4_b_r_s;
        bfy4_b_i_q <= bfy4_b_i_s;
        bfy4_t1_r_q <= twiddle_word[15:0];
        bfy4_t1_i_q <= twiddle_word[31:16];
      end

      if (bfy4_s1) begin
        bfy4_c_r_q <= bfy4_c_r_s;
        bfy4_c_i_q <= bfy4_c_i_s;
        bfy4_d_r_q <= bfy4_d_r_s;
        bfy4_d_i_q <= bfy4_d_i_s;
        bfy4_t2_r_q <= twiddle_word[15:0];
        bfy4_t2_i_q <= twiddle_word[31:16];
      end

      if (bfy4_s2) begin
        bfy4_t3_r_q <= twiddle_word[15:0];
        bfy4_t3_i_q <= twiddle_word[31:16];
        bfy4_o0_q <= bfy4_o0_d;
        bfy4_o1_q <= bfy4_o1_d;
        bfy4_o2_q <= bfy4_o2_d;
        bfy4_o3_q <= bfy4_o3_d;
      end

      if (bfy4_s_valid || bfy4_o0 || bfy4_o1 || bfy4_o2 || bfy4_o3 || twld_issue || twcfg_issue) begin
        bfy4_o_valid <= 1'b1;
        bfy4_trans_id <= fu_data_i.trans_id;
        if (bfy4_s_valid || twld_issue || twcfg_issue) begin
          bfy4_result <= '0;
        end else begin
          unique case (1'b1)
            bfy4_o0: bfy4_result <= bfy4_o0_q;
            bfy4_o1: bfy4_result <= bfy4_o1_q;
            bfy4_o2: bfy4_result <= bfy4_o2_q;
            bfy4_o3: bfy4_result <= bfy4_o3_q;
            default: bfy4_result <= '0;
          endcase
        end
      end
    end
  end

  // Ready logic: keep divider as the global ready signal (baseline behavior)
  assign mult_ready_o = mult_ready_div;

  // ---------------------
  // Division
  // ---------------------
  logic [CVA6Cfg.XLEN-1:0]
      operand_b,
      operand_a;  // input operands after input MUX (input silencing, word operations or full inputs)
  logic [CVA6Cfg.XLEN-1:0] result;  // result before result mux

  logic                    div_signed;  // signed or unsigned division
  logic                    rem;  // is it a reminder (or not a reminder e.g.: a division)
  logic word_op_d, word_op_q;  // save whether the operation was signed or not

  // is this a signed op?
  assign div_signed = fu_data_i.operation inside {DIV, DIVW, REM, REMW};
  // is this a modulo?
  assign rem        = fu_data_i.operation inside {REM, REMU, REMW, REMUW};

  // prepare the input operands and control divider
  always_comb begin
    // silence the inputs
    operand_a = '0;
    operand_b = '0;
    // control signals
    word_op_d = word_op_q;

    // we've go a new division operation
    if (mult_valid_i && fu_data_i.operation inside {DIV, DIVU, DIVW, DIVUW, REM, REMU, REMW, REMUW}) begin
      // is this a word operation?
      if (CVA6Cfg.IS_XLEN64 && (fu_data_i.operation == DIVW || fu_data_i.operation == DIVUW || fu_data_i.operation == REMW || fu_data_i.operation == REMUW)) begin
        // yes so check if we should sign extend this is only done for a signed operation
        if (div_signed) begin
          operand_a = sext32to64(fu_data_i.operand_a[31:0]);
          operand_b = sext32to64(fu_data_i.operand_b[31:0]);
        end else begin
          operand_a = fu_data_i.operand_a[31:0];
          operand_b = fu_data_i.operand_b[31:0];
        end

        // save whether we want sign extend the result or not, this is done for all word operations
        word_op_d = 1'b1;
      end else begin
        // regular op
        operand_a = fu_data_i.operand_a;
        operand_b = fu_data_i.operand_b;
        word_op_d = 1'b0;
      end
    end
  end

  // ---------------------
  // Serial Divider
  // ---------------------
  serdiv #(
      .CVA6Cfg(CVA6Cfg),
      .WIDTH  (CVA6Cfg.XLEN)
  ) i_div (
      .clk_i    (clk_i),
      .rst_ni   (rst_ni),
      .id_i     (fu_data_i.trans_id),
      .op_a_i   (operand_a),
      .op_b_i   (operand_b),
      .opcode_i ({rem, div_signed}),   // 00: udiv, 10: urem, 01: div, 11: rem
      .in_vld_i (div_valid_op),
      .in_rdy_o (mult_ready_div),
      .flush_i  (flush_i),
      .out_vld_o(div_valid),
      .out_rdy_i(div_ready_i),
      .id_o     (div_trans_id),
      .res_o    (result)
  );

  // Result multiplexer
  // if it was a signed word operation the bit will be set and the result will be sign extended accordingly
  assign div_result = (CVA6Cfg.IS_XLEN64 && word_op_q) ? sext32to64(result) : result;

  // ---------------------
  // Registers
  // ---------------------
  always_ff @(posedge clk_i or negedge rst_ni) begin
    if (~rst_ni) begin
      word_op_q <= '0;
    end else begin
      word_op_q <= word_op_d;
    end
  end
endmodule
