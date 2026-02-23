(* use_dsp = "no" *)
module adder_tree #(
    parameter integer IN_WIDTH  = 5,             // bit-width of each term
    parameter integer N_INPUTS  = 65,            // number of terms
    parameter integer OUT_WIDTH = IN_WIDTH + $clog2(N_INPUTS)
)(
    input  wire signed [IN_WIDTH*N_INPUTS-1:0] in_terms_flat,
    output wire signed [OUT_WIDTH-1:0]             sum_out
);

  localparam LEVELS = $clog2((N_INPUTS + 1)/2) + 1;
  wire signed [OUT_WIDTH-1:0] level_terms [0:LEVELS][0:N_INPUTS-1];

  genvar k, lvl;
  generate
    for (k = 0; k < N_INPUTS; k = k + 1) begin : lvl0
      wire signed [IN_WIDTH-1:0] term_k =
        in_terms_flat[(k+1)*IN_WIDTH-1 : k*IN_WIDTH];
      assign level_terms[0][k] = {{(OUT_WIDTH-IN_WIDTH){term_k[IN_WIDTH-1]}}, term_k};
    end

    for (lvl = 1; lvl <= LEVELS; lvl = lvl + 1) begin : tree_levels
      localparam PREV_CNT = (N_INPUTS + (1 << (lvl-1)) - 1) >> (lvl-1);
      localparam CURR_CNT = (PREV_CNT + 1) >> 1;

      for (k = 0; k < CURR_CNT; k = k + 1) begin : sum_pairs
        if (2*k+1 < PREV_CNT) begin
          assign level_terms[lvl][k] =
            level_terms[lvl-1][2*k] + level_terms[lvl-1][2*k+1];
        end else begin
          assign level_terms[lvl][k] = level_terms[lvl-1][2*k];
        end
      end
    end
  endgenerate

  assign sum_out = level_terms[LEVELS][0];

endmodule
