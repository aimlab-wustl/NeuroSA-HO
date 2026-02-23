`timescale 1ns/1ps

module top_solver_axi #(
    parameter AXI_DATA_WIDTH = 32,
    parameter STREAM_WIDTH   = 16,
    parameter CONFIG_WIDTH   = 4,
    parameter SAT_WIDTH      = 16,
    parameter FIFO_DEPTH     = 4
)(
    input  wire                       axi_clk,
    input  wire                       axi_reset_n,

    input  wire                       s_axis_valid,
    input  wire [AXI_DATA_WIDTH-1:0]  s_axis_data,
    input  wire                       s_axis_last,
    output wire                       s_axis_ready,

    input  wire                       m_axis_ready,
    output reg                        m_axis_valid,
    output reg  [AXI_DATA_WIDTH-1:0]  m_axis_data,
    output reg                        m_axis_last
);

  wire fifo_s_ready;
  wire fifo_m_valid;
  wire [AXI_DATA_WIDTH-1:0] fifo_m_data;
  wire fifo_m_last;

  wire fifo_m_ready = 1'b1;

  axis_data_fifo_0 u_axis_fifo (
    .s_axis_aresetn (axi_reset_n),
    .s_axis_aclk    (axi_clk),
    .s_axis_tvalid  (s_axis_valid),
    .s_axis_tready  (fifo_s_ready),
    .s_axis_tdata   (s_axis_data),
    .s_axis_tlast   (s_axis_last),
    .m_axis_tvalid  (fifo_m_valid),
    .m_axis_tready  (fifo_m_ready),
    .m_axis_tdata   (fifo_m_data),
    .m_axis_tlast   (fifo_m_last)
  );
  assign s_axis_ready = fifo_s_ready;

  wire                      s1_valid = fifo_m_valid;
  wire [AXI_DATA_WIDTH-1:0] s1_data  = fifo_m_data;
  wire                      s1_last  = fifo_m_last;

  reg                        s2_valid;
  reg  [AXI_DATA_WIDTH-1:0]  s2_data;
  reg                        s2_last;

  always @(posedge axi_clk or negedge axi_reset_n) begin
    if (!axi_reset_n) begin
      s2_valid <= 1'b0;
      s2_data  <= {AXI_DATA_WIDTH{1'b0}};
      s2_last  <= 1'b0;
    end else begin
      s2_valid <= s1_valid;
      s2_data  <= s1_data;
      s2_last  <= s1_last;
    end
  end

  wire [CONFIG_WIDTH-1:0]      cfg_in   = s2_data[STREAM_WIDTH +: CONFIG_WIDTH];
  wire [STREAM_WIDTH-1:0]      stream_in= s2_data[STREAM_WIDTH-1:0];

  wire        SOLVED;
  wire        T_READOUT;
  wire signed [SAT_WIDTH-1:0] sat_out;

  solver_readout u_readout (
    .clk        (axi_clk),
    .RESET      (~axi_reset_n),
    .ENABLE     (s2_valid),
    .LOAD_T     (cfg_in[0]),
    .LOAD_NOISE (cfg_in[1]),
    .RUN        (cfg_in[2]),
    .READ       (cfg_in[3]),
    .NOISE_IN   (stream_in),
    .T_IN       (stream_in[0]),
    .SOLVED     (SOLVED),
    .T_READOUT  (T_READOUT),
    .sat_out    (sat_out)
  );

  always @(posedge axi_clk or negedge axi_reset_n) begin
    if (!axi_reset_n) begin
      m_axis_valid <= 1'b0;
      m_axis_last  <= 1'b0;
      m_axis_data  <= {AXI_DATA_WIDTH{1'b0}};
    end else if (s2_valid && m_axis_ready) begin
      m_axis_valid <= 1'b1;
      m_axis_data  <= {{(AXI_DATA_WIDTH-SAT_WIDTH-2){1'b0}},
                       SOLVED, T_READOUT, sat_out};
      m_axis_last  <= s2_last;
    end else begin
      m_axis_valid <= 1'b0;
      m_axis_last  <= 1'b0;
    end
  end

endmodule
