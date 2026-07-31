clear; clc;

% =========================================================================
% 1. Load Look-Up Table (LUT)
% =========================================================================
load('180nch_1V8.mat');

% =========================================================================
% 2. Macro Component Specifications (Design Specs)
% =========================================================================
VDD  = 1.8;        % Supply voltage (V)
c_pd = 250e-15;    % Photodiode parasitic capacitance at input node (250 fF)
c_L  = 30e-15;     % Load capacitance at output node (30 fF)

% =========================================================================
% 3. Core Design Choices
% =========================================================================
Rin_target = 100;  % Target input resistance of the CG TIA stage (Ohm)
rd         = 1e3;  % Drain load resistance (1 kOhm)

% --- Degrees of Freedom for M1 ---
gm_id1     = 15;      % Target current efficiency M1 (S/A)
L1         = 0.4e-6;  % Channel length M1 (0.4 um)

% --- Degrees of Freedom for M2 (Tail Current Source) ---
gm_id2     = 10;      % Chosen in Strong Inversion to minimize width and Cdd2
L2         = 0.6e-6;  % Slightly longer channel to boost ro2 (decrease gds2)

% =========================================================================
% 4. Transistor Sizing Procedure (M1 and M2)
% =========================================================================
% Absolute parameters derived from specs
gm1 = 1 / Rin_target; % Target transconductance for M1 (~10 mS)
id  = gm1 / gm_id1;   % Required drain bias current for both M1 and M2 (~666.7 uA)

% Actual bias voltages from Cadence operating point
%VDS_op1   = 1.8-(id*rd);  % VDS for M1 (V)
%VSB_op1   = 0.067;  % VSB for M1 (V) -> This forms the drain voltage for M2!
%VDS_op2   = 0.067;  % VDS for M2 = Node X DC voltage (1.073 V)
%VSB_op2   = 0;      % M2 Source is tied to Ground (0 V)

% --- Sizing M1 ---
VGS1 = look_upVGS(nch, 'GM_ID', gm_id1, 'VDS', 0.93, 'VSB', 0.2, 'L', L1);
VG1 = VGS1 + 0.2;
jd1  = look_up(nch, 'ID_W',  'GM_ID', gm_id1, 'L', L1);
w1   = id / jd1;

% --- Sizing M2 ---
%VGS_ref = look_upVGS(nch, 'GM_ID', gm_id2, 'VGS', 'L', L2);
VGS2 = look_upVGS(nch, 'GM_ID', gm_id2, 'VDS', 0.2, 'L', L2);
jd_ref = look_up(nch, 'ID_W',  'GM_ID', gm_id2, 'VDS', VGS2, 'L', L2);
w2_ref = id / jd_ref;
jd2  = look_up(nch, 'ID_W',  'GM_ID', gm_id2, 'VDS', 0.2, 'L', L2);
w2   = id / jd2;

% =========================================================================
% 5. Intrinsic Technology Ratio Extractions from LUT
% =========================================================================
% --- M1 Parameters ---
gmb_gm1 = look_up(nch, 'GMB_GM', 'GM_ID', gm_id1, 'L', L1);
gds_gm1 = look_up(nch, 'GDS_GM', 'GM_ID', gm_id1, 'L', L1);
gm_css1 = look_up(nch, 'GM_CSS', 'GM_ID', gm_id1, 'L', L1);
gm_cdd1 = look_up(nch, 'GM_CDD', 'GM_ID', gm_id1, 'L', L1);

gds1 = gds_gm1 * gm1;

% --- M2 Parameters ---
% Extract M2 ratios using M2 target gm/ID, L2, and its distinct bias conditions
gds_gm2 = look_up(nch, 'GDS_GM', 'GM_ID', gm_id2, 'L', L2);
gm_cdd2 = look_up(nch, 'GM_CDD', 'GM_ID', gm_id2, 'L', L2);

% Denormalize M2 parameters using its own absolute gm2
gm2  = gm_id2 * id;
gds2 = gds_gm2 * gm2;

% =========================================================================
% 6. Updated Small-Signal & Node Pole Calculations (With M2 Loaded)
% =========================================================================
RT_ideal = rd;

% Update A: Transimpedance Gain updated with gds2 loading at node X
RT_actual = ((1 + gmb_gm1) .* rd) ./ (1 + gmb_gm1 + gds_gm1 + (gds2 / gm1));
RT_dB_ohm = 20 * log10(RT_actual);

% Update B: Rin updated with parallel gds2 contribution at the denominator
Rin_actual = (1 + rd * gds1) / (gm1 * (1 + gmb_gm1) + gds1 + gds2);

% Update C: Input capacitance extended to include M2 drain parasitics (Cdd2)
css1_actual = gm1 / gm_css1;
cdd2_actual = gm2 / gm_cdd2;
Cin_total   = c_pd + css1_actual + cdd2_actual;

fp1 = (1 / (2 * pi)) * (gm1 * (1 + gmb_gm1) + gds1 + gds2) / Cin_total;

% Output pole remains unchanged as M2 is attached to input node X
fp2 = (1 / (2 * pi)) / rd / (c_L + (gm1 / gm_cdd1));

% Overall combined -3dB Bandwidth
f_3dB = 1 / (2 * pi * sqrt((1/(2*pi*fp1))^2 + (1/(2*pi*fp2))^2));

%noise
k_B = 1.38e-23;
T = 300;
gamma_n = 0.6247;

In_in_square = 4*k_B*T*(gamma_n*gm2+(1/rd));
In_in_rms = sqrt(In_in_square);

% =========================================================================
% 7. Design Summary Printout
% =========================================================================
fprintf('\n================== TIA SIZING & BANDWIDTH REPORT (WITH M2) ==================\n');
fprintf('M1 Transistor Width (W1)     : %.2f um\n', w1 * 1e6);
fprintf('M2 Transistor Width (W2)     : %.2f um\n', w2 * 1e6);
fprintf('Bias Current Consumption (ID): %.2f uA\n', id * 1e6);
fprintf('Actual Input Resistance (Rin): %.2f Ohm (Target: %.2f Ohm)\n', Rin_actual, Rin_target);
fprintf('---------------------------------------------------------------------------\n');
fprintf('Ideal Transimpedance Gain RT : %.2f Ohm\n', RT_ideal);
fprintf('Actual Transimpedance Gain RT: %.2f Ohm (%.2f dB*Ohm)\n', RT_actual, RT_dB_ohm);
fprintf('---------------------------------------------------------------------------\n');
fprintf('Input Node Pole (fp1)        : %.2f MHz\n', fp1 / 1e6);
fprintf('Output Node Pole (fp2)       : %.2f MHz\n', fp2 / 1e6);
fprintf('Estimated -3dB Bandwidth     : %.2f MHz\n', f_3dB / 1e6);
fprintf('===========================================================================\n');
