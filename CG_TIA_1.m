clear; clc;

% =========================================================================
% 1. Load Look-Up Table (LUT)
% =========================================================================
load('180nch_1V8.mat');

% =========================================================================
% 2. Macro Component Specifications (Design Specs)
% =========================================================================
VDD  = 1.8;        % Supply voltage (V)
c_pd = 50e-15;     % Photodiode parasitic capacitance at input node (50 fF)
c_L  = 100e-15;    % Load capacitance at output node (100 fF)

% =========================================================================
% 3. Core Design Choices
% =========================================================================
Rin_target = 400;  % Target input resistance of the CG TIA stage (Ohm)
rd         = 1.2e3;% Drain load resistance (1.2 kOhm)

% --- Degrees of Freedom for M1 ---
gm_id1     = 15;      % Target current efficiency M1 (S/A)
L1         = 0.6e-6;  % Channel length M1 (0.6 um)

% --- Degrees of Freedom for M2 (Tail Current Source) ---
gm_id2     = 10;      % Chosen in Strong Inversion to minimize width and Cdd2
L2         = 0.8e-6;  % Slightly longer channel to boost ro2 (decrease gds2)

% =========================================================================
% 4. Transistor Sizing Procedure (M1 and M2)
% =========================================================================
% Absolute parameters derived from specs
gm1 = 1 / Rin_target; % Target transconductance for M1
id  = gm1 / gm_id1;   % Required drain bias current for both M1 and M2

% Actual bias voltages dynamically calculated
VSB_op1 = 0.2;                        % Target VSB for M1 (Assume input node at 0.2V)
VDS_op1 = VDD - (id * rd) - VSB_op1;  % Dynamic VDS for M1
VDS_op2 = VSB_op1;                    % VDS for M2 is the input node voltage

% --- Sizing M1 ---
VGS1 = look_upVGS(nch, 'GM_ID', gm_id1, 'VDS', VDS_op1, 'VSB', VSB_op1, 'L', L1);
VG1  = VGS1 + VSB_op1;
jd1  = look_up(nch, 'ID_W',  'GM_ID', gm_id1, 'L', L1);
w1   = id / jd1;

% --- Sizing M2 ---
VGS2 = look_upVGS(nch, 'GM_ID', gm_id2, 'VDS', VDS_op2, 'L', L2);
jd2  = look_up(nch, 'ID_W',  'GM_ID', gm_id2, 'VDS', VDS_op2, 'L', L2);
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
gds_gm2 = look_up(nch, 'GDS_GM', 'GM_ID', gm_id2, 'L', L2);
gm_cdd2 = look_up(nch, 'GM_CDD', 'GM_ID', gm_id2, 'L', L2);

gm2  = gm_id2 * id;
gds2 = gds_gm2 * gm2;

% =========================================================================
% 6. Updated Small-Signal & Node Pole Calculations (With M2 Loaded)
% =========================================================================
RT_ideal = rd;

% Update A: Transimpedance Gain updated with gds2 loading at node X
RT_actual = ((1 + gmb_gm1) * rd) / (1 + gmb_gm1 + gds_gm1 + (gds2 / gm1));
RT_dB_ohm = 20 * log10(RT_actual);

% Update B: Rin updated with parallel gds2 contribution at the denominator
Rin_actual = (1 + rd * gds1) / (gm1 * (1 + gmb_gm1) + gds1 + gds2);

% Update C: Input capacitance extended to include M2 drain parasitics (Cdd2)
css1_actual = gm1 / gm_css1;
cdd2_actual = gm2 / gm_cdd2;
Cin_total   = c_pd + css1_actual + cdd2_actual;

% Polarity calculation improved by using Rin_actual directly
fp1 = 1 / (2 * pi * Rin_actual * Cin_total);
fp2 = 1 / (2 * pi * rd * (c_L + (gm1 / gm_cdd1)));

% Overall combined -3dB Bandwidth
f_3dB = 1 / (2 * pi * sqrt((1/(2*pi*fp1))^2 + (1/(2*pi*fp2))^2));

% =========================================================================
% 7. Advanced Noise Calculations
% =========================================================================
k_B = 1.38e-23;
T = 300;
gamma_n = 0.6247;

% Low-Frequency (DC) Noise contribution from M2 and rd
In_in_square_LF = 4 * k_B * T * (gamma_n * gm2 + (1/rd));
In_in_rms_LF_pA = sqrt(In_in_square_LF) * 1e12; % Convert to pA/sqrt(Hz)

% High-Frequency noise peaking contribution from M1 at Bandwidth Edge
In_in_square_HF = 4 * k_B * T * (gamma_n / gm1) * (2 * pi * f_3dB * Cin_total)^2;
In_in_rms_3dB_pA = sqrt(In_in_square_LF + In_in_square_HF) * 1e12;

% =========================================================================
% 8. Design Summary Printout
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
fprintf('---------------------------------------------------------------------------\n');
fprintf('Estimated Noise Density (DC) : %.2f pA/sqrt(Hz)\n', In_in_rms_LF_pA);
fprintf('Estimated Noise Density(3dB) : %.2f pA/sqrt(Hz)\n', In_in_rms_3dB_pA);
fprintf('===========================================================================\n');
