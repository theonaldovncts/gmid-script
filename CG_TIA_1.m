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

% Actual bias voltages dynamically calculated
VSB_op1 = 0.2; % Target DC voltage at input node (Source M1 / Drain M2)
VDS_op1 = VDD - (id * rd) - VSB_op1; % Dynamic VDS for M1
VDS_op2 = VSB_op1;                   % VDS for M2 is the input node voltage

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

% Denormalize M2 parameters using its own absolute gm2
gm2  = gm_id2 * id;
gds2 = gds_gm2 * gm2;

% =========================================================================
% 6. Updated Small-Signal & Node Pole Calculations
% =========================================================================
RT_ideal = rd;

% Update A: Transimpedance Gain updated with gds2 loading at node X
RT_actual = ((1 + gmb_gm1) .* rd) ./ (1 + gmb_gm1 + gds_gm1 + (gds2 / gm1));
RT_dB_ohm = 20 * log10(RT_actual);

% Update B: Rin updated with parallel gds2 contribution
Rin_actual = (1 + rd * gds1) / (gm1 * (1 + gmb_gm1) + gds1 + gds2);

% Update C: Input capacitance extended to include M2 drain parasitics (Cdd2)
css1_actual = gm1 / gm_css1;
cdd2_actual = gm2 / gm_cdd2;
Cin_total   = c_pd + css1_actual + cdd2_actual;

% Pole 1: Calculated using exact Rin_actual
fp1 = 1 / (2 * pi * Rin_actual * Cin_total);

% Pole 2: Output node includes M1 output conductance (gds1) loading
Rout_actual = 1 / (1/rd + gds1); 
fp2 = 1 / (2 * pi * Rout_actual * (c_L + (gm1 / gm_cdd1)));

% Overall combined -3dB Bandwidth
f_3dB = 1 / (2 * pi * sqrt((1/(2*pi*fp1))^2 + (1/(2*pi*fp2))^2));

% =========================================================================
% 7. Advanced Noise Analysis (Input-Referred Noise Current)
% =========================================================================
k_B = 1.38e-23;
T = 300;
gamma_n = 0.6247;

% Effective Noise Bandwidth (ENBW) - 1st order approx
ENBW = f_3dB * (pi / 2);

% A. Low-Frequency (Flat-band) Noise Density (A^2/Hz)
% Dominated by Thermal Noise from Rd and M2
In2_LF = 4 * k_B * T * (gamma_n * gm2 + (1 / rd));
In_density_LF_pA = sqrt(In2_LF) * 1e12; % Convert to pA/sqrt(Hz)

% B. High-Frequency Noise Peaking Contribution at f_3dB (A^2/Hz)
% M1 channel noise amplified by total input capacitance
In2_HF_peak = 4 * k_B * T * (gamma_n / gm1) * (2 * pi * f_3dB * Cin_total)^2;
In_density_3dB_pA = sqrt(In2_LF + In2_HF_peak) * 1e12; % Density at f_3dB

% C. Total Integrated RMS Noise Current (A_rms)
% Integrated from DC to ENBW including the f^2 high-frequency term
Total_In2_int = (In2_LF * ENBW) + (4 * k_B * T * (gamma_n / gm1) * (2 * pi * Cin_total)^2 * (1/3) * ENBW^3);
In_rms_total_nA = sqrt(Total_In2_int) * 1e9; % Convert to nA_rms

% =========================================================================
% 8. Design Summary Printout
% =========================================================================
fprintf('\n================== TIA SIZING & BANDWIDTH REPORT ==================\n');
fprintf('M1 Transistor Width (W1)     : %.2f um\n', w1 * 1e6);
fprintf('M2 Transistor Width (W2)     : %.2f um\n', w2 * 1e6);
fprintf('Bias Current Consumption (ID): %.2f uA\n', id * 1e6);
fprintf('Actual Input Resistance (Rin): %.2f Ohm (Target: %.2f Ohm)\n', Rin_actual, Rin_target);
fprintf('-------------------------------------------------------------------\n');
fprintf('Ideal Transimpedance Gain RT : %.2f Ohm\n', RT_ideal);
fprintf('Actual Transimpedance Gain RT: %.2f Ohm (%.2f dB*Ohm)\n', RT_actual, RT_dB_ohm);
fprintf('-------------------------------------------------------------------\n');
fprintf('Input Node Pole (fp1)        : %.2f MHz\n', fp1 / 1e6);
fprintf('Output Node Pole (fp2)       : %.2f MHz\n', fp2 / 1e6);
fprintf('Estimated -3dB Bandwidth     : %.2f MHz\n', f_3dB / 1e6);
fprintf('-------------------------------------------------------------------\n');
fprintf('Input Noise Density (DC)     : %.2f pA/sqrt(Hz)\n', In_density_LF_pA);
fprintf('Input Noise Density (f_3dB)  : %.2f pA/sqrt(Hz)\n', In_density_3dB_pA);
fprintf('Total Integrated RMS Noise   : %.2f nA_rms\n', In_rms_total_nA);
fprintf('===================================================================\n');
