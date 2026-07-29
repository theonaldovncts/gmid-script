%% =========================================================================
%  SIMULATION & DESIGN OPTIMIZATION OF DIFFERENTIAL TIA / AMPLIFIER
%  =========================================================================
clear; clc; close all;

% 1. Load Transistor LUT Database
load('180nch_1V8.mat');

% 2. Circuit & Load Specifications
A_v0      = 4.6;         % Target Open-Loop Gain (V/V)
C_L_ext   = 50e-15;      % External Load Capacitance (50 fF)
C_in_ext  = 50e-15;      % External Input Capacitance / Photodiode (50 fF)
R_s       = 0;           % Source Resistance (0 Ohm = Ideal Current Source)

% 3. Transistor Selection Parameters & Margin
L_channel   = 0.2e-6;    % Transistor Length (200 nm)
gm_ID_target = 9;        % Initial gm/ID selection
% C_t_sim     = 71.48e-15; % Total capacitance from simulation
% margin_Ct   = 0.35;
% C_t_choice  = C_t_sim * (1 + margin_Ct);
% R_f         = Rf_Ct / C_t_choice;
R_f         = 5400;

% 4. Target Bandwidth & Gain Calculation
data_rate   = 2e9;       % 2 Gbps
bw_min      = 0.7 * data_rate; % Min Bandwidth (1.4 GHz)
Rf_Ct       = sqrt(2 * A_v0 * (A_v0 + 1)) / (2 * pi * bw_min);
bw_req      = (2 * A_v0) / (Rf_Ct * pi * 2);

margin_bw   = 0.1;
target_bw   = bw_req * (1 + margin_bw); % Final Target Bandwidth

GBW         = A_v0 * target_bw;
A_v0_req    = sqrt(2) * GBW / (bw_min * 2);

% 5. Parameter Sweep Ranges
gm_ID_sweep = 8 : 0.1 : 15;        % Sweep gm/ID from 8 to 15 V^-1
R_d_sweep   = 1000 : 100 : 1500;    % Sweep Drain Resistance (1000 Ohm - 1500 Ohm)

% Check fT from LUT
fT_sweep    = look_up(nch, 'GM_CGG', 'GM_ID', gm_ID_sweep, 'L', 180e-9:10e-9:200e-9) / (2*pi);
fT_selected = interp1(gm_ID_sweep, fT_sweep(3,:), gm_ID_target);

gm_ID_sweep = gm_ID_sweep';        % Convert to column matrix

% Pre-allocate Result Matrices
mat_ID   = zeros(length(gm_ID_sweep), length(R_d_sweep));
mat_BW   = zeros(length(gm_ID_sweep), length(R_d_sweep));
mat_Zeta = zeros(length(gm_ID_sweep), length(R_d_sweep));

%% =========================================================================
%  PARAMETER SWEEP LOOP (R_d & gm/ID)
%  =========================================================================
for j = 1:length(R_d_sweep)
    R_d = R_d_sweep(j);
    
    % Extract Intrinsic Gain (gm * rds) from LUT
    gm_rds = look_up(nch, 'GM_GDS', 'GM_ID', gm_ID_sweep, 'L', L_channel);
    
    % Calculate gm, ID, and Transistor Width W
    gm    = (1 / R_d) * ((1 / A_v0) - (1 ./ gm_rds)) .^ -1;
    I_D   = gm ./ gm_ID_sweep;
    W_dev = I_D ./ look_up(nch, 'ID_W', 'GM_ID', gm_ID_sweep, 'L', L_channel);
    
    % Calculate Parasitic Capacitances
    C_gs = W_dev .* look_up(nch, 'CGS_W', 'GM_ID', gm_ID_sweep, 'L', L_channel);
    C_gd = W_dev .* look_up(nch, 'CGD_W', 'GM_ID', gm_ID_sweep, 'L', L_channel);
    C_dd = W_dev .* look_up(nch, 'CDD_W', 'GM_ID', gm_ID_sweep, 'L', L_channel);
    C_db = C_dd - abs(C_gd);
    
    % Total Input & Output Node Capacitances
    C_L_tot  = C_L_ext + C_db + abs(C_gd);
    C_in_tot = C_in_ext + abs(C_gs) + abs(C_gd) * (1 + A_v0);
    
    % Combined Output Resistance
    r_ds    = gm_rds ./ gm;
    R_d_tot = (r_ds .* R_d) ./ (r_ds + R_d);
    
    % Bandwidth Calculation Based on R_s Case
    if R_s == 0
        % Ideal Case: Single Dominant Pole at Output Node
        tau_out       = R_d_tot .* C_L_tot;
        mat_BW(:, j) = 1 ./ (2 * pi * tau_out);
    else
        % Non-Ideal Case (Rs > 0): 2-Pole System
        tau1 = R_s * C_in_tot + R_d_tot * (C_L_tot + abs(C_gd));
        tau2 = R_s * R_d_tot * (C_in_tot .* C_L_tot + C_in_tot .* abs(C_gd) + C_L_tot .* abs(C_gd));
        
        % Pole locations
        p1 = 1 ./ (2 * pi * tau1);
        p2 = tau1 ./ (tau2 * 2 * pi);
        
        % -3dB Bandwidth Estimation for 2-Pole System (Dominant Pole / Equivalent)
        mat_BW(:, j) = 1 ./ (2 * pi * sqrt(tau1.^2 - 2*tau2));
    end
    
    mat_ID(:, j) = I_D;
    
    C_in_arr = C_in_ext + abs(C_gs);
    C_f_arr  = abs(C_gd);
    C_L_arr  = C_db + C_L_ext;
    
    num_z = R_f .* (1 + gm .* R_d) .* C_f_arr + R_d .* C_L_arr + (R_f + R_d) .* C_in_arr;
    den_z = 2 .* sqrt(R_f .* R_d .* (1 + gm .* R_d) .* (C_f_arr .* C_L_arr + C_L_arr .* C_in_arr + C_in_arr .* C_f_arr));
    
    mat_Zeta(:, j) = num_z ./ den_z;
end

% Brief Status Message to Console
if R_s == 0
    fprintf('=== Analysis Complete: Ideal Case (Rs = 0 Ohm) ===\n');
else
    fprintf('=== Analysis Complete: Non-Ideal Case (Rs = %d Ohm) ===\n', R_s);
end

%% =========================================================================
%  SWEEP CURVE PLOTTING
%  =========================================================================
figure;
plot(mat_ID / 1e-6, mat_BW / 1e6, 'LineStyle', '-');
xlabel('Bias Current I_D (\muA)', 'FontSize', 14);
ylabel('Bandwidth (MHz)', 'FontSize', 14);
title('Bandwidth Sweep vs Bias Current');
set(gca, 'FontSize', 12);
grid on;

%% =========================================================================
%  AUTOMATIC EXTRACTION OF OPTIMAL POINT (Lowest Power / Min I_D)
%  =========================================================================
target_zeta = 0.707;
tol_zeta    = 0.05;

valid_indices = find(mat_BW >= target_bw & mat_Zeta >= (target_zeta - tol_zeta) & mat_Zeta <= (target_zeta + tol_zeta));

if isempty(valid_indices)
    fprintf('\n [Warning] No point with Zeta ~0.707 satisfies the BW target. \n');
    fprintf('Searching for points achieving target BW with Zeta closest to ~0.707...\n');
    
    valid_bw_indices = find(mat_BW >= target_bw);
    if isempty(valid_bw_indices)
        error('BW Specification %.2f MHz is not reached at any design point.', target_bw / 1e6);
    end
    
    % Find index with minimum Zeta error from 0.707
    [~, closest_zeta_idx] = min(abs(mat_Zeta(valid_bw_indices) - target_zeta));
    valid_indices = valid_bw_indices(closest_zeta_idx);
end

% Find minimum ID among points satisfying bandwidth criteria
[~, min_loc_idx]              = min(mat_ID(valid_indices));
best_linear_idx               = valid_indices(min_loc_idx);
[best_gmid_idx, best_rd_idx]  = ind2sub(size(mat_BW), best_linear_idx);

% Selected Optimal Parameters
opt_ID   = mat_ID(best_gmid_idx, best_rd_idx);
opt_BW   = mat_BW(best_gmid_idx, best_rd_idx);
opt_GBW  = A_v0 * opt_BW;
opt_gmID = gm_ID_sweep(best_gmid_idx);
opt_Rd   = R_d_sweep(best_rd_idx);

% Optimal Transistor Sizing & Parameter Extraction for Cadence Virtuoso
opt_W  = opt_ID / look_up(nch, 'ID_W', 'GM_ID', opt_gmID, 'L', L_channel);
opt_gm = opt_ID * opt_gmID;

opt_Cgs = opt_W * look_up(nch, 'CGS_W', 'GM_ID', opt_gmID, 'L', L_channel);
opt_Cgd = opt_W * look_up(nch, 'CGD_W', 'GM_ID', opt_gmID, 'L', L_channel);
opt_Cdd = opt_W * look_up(nch, 'CDD_W', 'GM_ID', opt_gmID, 'L', L_channel);
opt_Cdb = opt_Cdd - abs(opt_Cgd);

% Equivalent RLC Components for Butterworth Response Verification
C_in_opt = C_in_ext + abs(opt_Cgs);
C_f_opt  = abs(opt_Cgd);
C_L_opt  = opt_Cdb + C_L_ext;

% Damping Factor (Zeta) & Natural Frequency (Wn) Calculation
num_zeta = R_f * (1 + opt_gm * opt_Rd) * C_f_opt + opt_Rd * C_L_opt + (R_f + opt_Rd) * C_in_opt;
den_zeta = 2 * sqrt(R_f * opt_Rd * (1 + opt_gm * opt_Rd) * (C_f_opt * C_L_opt + C_L_opt * C_in_opt + C_in_opt * C_f_opt));
zeta_opt = num_zeta / den_zeta;

omega_n          = sqrt(opt_gm / (R_f * (C_f_opt * C_L_opt + C_L_opt * C_in_opt + C_in_opt * C_f_opt)));
f_3dB_theoretical = omega_n / (2 * pi);

C_I      = abs(opt_Cgs) + abs(opt_Cgd);
C_T      = C_in_ext + C_I;
CI_ratio = C_I / C_in_ext;

% f_3dB_theoretical2 = 1 / (sqrt(2) * A_v0 * R_f * C_T);

%% =========================================================================
%  INTEGRATED NOISE ANALYSIS
%  =========================================================================
kB      = 1.3806488e-23;
Tn      = nch.TEMP;
sth_tgt_n = look_up(nch, 'STH_GM', 'GM_ID', opt_gmID, 'L', L_channel);
gamma_n   = sth_tgt_n / (4 * kB * Tn);

Gamma   = gamma_n;

K_noise  = 17.6 * pi * kB * Tn * (f_3dB_theoretical^3);
term_gbw = C_T / opt_GBW;
term_f2  = Gamma * (C_T^2) / (fT_selected * C_I);

I_n_in_sq  = 2 * K_noise * (term_gbw + term_f2);
I_n_in_rms = sqrt(I_n_in_sq);

% I_n_out_sq         = 8 * kB * Tn * (gamma_n * opt_gm + (1 / opt_Rd));
% V_n_in_sq          = I_n_out_sq / (opt_gm^2);
% integrated_noise_sq = V_n_in_sq * ((pi / 2) * opt_BW);
% V_n_in_rms         = sqrt(integrated_noise_sq);

%% =========================================================================
%  DISPLAY DESIGN RESULTS TO COMMAND WINDOW
%  =========================================================================
fprintf('\n================== TIA DESIGN OPTIMIZATION RESULTS ==================\n');
fprintf('Target BW Specification   : %.2f MHz\n', target_bw / 1e6);
fprintf('Calculated Bandwidth (BW) : %.2f MHz\n', opt_BW / 1e6);
fprintf('Minimum Bias Current (ID) : %.2f uA\n', opt_ID / 1e6);
% fprintf('Integrated Noise Sq       : %.3e V^2\n', integrated_noise_sq);
% fprintf('Input RMS Noise Voltage   : %.3e V rms\n', V_n_in_rms);
fprintf('Damping Factor (Zeta)     : %.3f (Target Butterworth = 0.707)\n', zeta_opt);
fprintf('---------------------------------------------------------------------\n');
fprintf('DESIGN PARAMETERS FOR CADENCE VIRTUOSO:\n');
fprintf('Transistor Length (L)     : %.2f nm\n', L_channel / 1e-9);
fprintf('Transistor Width (W)      : %.2f um\n', opt_W / 1e-6);
fprintf('Selected gm/ID Value      : %.2f V^-1\n', opt_gmID);
fprintf('Transconductance (gm)     : %.2f mS\n', opt_gm * 1e3);
fprintf('Drain Resistor (RD)       : %.2f Ohm\n', opt_Rd);
fprintf('Feedback Resistor (RF)    : %.2f Ohm\n', R_f);
fprintf('=====================================================================\n');

% Highlight Optimal Design Point on Plot
hold on;
plot(opt_ID / 1e-6, opt_BW / 1e6, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'r');
legend('Sweep Curve', 'Optimal Design Point', 'Location', 'best');
hold off;
