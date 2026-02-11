%==========================================================================
% mainMAG.m - Magnetometer (MAG) Simulation & Calibration Framework
%==========================================================================
%
% DESCRIPTION:
%   Simulation framework for generating high-fidelity synthetic 3-axis
%   Magnetometer data for spacecraft AOCS/GNC. This script models
%   environmental magnetic fields (IGRF) and realistic sensor error sources
%   including hard iron (bias), soft iron (scale/non-orthogonality), and
%   stochastic noise (wideband + bias instability).
%
% PURPOSE:
%   1. Define hardware specifications for typical space-grade Fluxgates/AMR.
%   2. Generate "Truth" magnetic field vectors (Orbit + IGRF model).
%   3. Simulate sensor measurement corruption.
%   4. Validate magnetometer calibration algorithms (e.g., TWO-STEP).
%
% ALGORITHM:
%   The measurement model follows the standard equation:
%     B_meas = T_det * B_true + b_hard_iron + b_temp + noise
%   
%   where:
%     T_det       = M_SoftIron * M_NonOrth * DCM_Mounting
%     M_SoftIron  = Scale Factor / Soft Iron matrix (near identity)
%     M_NonOrth   = Non-orthogonality matrix
%     b_hard_iron = Constant offset (Hard Iron + Electronic Bias)
%     b_temp      = Temperature-dependent slowly varying bias
%     noise       = Wideband Gaussian noise
%
%==========================================================================

clear;
close all;
clc;

%% ========================================================================
% 1. MAG HARDWARE PARAMETERS (typical values from Space-Grade Fluxgate/AMR datasheets)
% =========================================================================

% --- General Configuration ---
MAG.rate = 10;                         % Sampling frequency [Hz]
MAG.dt   = 1 / MAG.rate;       

% Dynamic Range & Resolution       
MAG.range      = 100000;               % [nT] +/- full scale
MAG.resolution = 0.5;                  % [nT] ADC1ae4ffw quantization step

% Mounting Misalignment: Error aligning the IMU case to the Body frame.
MAG.mountingError = [0.2; 0.2; 0.5];   % [deg] roll, pitch, yaw

% --- Stochastic Noise Parameters ---
MAG.noiseDensity = 0.1;                % [nT/sqrt(Hz)] white noise spectral density

% Anti-alias / measurement bandwidth model for noise integration
MAG.antiAlias.type = "onePole";        % "idealNyquist" | "onePole"
MAG.antiAlias.fc   = 0.25 * MAG.rate;  % [Hz] 1-pole cutoff used for ENBW (if onePole)

% Bias drift as bounded 1st-order Gauss-Markov process
MAG.biasInstability = 5.0;             % [nT] steady-state 1-sigma (per axis)
MAG.biasTau         = 3600;            % [s] correlation time (time constant)

% --- Deterministic Errors ---
MAG.hardIronLimit = 500.0;             % [nT] static hard-iron magnitude limit (per axis, uniform)
MAG.SF = 1000;                         % [ppm] diagonal scale errors (1-sigma)
MAG.softIronCoupling = 500;            % [ppm] off-diagonal symmetric coupling (1-sigma)
MAG.nonOrthogonality  = 0.1;           % [deg] small non-orth angles (1-sigma)


% --- Derived / Initialization ---
fprintf('\n=== Initializing Magnetometer Model ===\n');

% 1) Mounting misalignment (Body -> Sensor) true DCM
roll_err  = deg2rad(MAG.mountingError(1));
pitch_err = deg2rad(MAG.mountingError(2));
yaw_err   = deg2rad(MAG.mountingError(3));
MAG.DCM_B2S_true = angle2dcm(yaw_err, pitch_err, roll_err, 'ZYX');

% 2) Soft iron / scale factor (symmetric matrix near identity)
SF_err = (MAG.SF * 1e-6) .* randn(3,1);
C_ppm  = (MAG.softIronCoupling * 1e-6) .* randn(3,1); % xy, xz, yz couplings
MAG.M_SoftIron = eye(3) + [ SF_err(1),  C_ppm(1),  C_ppm(2);
                             C_ppm(1), SF_err(2),  C_ppm(3);
                             C_ppm(2),  C_ppm(3), SF_err(3)];
 
fprintf(' Soft-Iron / Scale-Factor Matrix: [%+.6f  %+.6f  %+.6f]\n', MAG.M_SoftIron(1,1), MAG.M_SoftIron(1,2), MAG.M_SoftIron(1,3));
fprintf('                                  [%+.6f  %+.6f  %+.6f]\n', MAG.M_SoftIron(2,1), MAG.M_SoftIron(2,2), MAG.M_SoftIron(2,3));
fprintf('                                  [%+.6f  %+.6f  %+.6f]\n', MAG.M_SoftIron(3,1), MAG.M_SoftIron(3,2), MAG.M_SoftIron(3,3));

% 3) Non-orthogonality (small-angle upper triangular correction)
nonorth_err = deg2rad(MAG.nonOrthogonality) .* randn(3,1);
MAG.M_NonOrth = [              1, -nonorth_err(3),  nonorth_err(2);
                 -nonorth_err(3),               1, -nonorth_err(1);
                  nonorth_err(2), -nonorth_err(1),               1];

% Deterministic sensor matrix acting in Sensor axes (after mounting)
MAG.M_Deterministic = MAG.M_SoftIron * MAG.M_NonOrth;

% 4) Static hard-iron bias (Sensor axes)
MAG.biasHardIron = MAG.hardIronLimit * (2*rand(3,1) - 1);

fprintf(' Hard Iron Bias (S):              [%.2f, %.2f, %.2f] nT\n', MAG.biasHardIron(1), MAG.biasHardIron(2), MAG.biasHardIron(3));

% 5) White noise RMS per sample from density and assumed bandwidth
if string(MAG.antiAlias.type) == "idealNyquist"
    BW = MAG.antiAlias.fc/2;
elseif string(MAG.antiAlias.type) == "onePole"
    BW = 1.57 * MAG.antiAlias.fc; % ENBW for 1-pole low-pass.
else
    BW = MAG.antiAlias.fc/2;
end
MAG.sigmaWhite = MAG.noiseDensity * sqrt(BW);

fprintf(' White Noise RMS:                 %.4f nT (per axis)\n', MAG.sigmaWhite);

fprintf('=== MAG Initialization Complete ===\n');

%% ========================================================================
% 2. SIMULATION SETTINGS & GROUND TRUTH GENERATION
% =========================================================================

% --- Simulation time settings ---
startEpoch = datetime(2026,2,10,12,0,0,'TimeZone','UTC');
T_total = 60*20;      % Total simulation time [s]
t = 0:MAG.dt:T_total; % Time vector [s]
N = numel(t);

fprintf('\n=== Generating Truth Trajectory ===\n');
fprintf(' Simulation time: %.1f s (%d samples)\n', T_total, N);

% --- Simple circular orbit model (ECI position) ---

fprintf(' Propagating orbit (J2 + Drag + SRP) ... \n');

% Spacecraft Parameters
scParams.Mass      = 100.0; % [kg]
scParams.Area_drag = 1.0;   % [m^2]
scParams.Area_srp  = 1.2;   % [m^2] (Usually larger due to panels)
scParams.Cd        = 2.2;   % [-] Drag Coeff
scParams.Cr        = 1.5;   % [-] Radiation Pressure Coeff (1=Absorb, 2=Reflect)

% Orbital Elements
orbElements.sma      = 6378e3 + 500e3; 
orbElements.ecc      = 0.001;
orbElements.inc_deg  = 97.0;
orbElements.raan_deg = 40.0;
orbElements.argp_deg = 0.0;
orbElements.nu_deg   = 0.0;
orbElements.epoch    = startEpoch;

[r_ECI, v_ECI_unused] = simulateOrbit(t, orbElements, scParams);

% --- ECI to ECEF and Geodetic Coordinates (LLA) ---

% Calculate DCMs for all time steps at once (dcmeci2ecef handles vector time)
DCM_ECI2ECEF_all = dcmeci2ecef('IAU-2000/2006', startEpoch + seconds(t'));

% Rotate ECI positions to ECEF
r_ECEF = zeros(3,N);
for k = 1:N
    % Transpose the 3x3 page for the k-th sample if necessary, 
    % but dcmeci2ecef output format is 3x3xN.
    r_ECEF(:,k) = DCM_ECI2ECEF_all(:,:,k) * r_ECI(:,k);
end

% Convert ECEF to Geodetic (LLA) using WGS84
% ecef2lla expects Nx3 inputs [x, y, z] and returns [lat, lon, alt]
lla = ecef2lla(r_ECEF', 'WGS84'); 

lat = lla(:,1)'; % [deg] 1xN
lon = lla(:,2)'; % [deg] 1xN
alt = lla(:,3)'; % [m]   1xN

% --- Magnetic Field Truth ---

fprintf(' Using IGRF model ... \n');

B_true_ECI   = zeros(3,N);
decimal_year = 2025 * ones(1,N);

XYZ_NED                        = igrfmagm(alt, lat, lon, decimal_year, 13);
[B_ECEF_x, B_ECEF_y, B_ECEF_z] = ned2ecefv(XYZ_NED(:,1), XYZ_NED(:,2), XYZ_NED(:,3), ...
                                           lat(:), lon(:));
B_ECEF_all                     = [B_ECEF_x'; B_ECEF_y'; B_ECEF_z'];

DCM_ECEF2ECI_all = permute(DCM_ECI2ECEF_all, [2, 1, 3]);
for k = 1:N
    B_true_ECI(:,k) = DCM_ECEF2ECI_all(:,:,k) * B_ECEF_all(:,k);
end

% --- Truth Attitude & Body Frame Field ---
angle_true = deg2rad(30);                 % Rotation angle [rad]
axis_true = [1; 1; 1] / sqrt(3);          % Unit rotation axis

q_true      = zeros(4,N);
q_true(:,1) = [cos(angle_true/2); ...          % Scalar part
               sin(angle_true/2) * axis_true]; % Vector part

omega_true_b = zeros(3,N); % [rad/s] body rates
omega_true_b(1,:) = deg2rad(0.15);
omega_true_b(2,:) = deg2rad(0.10) * sin(2*pi*(1/300)*t);
omega_true_b(3,:) = deg2rad(0.08) * cos(2*pi*(1/500)*t);

for k = 2:N
    q_dot = 1/2 * skewOmega_custom(omega_true_b(:,k-1)) * q_true(:,k-1);
    q_true(:,k) = q_true(:,k-1) + q_dot * MAG.dt;
    q_true(:,k) = q_true(:,k) / norm(q_true(:,k));
end

B_true = zeros(3,N);
for k = 1:N
    DCM_ECI2B = quat2dcm_custom(q_true(:,k));
    B_true(:,k) = DCM_ECI2B * B_true_ECI(:,k);
end

fprintf('\n=== MAG Ground Truth Complete (Mean ||B||: %.0f nT) ===\n', mean(vecnorm(B_true)));

%% ========================================================================
% 3. GENERATE SYNTHETIC MAGNETOMETER MEASUREMENTS
% =========================================================================

fprintf('\n=== Generating Synthetic Magnetometer Measurements ===\n');

mag_meas = generateMAGMeasurements(t, B_true, MAG);

%% ========================================================================
% 4. VALIDATION AND ERROR ANALYSIS (CONSISTENT FRAMES)
% =========================================================================
fprintf('\n=== MAG Validation and Error Analysis ===\n');

% Truth in Sensor frame (before deterministic sensor matrix)
B_true_S = MAG.DCM_B2S_true * B_true;

% Error against "ideal sensor axes truth" (shows total corruption)
B_err = mag_meas.B_clean - B_true_S;
rms_B_nT = sqrt(mean(B_err.^2, 2));

fprintf('\n--- RMS Error vs. Ideal Sensor Truth ---\n');
fprintf('  RMS(B_sx):  %.2f nT\n', rms_B_nT(1));
fprintf('  RMS(B_sy):  %.2f nT\n', rms_B_nT(2));
fprintf('  RMS(B_sz):  %.2f nT\n', rms_B_nT(3));

% Separate white noise estimate by subtracting known deterministic and bias terms
B_det_S = mag_meas.B_det;
b_total = mag_meas.bias_total;
w_hat   = mag_meas.B_clean - (B_det_S + b_total); % should be ~white (pre-quantization)
white_std = std(w_hat, 0, 2);

fprintf('\n--- Empirical White Noise (Pre-Quant) ---\n');
fprintf('  Std: [%.3f, %.3f, %.3f] nT\n', white_std(1), white_std(2), white_std(3));
fprintf('  Spec RMS (theory): %.3f nT\n', MAG.sigmaWhite);

% Direction error (unit vectors) vs ideal sensor truth
u_true = B_true_S ./ vecnorm(B_true_S, 2, 1);
u_meas = mag_meas.B_meas ./ vecnorm(mag_meas.B_meas, 2, 1);
cang = sum(u_true .* u_meas, 1);
cang = max(min(cang, 1), -1);
angular_err_deg = rad2deg(acos(cang));

fprintf('\n--- Angular Alignment Error (Measured vs Ideal Sensor Truth) ---\n');
fprintf('  Mean:   %.4f deg\n', mean(angular_err_deg));
fprintf('  Std:    %.4f deg\n', std(angular_err_deg));
fprintf('  Max:    %.4f deg\n', max(angular_err_deg));

%% ========================================================================
% 5. RESULTS VISUALIZATION
% =========================================================================

% Generate comprehensive validation plots
saveFlag = 1;
plotMAGResults(t, B_true, mag_meas, MAG, saveFlag);

fprintf('\n=== All tasks completed successfully ===\n');