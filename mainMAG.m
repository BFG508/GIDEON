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
% 1. MAG HARDWARE PARAMETERS (Typical Space-Grade Fluxgate)
% =========================================================================
% Representative of sensors like Billingsley TFM100S, ZARM FGM, or NewSpace.

% --- General Configuration ---
MAG.rate = 10;       % Sampling frequency [Hz] (Typical AOCS: 10-50 Hz)
MAG.dt   = 1/MAG.rate;

% Dynamic Range & Resolution
% Earth field varies approx 25,000 to 65,000 nT (LEO).
MAG.range      = 100000; % [nT] Dynamic range (+/-)
MAG.resolution = 0.5;    % [nT] Digital resolution (e.g. 18-bit ADC over range)

% Mounting Misalignment (Mechanical)
% Error aligning the magnetometer bracket to the spacecraft Body frame.
MAG.mountingError = [0.2; 0.2; 0.5]; % [deg] Roll, Pitch, Yaw

% --- 1. Stochastic Noise Parameters ---
% Wideband Noise (White Noise Density)
% Typical Fluxgate: < 20 pT/sqrt(Hz). AMR/Commercial: ~1-10 nT/sqrt(Hz).
% We choose a conservative AOCS value:
MAG.noiseDensity = 0.1;  % [nT/sqrt(Hz)] Spectral density
MAG.sigma_white  = MAG.noiseDensity * sqrt(MAG.rate/2); % [nT] RMS (approx BW = rate/2)

% Bias Stability (Random Walk / Flicker)
% Drift over temperature and time.
MAG.biasInstability = 5.0; % [nT] 1-sigma over orbital period

% --- 2. Deterministic Errors (Hard Iron / Bias) ---
% "Hard Iron" is the combination of sensor zero-offset and the spacecraft's
% own static magnetic remanence. This is usually the largest error source.
MAG.hardIronLimit = 500.0; % [nT] Max expected offset (uncalibrated)

% --- 3. Scale Factor & Soft Iron Errors ---
% "Soft Iron" distorts the field shape (ellipsoid effect).
% Modeled as Scale Factor error + Off-diagonal coupling.
MAG.ScaleFactorErr = 1000; % [ppm] (0.1% linearity error)

% --- 4. Geometric Errors (Non-Orthogonality) ---
% Intrinsic sensor axis misalignment (deviation from 90 deg).
MAG.nonOrthogonality = 0.1; % [deg] (1-sigma)

% --- Derived Matrices & Pre-Allocation ---
fprintf('\n=== Initializing Magnetometer Model ===\n');

    % 1. Construct Mounting Misalignment Matrix (DCM Body -> Case)
    roll_err  = deg2rad(MAG.mountingError(1));
    pitch_err = deg2rad(MAG.mountingError(2));
    yaw_err   = deg2rad(MAG.mountingError(3));
    MAG.DCM_mounting = angle2dcm(yaw_err, pitch_err, roll_err, 'ZYX');

    % 2. Construct Soft Iron / Scale Factor Matrix (Diagonal + Symm)
    % S = I + diag(sf_err)
    sf_err = (MAG.ScaleFactorErr * 1e-6) * randn(3,1);
    MAG.M_SoftIron = eye(3) + diag(sf_err);
    
    fprintf(' Scale Factor Errors: [%.1f, %.1f, %.1f] ppm\n', ...
            sf_err(1)*1e6, sf_err(2)*1e6, sf_err(3)*1e6);

    % 3. Construct Non-Orthogonality Matrix (Upper Triangular)
    % Corrects the sensor frame to be orthogonal
    no_err = deg2rad(MAG.nonOrthogonality) * randn(3,1);
    MAG.M_NonOrth = [1, -no_err(3),  no_err(2);
                     0,          1, -no_err(1);
                     0,          0,         1];

    % 4. Total Deterministic Transform (Body True -> Sensor Meas)
    % T_mag = M_SoftIron * M_NonOrth * DCM_Mounting
    MAG.T_Deterministic = MAG.M_SoftIron * MAG.M_NonOrth * MAG.DCM_mounting;

    % 5. Generate Static Hard Iron Bias (Constant for this run)
    MAG.biasHardIron = MAG.hardIronLimit * (2*rand(3,1) - 1); % Uniform +/- limit
    
    fprintf(' Hard Iron Bias:      [%.2f, %.2f, %.2f] nT\n', ...
            MAG.biasHardIron(1), MAG.biasHardIron(2), MAG.biasHardIron(3));
    fprintf(' Noise RMS (White):   %.3f nT\n', MAG.sigma_white);

fprintf('=== MAG Initialization Complete ===\n');

%% ========================================================================
% 2. SIMULATION SETTINGS & GROUND TRUTH GENERATION
% =========================================================================
% This section generates the magnetic field "truth" in BODY frame:
%   1) Define time vector and orbital parameters
%   2) Propagate orbit (ECI) and convert to Geodetic coordinates (LLA)
%   3) Compute Earth's magnetic field using IGRF-13/14 (if available)
%      via MATLAB's 'igrfmagm' function (Aerospace Toolbox).
%   4) Rotate B-field: NED (IGRF output) -> ECEF -> ECI -> BODY
%
% Notes:
%   - 'igrfmagm' returns B-field in local NED frame [North, East, Down].
%   - Fallback to Dipole model provided if Aerospace Toolbox is missing.
% =========================================================================

fprintf('\n=== Generating MAG Ground Truth ===\n');

% -------------------------------------------------------------------------
% 2.1 Simulation time settings
% -------------------------------------------------------------------------
T_total = 60*20;               % Total simulation time [s] (20 minutes)
t = 0:MAG.dt:T_total;          % Time vector [s]
N = numel(t);

fprintf(' Simulation time: %.1f s (%d samples @ %.1f Hz)\n', T_total, N, MAG.rate);

% -------------------------------------------------------------------------
% 2.2 Simple circular orbit model (ECI position)
% -------------------------------------------------------------------------
mu_E = 3.986004418e14;         % Earth's gravitational parameter [m^3/s^2]
R_E  = 6378137.0;              % Earth mean equatorial radius [m]
w_E  = 7.2921150e-5;           % Earth rotation rate [rad/s]

orb.altitude_m      = 500e3;   % [m] Typical LEO altitude
orb.inclination_deg = 97.0;    % [deg] Sun-Synchronous-like
orb.RAAN_deg        = 40.0;    % [deg]
orb.argp_deg        = 0.0;     % [deg]
orb.nu0_deg         = 0.0;     % [deg] Initial true anomaly

a = R_E + orb.altitude_m;      % Semi-major axis [m]
n = sqrt(mu_E / a^3);          % Mean motion [rad/s]

% Perifocal position: r_pf = a*[cos(nu); sin(nu); 0]
nu = deg2rad(orb.nu0_deg) + n*t;
r_pf = [a*cos(nu); a*sin(nu); zeros(1,N)];

% Perifocal -> ECI rotation
RAAN = deg2rad(orb.RAAN_deg);
inc  = deg2rad(orb.inclination_deg);
argp = deg2rad(orb.argp_deg);

R3_RAAN = [ cos(RAAN) sin(RAAN) 0; -sin(RAAN) cos(RAAN) 0; 0 0 1];
R1_inc  = [ 1 0 0; 0 cos(inc) sin(inc); 0 -sin(inc) cos(inc)];
R3_argp = [ cos(argp) sin(argp) 0; -sin(argp) cos(argp) 0; 0 0 1];

DCM_P2ECI = (R3_RAAN' * R1_inc' * R3_argp')'; 
r_ECI = DCM_P2ECI * r_pf;      % [m], 3xN

% -------------------------------------------------------------------------
% 2.3 ECI -> ECEF and Geodetic Coordinates (LLA)
% -------------------------------------------------------------------------
theta = w_E * t; % Greenwich Hour Angle (simplified)

r_ECEF  = zeros(3,N);
lat_deg = zeros(1,N);
lon_deg = zeros(1,N);
alt_m   = zeros(1,N);

for k = 1:N
    cth = cos(theta(k)); sth = sin(theta(k));
    DCM_ECI2ECEF = [ cth sth 0; -sth cth 0; 0 0 1];
    
    r_ECEF(:,k) = DCM_ECI2ECEF * r_ECI(:,k);
    
    % Spherical approximation for LLA (sufficient for simulation inputs)
    % For high precision, use 'ecef2lla' if available.
    rk  = norm(r_ECEF(:,k));
    lat = asin(r_ECEF(3,k) / rk);
    lon = atan2(r_ECEF(2,k), r_ECEF(1,k));
    
    lat_deg(k) = rad2deg(lat);
    lon_deg(k) = rad2deg(lon);
    alt_m(k)   = rk - R_E;
end

% -------------------------------------------------------------------------
% 2.4 Magnetic Field Truth (IGRF-13 via 'igrfmagm')
% -------------------------------------------------------------------------
B_true_ECI_nT = zeros(3,N);
decimal_year = 2025.0;

% Check for Aerospace Toolbox function 'igrfmagm'
if exist('igrfmagm', 'file') == 2
    fprintf(' using IGRF model (igrfmagm) ... ');
    
    % Expand decimal_year to match input vector sizes
    decimal_year_vec = decimal_year * ones(1, N);
    
    % Call igrfmagm (vectorized)
    % Output XYZ is in North-East-Down (NED) frame in nT
    [XYZ_NED, ~, ~, ~, ~] = igrfmagm(alt_m, lat_deg, lon_deg, decimal_year_vec);
    
    for k = 1:N
        % Extract NED vector for this sample
        B_NED = XYZ_NED(k, :)'; % [North; East; Down]
        
        % Convert NED -> ECEF
        % Transformation depends on Latitude (phi) and Longitude (lambda)
        phi = deg2rad(lat_deg(k));
        lam = deg2rad(lon_deg(k));
        
        cphi = cos(phi); sphi = sin(phi);
        clam = cos(lam); slam = sin(lam);
        
        % Rotation matrix from NED to ECEF (local tangent plane)
        M_NED2ECEF = [ -sphi*clam, -slam, -cphi*clam;
                       -sphi*slam,  clam, -cphi*slam;
                        cphi,       0,    -sphi     ];
                   
        B_ECEF = M_NED2ECEF * B_NED;
        
        % Convert ECEF -> ECI
        cth = cos(theta(k)); sth = sin(theta(k));
        DCM_ECEF2ECI = [ cth -sth 0; sth cth 0; 0 0 1];
        
        B_true_ECI_nT(:,k) = DCM_ECEF2ECI * B_ECEF;
    end
    fprintf('Done.\n');
    
else
    fprintf(' warning: "igrfmagm" not found. Using Dipole Model fallback.\n');
    B_equator_nT = 30000;
    K_dipole = B_equator_nT * R_E^3;
    m_hat = [0;0;1];
    
    for k = 1:N
        r_vec = r_ECEF(:,k); 
        r_n = norm(r_vec); 
        u_r = r_vec/r_n;
        
        B_ECEF = (K_dipole/r_n^3) * (3*dot(m_hat, u_r)*u_r - m_hat);
        
        cth = cos(theta(k)); sth = sin(theta(k));
        DCM_ECEF2ECI = [ cth -sth 0; sth cth 0; 0 0 1];
        B_true_ECI_nT(:,k) = DCM_ECEF2ECI * B_ECEF;
    end
end

% -------------------------------------------------------------------------
% 2.5 Truth Attitude & Body Frame Field
% -------------------------------------------------------------------------
% Define attitude profile (slow tumble for sensor excitation)
angle0 = deg2rad(20);
axis0  = [1; 1; 1] / sqrt(3);
q_true = zeros(4,N); 
q_true(:,1) = [cos(angle0/2); sin(angle0/2)*axis0];

omega_true_b = zeros(3,N);
omega_true_b(1,:) = deg2rad(0.15);
omega_true_b(2,:) = deg2rad(0.10) * sin(2*pi*(1/300)*t);
omega_true_b(3,:) = deg2rad(0.08) * cos(2*pi*(1/500)*t);

% Propagate attitude
for k = 2:N
    wx = omega_true_b(1,k-1); wy = omega_true_b(2,k-1); wz = omega_true_b(3,k-1);
    Omega = [0 -wx -wy -wz; wx 0 wz -wy; wy -wz 0 wx; wz wy -wx 0];
    q_true(:,k) = q_true(:,k-1) + 0.5 * Omega * q_true(:,k-1) * MAG.dt;
    q_true(:,k) = q_true(:,k) / norm(q_true(:,k));
end

% Rotate Field to Body Frame
B_true_B_nT = zeros(3,N);
for k = 1:N
    DCM_ECI2B = quat2dcm_custom(q_true(:,k)); 
    B_true_B_nT(:,k) = DCM_ECI2B * B_true_ECI_nT(:,k);
end

fprintf('=== MAG Ground Truth Complete (Mean Field: %.0f nT) ===\n', mean(vecnorm(B_true_B_nT)));

%% ========================================================================
% 3. GENERATE SYNTHETIC MAGNETOMETER MEASUREMENTS
% =========================================================================
fprintf('\n=== Generating Synthetic Magnetometer Measurements ===\n');

mag_meas = generateMAGMeasurements(t, B_true_B_nT, MAG);

fprintf('=== Measurement Generation Complete ===\n');

%% ========================================================================
% 4. VALIDATION AND ERROR ANALYSIS
% =========================================================================
fprintf('\n=== MAG Validation and Error Analysis ===\n');

% Measurement error (excluding quantization for analysis)
B_err_nT = mag_meas.B_clean_nT - B_true_B_nT;  % [nT]

% RMS error per axis (total corruption before quantization)
rms_B_nT = sqrt(mean(B_err_nT.^2, 2));         % [nT]

fprintf('\n--- Magnetic Field RMS Error ---\n');
fprintf('  RMS(B_x):  %.2f nT\n', rms_B_nT(1));
fprintf('  RMS(B_y):  %.2f nT\n', rms_B_nT(2));
fprintf('  RMS(B_z):  %.2f nT\n', rms_B_nT(3));
fprintf('  RMS(||B||): %.2f nT\n', sqrt(sum(rms_B_nT.^2)));

% Field magnitude error
B_true_mag = vecnorm(B_true_B_nT, 2, 1);
B_meas_mag = vecnorm(mag_meas.B_meas_nT, 2, 1);
mag_error_nT = B_meas_mag - B_true_mag;

fprintf('\n--- Magnetic Field Magnitude Error ---\n');
fprintf('  Mean error:  %.2f nT\n', mean(mag_error_nT));
fprintf('  Std dev:     %.2f nT\n', std(mag_error_nT));
fprintf('  Max error:   %.2f nT\n', max(abs(mag_error_nT)));

% Angular error (unit vector misalignment)
% Angle between true and measured field directions
u_true = B_true_B_nT ./ vecnorm(B_true_B_nT, 2, 1);
u_meas = mag_meas.B_meas_nT ./ vecnorm(mag_meas.B_meas_nT, 2, 1);

angular_err_rad = acos(max(min(dot(u_true, u_meas), 1), -1));
angular_err_deg = rad2deg(angular_err_rad);

fprintf('\n--- Angular Alignment Error (Direction) ---\n');
fprintf('  Mean:   %.4f deg\n', mean(angular_err_deg));
fprintf('  Std:    %.4f deg\n', std(angular_err_deg));
fprintf('  Max:    %.4f deg  ', max(angular_err_deg));

% Performance assessment for AOCS attitude determination
% Typical requirement: < 0.5 deg for fine pointing
if max(angular_err_deg) < 0.1
    fprintf('✓ Excellent (better than 0.1°)\n');
elseif max(angular_err_deg) < 0.5
    fprintf('✓ Good (within 0.5° AOCS spec)\n');
elseif max(angular_err_deg) < 2.0
    fprintf('○ Acceptable (coarse AOCS)\n');
else
    fprintf('⚠ High error (check calibration)\n');
end

% Dynamic bias drift statistics
bias_dyn_nT = mag_meas.bias_hard_dyn;  % [nT]
mean_bias_dyn = mean(bias_dyn_nT, 2);
std_bias_dyn  = std(bias_dyn_nT, 0, 2);
max_bias_dyn  = max(abs(bias_dyn_nT), [], 2);

fprintf('\n--- Dynamic Bias Drift Statistics ---\n');
fprintf('  Axis    Mean [nT]    Std [nT]     Max [nT]\n');
fprintf('  ───────────────────────────────────────────\n');
fprintf('  X       %+8.2f    %7.2f    %7.2f\n', mean_bias_dyn(1), std_bias_dyn(1), max_bias_dyn(1));
fprintf('  Y       %+8.2f    %7.2f    %7.2f\n', mean_bias_dyn(2), std_bias_dyn(2), max_bias_dyn(2));
fprintf('  Z       %+8.2f    %7.2f    %7.2f\n', mean_bias_dyn(3), std_bias_dyn(3), max_bias_dyn(3));

% Check consistency with specification
fprintf('\n--- Performance vs. Specifications ---\n');
fprintf('  Noise Density Spec:     %.2f nT/√Hz\n', MAG.noiseDensity);
fprintf('  Bias Instability Spec:  %.1f nT\n', MAG.biasInstability);

% Empirical noise estimation (high-frequency component after detrending)
B_err_detrended = B_err_nT - mean(B_err_nT, 2); % Remove mean
empirical_noise_std = std(B_err_detrended, 0, 2);

fprintf('  Empirical Noise (RMS):  [%.2f, %.2f, %.2f] nT\n', ...
        empirical_noise_std(1), empirical_noise_std(2), empirical_noise_std(3));

% Theoretical white noise RMS from spec
theoretical_noise_rms = MAG.sigma_white;
fprintf('  Theoretical Noise:      %.2f nT  ', theoretical_noise_rms);

% Validate consistency
noise_ratio = mean(empirical_noise_std) / theoretical_noise_rms;
if abs(noise_ratio - 1.0) < 0.2
    fprintf('✓ Consistent\n');
elseif abs(noise_ratio - 1.0) < 0.5
    fprintf('○ Acceptable\n');
else
    fprintf('⚠ Deviation detected (ratio = %.2f)\n', noise_ratio);
end

% Total error budget breakdown
total_rms = sqrt(sum(rms_B_nT.^2));
hard_iron_magnitude = norm(MAG.biasHardIron);
soft_iron_magnitude = norm(MAG.T_Deterministic - eye(3), 'fro') * mean(B_true_mag);
noise_magnitude = norm(empirical_noise_std);

fprintf('\n--- Error Budget Breakdown (Approximate) ---\n');
fprintf('  Hard Iron Bias:     %.1f nT  (%.1f%%)\n', ...
        hard_iron_magnitude, 100*hard_iron_magnitude/total_rms);
fprintf('  Soft Iron/Scale:    %.1f nT  (%.1f%%)\n', ...
        soft_iron_magnitude, 100*soft_iron_magnitude/total_rms);
fprintf('  White Noise:        %.1f nT  (%.1f%%)\n', ...
        noise_magnitude, 100*noise_magnitude/total_rms);
fprintf('  Total RMS:          %.1f nT\n', total_rms);

% Calibration quality indicator
if hard_iron_magnitude < 50
    fprintf('\n  Calibration status: ✓ Well-calibrated (Hard Iron < 50 nT)\n');
elseif hard_iron_magnitude < 200
    fprintf('\n  Calibration status: ○ Nominal (Hard Iron < 200 nT)\n');
else
    fprintf('\n  Calibration status: ⚠ Requires calibration (Hard Iron > 200 nT)\n');
end

fprintf('\n=== MAG Error Analysis Completed ===\n');

%% ========================================================================
% 5. RESULTS VISUALIZATION
% =========================================================================
% Generate comprehensive validation plots

saveFlag = 1;
%plotMAGResults(t, B_true_B_nT, mag_meas, MAG, mag_truth, saveFlag);

fprintf('\n=== All tasks completed successfully ===\n');
