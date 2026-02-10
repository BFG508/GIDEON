%==========================================================================
% mainIMU.m - Inertial Measurement Unit (IMU) Simulation & Calibration
%==========================================================================

% DESCRIPTION:
% Simulation framework for generating high-fidelity synthetic IMU data
% (Gyroscopes and Accelerometers) for AOCS/GNC development. This script
% models realistic error sources including deterministic errors (bias, 
% scale factor, misalignment) and stochastic noise processes (Angle/Velocity
% Random Walk, Bias Instability).

% PURPOSE:
% 1. Define hardware specifications for typical space-grade IMUs.
% 2. Generate "Truth" trajectory (angular velocity & acceleration).
% 3. Simulate sensor measurements corruption.
% 4. Validate IMU error models before integration into MEKF.

% ALGORITHM:
% The measurement model follows the standard equation:
%   u_meas = (I + M + S) * (u_true + b_static) + b_dynamic + noise
% where:
%   I = Identity matrix
%   M = Misalignment matrix (non-orthogonality + mounting errors)
%   S = Scale Factor error matrix
%   b_static = Turn-on bias (constant per run)
%   b_dynamic = Bias instability (Random Walk process)
%   noise = High-frequency white noise (ARW/VRW)

% REFERENCES:
% [1] "IEEE Standard for Inertial Sensor Terminology", IEEE Std 528-2001.
% [2] Woodman, O.J., "An Introduction to Inertial Navigation", 
%     University of Cambridge, Technical Report 696, 2007.
% [3] Grewal, M.S., et al., "Global Positioning Systems, Inertial 
%     Navigation, and Integration", Wiley, 2013.

%==========================================================================

clear;
close all;
clc;

%% ========================================================================
% 1. IMU HARDWARE PARAMETERS (Typical Space-Grade MEMS/FOG values)
% =========================================================================

% --- General Configuration ---
IMU.rate = 120;      % Sampling frequency [Hz] (Typical: 10-1000 Hz)
IMU.dt = 1/IMU.rate; % Sampling time step [s]

% Mounting Misalignment: Error aligning the IMU case to the Satellite Body frame.
IMU.mountingError = [0.1; 0.1; 0.1]; % [deg] Roll, Pitch, Yaw mounting error

% --- Gyroscope Parameters (Angular Rate) ---
% Representative of a high-end MEMS or low-end FOG (Fiber Optic Gyro)

    % 1. Stochastic Noise (Random processes)
    % ARW: High-frequency noise (white noise density). Determines short-term precision.
    IMU.gyro.ARW = 0.05;              % Angle Random Walk [deg/sqrt(h)]
    
    % Bias Instability: Long-term drift (flicker noise 1/f floor).
    IMU.gyro.biasInstability = 1.0;   % [deg/hr]
    
    % RRW: Rate Random Walk. Models the "wandering" of the bias over time.
    % Often not specified in datasheets, approximated from B.I. or set empirically.
    IMU.gyro.RRW = 0.1;               % [deg/h/sqrt(h)] (bias drift rate)
    
    % 2. Deterministic Errors (Calibration residuals)
    % Static Bias: Constant offset turned on/off (Turn-on to turn-on repeatability).
    IMU.gyro.biasStaticLimit = 10.0;  % [deg/h] (1-sigma)
    
    % Scale Factor Error: Sensitivity error (linear).
    IMU.gyro.SF = 300;                % [ppm] (parts per million)
    
    % 3. Geometric Errors
    % Internal Misalignment: Non-orthogonality of the sensor triad axes.
    IMU.gyro.nonorthogonality = 0.05; % [deg] (1-sigma)

% --- Accelerometer Parameters (Specific Force) ---
% Required for g-sensitivity corrections or Inertial Navigation (not just Attitude)

    % 1. Stochastic Noise
    IMU.accel.VRW = 50;                % Velocity Random Walk [micro-g/sqrt(Hz)]
    IMU.accel.biasInstability = 20;    % [micro-g]
      
    % 2. Deterministic Errors  
    IMU.accel.biasStaticLimit = 1.0;   % [mg]
    IMU.accel.SF = 300;                % [ppm]
    
    % 3. Geometric Errors
    IMU.accel.nonorthogonality = 0.05; % [deg]

%% ========================================================================
% 2. DERIVED MATRICES & PRE-ALLOCATION
% =========================================================================

fprintf('\n=== Initializing IMU Model ===\n');
fprintf(' Update Rate: %d Hz (dt = %.4f s)\n', IMU.rate, IMU.dt);

% Construct Mounting Misalignment Matrix (DCM Body -> Case)
% Represents the physical rotation due to imperfect mounting.
roll_err  = deg2rad(IMU.mountingError(1));
pitch_err = deg2rad(IMU.mountingError(2));
yaw_err   = deg2rad(IMU.mountingError(3));

IMU.DCM_Mounting = angle2dcm(yaw_err, pitch_err, roll_err, 'ZYX');

% Construct Scale Factor Matrix (Diagonal) - Gyro
% S = diag([sf_x, sf_y, sf_z]) where sf = scale_error_ppm * 1e-6
sf_err_gyro = (IMU.gyro.SF * 1e-6) * randn(3,1); 
IMU.gyro.M_Scale = eye(3) + diag(sf_err_gyro);

fprintf(' Gyro Scale Factor Errors: [%.1f, %.1f, %.1f] ppm\n', ...
    sf_err_gyro(1)*1e6, sf_err_gyro(2)*1e6, sf_err_gyro(3)*1e6);

% Construct Non-Orthogonality Matrix - Gyro
% Models small deviations of the axes from a perfect 90 deg triad.
% Upper triangular approximation for small angles.
nonortho_err_gyro = deg2rad(IMU.gyro.nonorthogonality) * randn(3,1);
IMU.gyro.M_NonOrth = [1, -nonortho_err_gyro(3),  nonortho_err_gyro(2);
                      0,                     1, -nonortho_err_gyro(1);
                      0,                     0,                    1];

% Total deterministic transform (Body True -> Sensor Meas)
% T_meas = M_Scale * M_NonOrth * DCM_Mounting
IMU.gyro.T_Deterministic = IMU.gyro.M_Scale * IMU.gyro.M_NonOrth * IMU.DCM_Mounting;

% Generate Static Turn-on Bias (Constant for this simulation run)
IMU.gyro.biasStatic = deg2rad(IMU.gyro.biasStaticLimit / 3600) * randn(3,1);

fprintf(' Gyro Static Turn-on Bias: [%.4f, %.4f, %.4f] deg/h\n', ...
    rad2deg(IMU.gyro.biasStatic(1))*3600, ...
    rad2deg(IMU.gyro.biasStatic(2))*3600, ...
    rad2deg(IMU.gyro.biasStatic(3))*3600);

fprintf('=== IMU Initialization Complete ===\n');

%% ========================================================================
% 3. SIMULATION SETTINGS & TRUTH GENERATION
% =========================================================================

% Simulation duration
T_total = 600;        % Total simulation time [s] (e.g. 10 minutes)
t = 0:IMU.dt:T_total; % Time vector [s]
N = numel(t);

fprintf('\n=== Generating Truth Trajectory ===\n');
fprintf(' Simulation time: %.1f s (%d samples)\n', T_total, N);

%-------------------------------------------------------------------------- 
% 3.1 TRUE ANGULAR RATE PROFILE (BODY FRAME)
%--------------------------------------------------------------------------
% Define a smooth, low-rate rotational profile representative of:
% - A spacecraft in fine-pointing mode with small residual rates
% - Useful to test IMU noise and bias behaviour without large dynamics

omega_true_b = zeros(3, N);        % [rad/s]

% Constant small roll rate (e.g. residual wheel imbalance)
omega_true_b(1, :) = deg2rad(0.03) * ones(1, N); % [rad/s]

% Slow sinusoidal pitch motion (e.g. thermal flexing / structural mode)
f_pitch = 1/300;                                            % [Hz]
omega_true_b(2, :) = deg2rad(0.02) * sin(2*pi*f_pitch * t); % [rad/s]

% Even slower yaw drift (e.g. long-term pointing drift)
f_yaw = 1/500;                                            % [Hz]
omega_true_b(3, :) = deg2rad(0.01) * cos(2*pi*f_yaw * t); % [rad/s]

%-------------------------------------------------------------------------- 
% 3.2 TRUE SPECIFIC FORCE PROFILE (BODY FRAME, OPTIONAL)
%--------------------------------------------------------------------------
% For pure attitude determination with star trackers + gyros, accelerometers
% are not strictly required (f_true_b can be set to []). Here we define a
% simple example assuming:
% - Spacecraft in LVLH-like frame
% - No translational dynamics modelled (for now)
%
% For now, we disable accelerometer usage by setting f_true_b = [].
% Later, this can be replaced by a full 6-DOF dynamics model.

use_accel_truth = true;

if use_accel_truth
    % Example: constant gravity-like specific force in body Z (e.g. ground test)
    g0 = 9.80665;                  % [m/s^2]
    f_true_b = repmat([0; 0; -g0], 1, N);
else
    f_true_b = [];                 % No accelerometer truth used in this phase
end

fprintf(' Maximum truth angular rates: [%.3f, %.3f, %.3f] deg/s\n', ...
    max(rad2deg(omega_true_b(1,:))), ...
    max(rad2deg(omega_true_b(2,:))), ...
    max(rad2deg(omega_true_b(3,:))));

%% ========================================================================
% 4. IMU MEASUREMENT GENERATION & QUICKLOOK PLOTS
% =========================================================================

fprintf('\n=== Generating Synthetic IMU Measurements ===\n');
imu_meas = generateIMUMeasurements(t, omega_true_b, f_true_b, IMU);
fprintf(' IMU measurement generation completed.\n');

fig_IMU = plotIMUResults(t, omega_true_b, f_true_b, imu_meas, use_accel_truth);
fprintf('\n=== IMU Simulation & Visualization Completed ===\n');

%% ========================================================================
% 5. VALIDATION AND ERROR ANALYSIS
% =========================================================================

fprintf('\n=== IMU Validation and Error Analysis ===\n');

% Gyro errors
gyro_err = imu_meas.gyro_meas_b - omega_true_b;      % [rad/s]
theta_err = cumsum(gyro_err, 2) * IMU.dt;            % [rad]

% 5.1 RMS rate error per axis
rms_rate_rad = sqrt(mean(gyro_err.^2, 2));           % [rad/s]
rms_rate_deg = rad2deg(rms_rate_rad);                % [deg/s]

fprintf('\n--- Gyro Rate RMS Error ---\n');
fprintf(' RMS(omega_x): %.4f deg/s\n', rms_rate_deg(1));
fprintf(' RMS(omega_y): %.4f deg/s\n', rms_rate_deg(2));
fprintf(' RMS(omega_z): %.4f deg/s\n', rms_rate_deg(3));

% 5.2 Orientation error growth (integrated small-angle error)
final_theta_deg = rad2deg(theta_err(:, end));

fprintf('\n--- Integrated Orientation Error (no correction) ---\n');
fprintf(' After %.1f s:\n', t(end));
fprintf('  Delta_theta_x: %.3f deg\n', final_theta_deg(1));
fprintf('  Delta_theta_y: %.3f deg\n', final_theta_deg(2));
fprintf('  Delta_theta_z: %.3f deg\n', final_theta_deg(3));

% 5.3 Gyro dynamic bias statistics
bias_dyn_deg_h = rad2deg(imu_meas.gyro_bias_dyn) * 3600; % [deg/h]

mean_bias_dyn = mean(bias_dyn_deg_h, 2);
std_bias_dyn  = std(bias_dyn_deg_h, 0, 2);

fprintf('\n--- Gyro Dynamic Bias Statistics ---\n');
fprintf(' Mean bias_x: %.3f deg/h, Std: %.3f deg/h\n', ...
    mean_bias_dyn(1), std_bias_dyn(1));
fprintf(' Mean bias_y: %.3f deg/h, Std: %.3f deg/h\n', ...
    mean_bias_dyn(2), std_bias_dyn(2));
fprintf(' Mean bias_z: %.3f deg/h, Std: %.3f deg/h\n', ...
    mean_bias_dyn(3), std_bias_dyn(3));

% 5.4 Accelerometer error statistics (if used)
if use_accel_truth && ~isempty(f_true_b) && isfield(imu_meas, 'accel_meas_b')
    accel_err = imu_meas.accel_meas_b - f_true_b;          % [m/s^2]
    rms_accel = sqrt(mean(accel_err.^2, 2));               % [m/s^2]

    fprintf('\n--- Accelerometer RMS Error ---\n');
    fprintf(' RMS(a_x): %.4e m/s^2\n', rms_accel(1));
    fprintf(' RMS(a_y): %.4e m/s^2\n', rms_accel(2));
    fprintf(' RMS(a_z): %.4e m/s^2\n', rms_accel(3));

    mean_accel_bias = mean(imu_meas.accel_bias_dyn, 2);
    std_accel_bias  = std(imu_meas.accel_bias_dyn, 0, 2);

    fprintf('\n--- Accelerometer Dynamic Bias Statistics ---\n');
    fprintf(' Mean bias_x: %.4e m/s^2, Std: %.4e m/s^2\n', ...
        mean_accel_bias(1), std_accel_bias(1));
    fprintf(' Mean bias_y: %.4e m/s^2, Std: %.4e m/s^2\n', ...
        mean_accel_bias(2), std_accel_bias(2));
    fprintf(' Mean bias_z: %.4e m/s^2, Std: %.4e m/s^2\n', ...
        mean_accel_bias(3), std_accel_bias(3));
end

fprintf('\n=== IMU Error Analysis Completed ===\n');