%==========================================================================
% mainIMU.m - Inertial Measurement Unit (IMU) Simulation & Calibration
%==========================================================================
%
% DESCRIPTION:
%   Simulation framework for generating high-fidelity synthetic IMU data
%   (Gyroscopes and Accelerometers) for AOCS/GNC development. This script
%   models realistic error sources including deterministic errors (bias, 
%   scale factor, misalignment) and stochastic noise processes (Angle/Velocity
%   Random Walk, Bias Instability, Rate Random Walk).
%
% PURPOSE:
%   1. Define hardware specifications for typical space-grade IMUs.
%   2. Generate "Truth" trajectory (angular velocity & specific force).
%   3. Simulate sensor measurement corruption with independent error models.
%   4. Validate IMU error models through statistical analysis and plots.
%
% ALGORITHM:
%   The measurement model follows the standard equation for each sensor type:
%     u_meas = T_det * u_true + b_static + b_dynamic + noise
%   
%   where:
%     T_det     = M_SF * M_NonOrth * DCM_Mounting (Deterministic Transform)
%     M_SF      = Scale Factor error matrix (Sensor specific)
%     M_NonOrth = Non-orthogonality matrix (Sensor specific)
%     b_static  = Turn-on bias (Constant per run)
%     b_dynamic = Dynamic bias drift (Random Walk / Gauss-Markov process)
%     noise     = High-frequency white noise (ARW for Gyro, VRW for Accel)
%
%   Note: The accelerometer and gyroscope have independent internal error
%   matrices (Scale Factor & Non-orthogonality) but share the external
%   mounting misalignment (DCM_Mounting).
%==========================================================================

clear;
close all;
clc;

%% ========================================================================
% 1. IMU HARDWARE PARAMETERS (typical values from Space-Grade MEMS/FOG datasheets)
% =========================================================================

% --- General Configuration ---
IMU.rate = 120;      % Sampling frequency [Hz] (Typical: 10-1000 Hz)
IMU.dt = 1/IMU.rate; % Sampling time step [s]

% Mounting Misalignment: Error aligning the IMU case to the Body frame.
IMU.mountingError = [0.1; 0.1; 0.1]; % [deg] Roll, Pitch, Yaw mounting error

% --- Gyroscope Parameters (Angular Rate) ---
% Representative of a high-end MEMS or low-end FOG
    % 1. Stochastic Noise (Random processes)
    % ARW: High-frequency noise (white noise density). Determines short-term precision.
    IMU.gyro.ARW = 0.05;              % Angle Random Walk [deg/sqrt(h)]
    
    % Bias Instability: Long-term drift (flicker noise 1/f floor).
    IMU.gyro.biasInstability = 1.0;   % [deg/h]
    
    % RRW: Rate Random Walk. Models the "wandering" of the bias over time.
    IMU.gyro.RRW = 0.1;               % [deg/h/sqrt(h)] (bias drift rate)
    
    % 2. Deterministic Errors (Calibration residuals)
    % Static Bias: Constant offset per power cycle (Turn-on to turn-on repeatability).
    IMU.gyro.biasStaticLimit = 10.0;  % [deg/h] (1-sigma)
    
    % Scale Factor Error: Sensitivity error (linear gain mismatch).
    IMU.gyro.SF = 300;                % [ppm]
    
    % 3. Geometric Errors
    % Internal Misalignment: Non-orthogonality of the sensor triad axes.
    IMU.gyro.nonorthogonality = 0.05; % [deg] (1-sigma)

% --- Accelerometer Parameters (Specific Force) ---
% Required for g-sensitivity corrections or Inertial Navigation
    % 1. Stochastic Noise
    IMU.accel.VRW = 50;                % Velocity Random Walk [μg/sqrt(Hz)]
    IMU.accel.biasInstability = 20;    % [μg]
    IMU.accel.ARW = 0.05;              % Acceleration Random Walk [m/s/sqrt(h)] (for bias drift)
      
    % 2. Deterministic Errors  
    IMU.accel.biasStaticLimit = 1.0;   % [mg]
    IMU.accel.SF = 300;                % [ppm]
    
    % 3. Geometric Errors
    IMU.accel.nonorthogonality = 0.05; % [deg]

% --- Derived Matrices & Pre-Allocation ---
fprintf('\n=== Initializing IMU Model ===\n');

    % Construct Mounting Misalignment Matrix (DCM Body -> Case)
    % This matrix is shared between gyro and accel (common mechanical mounting).
    roll_err  = deg2rad(IMU.mountingError(1));
    pitch_err = deg2rad(IMU.mountingError(2));
    yaw_err   = deg2rad(IMU.mountingError(3));
    
    IMU.DCM_mounting = angle2dcm(yaw_err, pitch_err, roll_err, 'ZYX');
    
    % --- Gyroscope Deterministic Transform ---
    % Construct Scale Factor Matrix (Diagonal)
    % S = diag([sf_x, sf_y, sf_z]) where sf = scale_error_ppm * 1e-6
    sf_err_gyro = (IMU.gyro.SF * 1e-6) * randn(3,1); 
    IMU.gyro.M_SF = eye(3) + diag(sf_err_gyro);
    
    fprintf(' Gyro Scale Factor Errors: [%.1f, %.1f, %.1f] ppm\n', ...
        sf_err_gyro(1)*1e6, sf_err_gyro(2)*1e6, sf_err_gyro(3)*1e6);
    
    % Construct Non-Orthogonality Matrix
    % Models small deviations of the axes from a perfect 90° triad.
    % Upper triangular approximation for small angles.
    nonortho_err_gyro = deg2rad(IMU.gyro.nonorthogonality) * randn(3,1);
    IMU.gyro.M_nonorth = [                    1, -nonortho_err_gyro(3),  nonortho_err_gyro(2);
                          -nonortho_err_gyro(3),                     1, -nonortho_err_gyro(1);
                           nonortho_err_gyro(2), -nonortho_err_gyro(1),                    1];
    
    % Total deterministic transform (Body True -> Gyro Sensor Frame)
    % T_gyro = M_Scale * M_NonOrth * DCM_Mounting
    IMU.gyro.T_Deterministic = IMU.gyro.M_SF * IMU.gyro.M_nonorth * IMU.DCM_mounting;
    
    % Generate Static Turn-on Bias (Constant for this simulation run)
    IMU.gyro.biasStatic = deg2rad(IMU.gyro.biasStaticLimit / 3600) * randn(3,1);
    
    fprintf(' Gyro Static Turn-on Bias: [%.4f, %.4f, %.4f] deg/h\n', ...
        rad2deg(IMU.gyro.biasStatic(1))*3600, ...
        rad2deg(IMU.gyro.biasStatic(2))*3600, ...
        rad2deg(IMU.gyro.biasStatic(3))*3600);
    
    % --- Accelerometer Deterministic Transform ---
    % Construct Scale Factor Matrix (independent from gyro)
    sf_err_accel = (IMU.accel.SF * 1e-6) * randn(3,1);
    IMU.accel.M_SF = eye(3) + diag(sf_err_accel);
    
    fprintf(' Accel Scale Factor Errors: [%.1f, %.1f, %.1f] ppm\n', ...
        sf_err_accel(1)*1e6, sf_err_accel(2)*1e6, sf_err_accel(3)*1e6);
    
    % Construct Non-Orthogonality Matrix (independent from gyro)
    nonortho_err_accel = deg2rad(IMU.accel.nonorthogonality) * randn(3,1);
    IMU.accel.M_nonorth = [                     1, -nonortho_err_accel(3),  nonortho_err_accel(2);
                           -nonortho_err_accel(3),                      1, -nonortho_err_accel(1);
                            nonortho_err_accel(2), -nonortho_err_accel(1),                     1];
    
    % Total deterministic transform (Body True -> Accel Sensor Frame)
    % Note: DCM_mounting is shared (common mechanical housing)
    IMU.accel.T_Deterministic = IMU.accel.M_SF * IMU.accel.M_nonorth * IMU.DCM_mounting;
    
    % Generate Static Turn-on Bias (Constant for this simulation run)
    g0 = 9.80665; % [m/s^2]
    IMU.accel.biasStatic = (IMU.accel.biasStaticLimit * 1e-3 * g0) * randn(3,1); % [mg] -> [m/s^2]
    
    fprintf(' Accel Static Turn-on Bias: [%.4f, %.4f, %.4f] mg\n', ...
        IMU.accel.biasStatic(1)/(1e-3*g0), ...
        IMU.accel.biasStatic(2)/(1e-3*g0), ...
        IMU.accel.biasStatic(3)/(1e-3*g0));

fprintf('=== IMU Initialization Complete ===\n');

%% ========================================================================
% 2. SIMULATION SETTINGS & GROUND TRUTH GENERATION
% =========================================================================

% Simulation time settings
T_total = 60*40;      % Total simulation time [s]
t = 0:IMU.dt:T_total; % Time vector [s]
N = numel(t);

fprintf('\n=== Generating Truth Trajectory ===\n');
fprintf(' Simulation time: %.1f s (%d samples)\n', T_total, N);

% --- True Angular Rate Profile (Body frame) ---
% Define a smooth, low-rate rotational profile representative of:
%  - A spacecraft in fine-pointing mode with small residual rates
%  - Useful to test IMU noise and bias behaviour without large dynamics
    omega_true_b = zeros(3, N);        % [rad/s]
    
    % Constant small roll rate (e.g. residual wheel imbalance)
    omega_true_b(1, :) = deg2rad(0.03) * ones(1, N); % [rad/s]
    
    % Slow sinusoidal pitch motion (e.g. thermal flexing / structural mode)
    f_pitch = 1/300;                                            % [Hz]
    omega_true_b(2, :) = deg2rad(0.02) * sin(2*pi*f_pitch * t); % [rad/s]
    
    % Even slower yaw drift (e.g. long-term pointing drift)
    f_yaw = 1/500;                                            % [Hz]
    omega_true_b(3, :) = deg2rad(0.01) * cos(2*pi*f_yaw * t); % [rad/s]

% --- True Specific Force Profile (Body frame) ---  
% For pure attitude determination with star trackers + gyros, accelerometers
% are not strictly required. Here we define a simple example assuming:
%  - Spacecraft in LVLH-like frame
%  - No translational dynamics modelled
    g0 = 9.80665;                         % [m/s^2]
    f_true_b = repmat([0; 0; -g0], 1, N);
    
    fprintf(' Maximum truth angular rates: [%.3f, %.3f, %.3f] deg/s\n', ...
        max(rad2deg(omega_true_b(1,:))), ...
        max(rad2deg(omega_true_b(2,:))), ...
        max(rad2deg(omega_true_b(3,:))));

%% ========================================================================
% 3. GENERATE SYNTHETIC IMU MEASUREMENTS
% =========================================================================

fprintf('\n=== Generating Synthetic IMU Measurements ===\n');
imu_meas = generateIMUMeasurements(t, omega_true_b, f_true_b, IMU);

%% ========================================================================
% 4. VALIDATION AND ERROR ANALYSIS
% =========================================================================

fprintf('\n=== IMU Validation and Error Analysis ===\n');

% Gyro errors
gyro_err  = imu_meas.gyro_meas_b - omega_true_b; % [rad/s]
theta_err = cumsum(gyro_err, 2) * IMU.dt;        % [rad]

% RMS rate error per axis
rms_rate_rad = sqrt(mean(gyro_err.^2, 2));       % [rad/s]
rms_rate_deg = rad2deg(rms_rate_rad);            % [deg/s]

fprintf('\n--- Gyro Rate RMS Error ---\n');
fprintf(' RMS(ω_x): %.4f deg/s\n', rms_rate_deg(1));
fprintf(' RMS(ω_y): %.4f deg/s\n', rms_rate_deg(2));
fprintf(' RMS(ω_z): %.4f deg/s\n', rms_rate_deg(3));

% Orientation error growth (integrated small-angle error)
final_theta_deg = rad2deg(theta_err(:, end));

fprintf('\n--- Integrated Orientation Error (no correction) ---\n');
fprintf(' After %.1f s:\n', t(end));
fprintf('  Δθ_x: %.3f deg\n', final_theta_deg(1));
fprintf('  Δθ_y: %.3f deg\n', final_theta_deg(2));
fprintf('  Δθ_z: %.3f deg\n', final_theta_deg(3));

% Gyro dynamic bias statistics
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

% Accelerometer error statistics
accel_err = imu_meas.accel_meas_b - f_true_b;          % [m/s²]
rms_accel = sqrt(mean(accel_err.^2, 2));               % [m/s²]

fprintf('\n--- Accelerometer RMS Error ---\n');
fprintf(' RMS(a_x): %.4e m/s²\n', rms_accel(1));
fprintf(' RMS(a_y): %.4e m/s²\n', rms_accel(2));
fprintf(' RMS(a_z): %.4e m/s²\n', rms_accel(3));

mean_accel_bias = mean(imu_meas.accel_bias_dyn, 2);
std_accel_bias  = std(imu_meas.accel_bias_dyn, 0, 2);

fprintf('\n--- Accelerometer Dynamic Bias Statistics ---\n');
fprintf(' Mean bias_x: %.4e m/s², Std: %.4e m/s²\n', ...
    mean_accel_bias(1), std_accel_bias(1));
fprintf(' Mean bias_y: %.4e m/s², Std: %.4e m/s²\n', ...
    mean_accel_bias(2), std_accel_bias(2));
fprintf(' Mean bias_z: %.4e m/s², Std: %.4e m/s²\n', ...
    mean_accel_bias(3), std_accel_bias(3));

% -------------------------------------------------------------------------
% Allan Deviation Analysis (Noise Characterization)
% -------------------------------------------------------------------------
    fprintf('\n--- Allan Deviation Analysis (X-axis representative) ---\n');
    
    % Compute Allan Deviation for Gyro (X-axis)
    [taus_g, sigma_allan_g] = computeAllanDeviation(gyro_err(1,:), IMU.dt);
    
    % Initialize empirical parameters (in case extraction fails)
    arw_empirical_deg_sqrt_h   = NaN;
    rrw_empirical_deg_h_sqrt_h = NaN;
    rrw_empirical_rad_s2       = NaN;
    
    % Find characteristic points in Allan Deviation curve
    [sigma_min_g, idx_min_g] = min(sigma_allan_g);
    tau_min_g = taus_g(idx_min_g);
    
    % Estimate Angle Random Walk (ARW) from left slope (tau = 1 s)
    idx_1s_g = find(taus_g >= 1.0, 1, 'first');
    if ~isempty(idx_1s_g) && idx_1s_g > 1
        arw_empirical_rad_s = sigma_allan_g(idx_1s_g);
        arw_empirical_deg_sqrt_h = rad2deg(arw_empirical_rad_s) * sqrt(3600);
    end
    
    % Bias Instability is the minimum value
    bias_inst_empirical_deg_h = rad2deg(sigma_min_g) * 3600;
    
    % Rate Random Walk (RRW) from right slope (tau = 10 s)
    idx_10s = find(taus_g >= 10, 1, 'first');
    if ~isempty(idx_10s)
        rrw_empirical_rad_s2 = sigma_allan_g(idx_10s) / sqrt(taus_g(idx_10s));
        rrw_empirical_deg_h_sqrt_h = rad2deg(rrw_empirical_rad_s2) * 3600 * sqrt(3600);
    end
    
    fprintf('\n GYRO Noise Parameters (Empirical from Allan Deviation):\n');
    if ~isnan(arw_empirical_deg_sqrt_h)
        fprintf('  Angle Random Walk (ARW):        %.4f deg/√h  (Spec: %.4f)\n', ...
                arw_empirical_deg_sqrt_h, IMU.gyro.ARW);
    end
    fprintf('  Bias Instability:               %.4f deg/h   (Spec: %.4f)\n', ...
            bias_inst_empirical_deg_h, IMU.gyro.biasInstability);
    if ~isnan(rrw_empirical_deg_h_sqrt_h)
        fprintf('  Rate Random Walk (RRW):         %.4f deg/h/√h (Spec: %.4f)\n', ...
                rrw_empirical_deg_h_sqrt_h, IMU.gyro.RRW);
    end
    fprintf('  Optimal averaging time (τ_min): %.2f s\n', tau_min_g);
    
    % Compute Allan Deviation for Accelerometer (X-axis)
    [taus_a, sigma_allan_a] = computeAllanDeviation(accel_err(1,:), IMU.dt);
    
    % Initialize
    vrw_empirical_ug_sqrt_Hz    = NaN;
    arw_a_empirical_m_s_sqrt_h  = NaN;
    arw_a_empirical_m_s2_sqrt_s = NaN;
    
    % Find characteristic points
    [sigma_min_a, idx_min_a] = min(sigma_allan_a);
    tau_min_a = taus_a(idx_min_a);
    
    % Estimate Velocity Random Walk (VRW) at tau = 1 s
    idx_1s_a = find(taus_a >= 1.0, 1, 'first');
    if ~isempty(idx_1s_a) && idx_1s_a > 1
        vrw_empirical_m_s2_sqrt_Hz = sigma_allan_a(idx_1s_a) / sqrt(2);
        vrw_empirical_ug_sqrt_Hz = vrw_empirical_m_s2_sqrt_Hz / (1e-6 * 9.80665);
    end
    
    % Bias Instability
    bias_inst_a_empirical_ug = sigma_min_a / (1e-6 * 9.80665);
    
    % Acceleration Random Walk from right slope
    idx_10s_a = find(taus_a >= 10, 1, 'first');
    if ~isempty(idx_10s_a)
        arw_a_empirical_m_s2_sqrt_s = sigma_allan_a(idx_10s_a) / sqrt(taus_a(idx_10s_a));
        arw_a_empirical_m_s_sqrt_h = arw_a_empirical_m_s2_sqrt_s * sqrt(3600);
    end
    
    fprintf('\n ACCELEROMETER Noise Parameters (Empirical from Allan Deviation):\n');
    if ~isnan(vrw_empirical_ug_sqrt_Hz)
        fprintf('  Velocity Random Walk (VRW):     %.2f μg/√Hz  (Spec: %.2f)\n', ...
                vrw_empirical_ug_sqrt_Hz, IMU.accel.VRW);
    end
    fprintf('  Bias Instability:               %.2f μg      (Spec: %.2f)\n', ...
            bias_inst_a_empirical_ug, IMU.accel.biasInstability);
    if ~isnan(arw_a_empirical_m_s_sqrt_h)
        fprintf('  Accel Random Walk (ARW):        %.4f m/s/√h (Spec: %.4f)\n', ...
                arw_a_empirical_m_s_sqrt_h, IMU.accel.ARW);
    end
    fprintf('  Optimal averaging time (τ_min): %.2f s\n', tau_min_a);
    
    % -------------------------------------------------------------------------
    % Implications for MEKF Filter Design
    % -------------------------------------------------------------------------
    fprintf('\n--- Implications for MEKF Design ---\n');
    fprintf(' Based on Allan Deviation analysis:\n\n');
    
    fprintf(' 1. MEASUREMENT NOISE (Matrix R):\n');
    fprintf('    Gyro white noise dominates at high rates (< %.1f s).\n', tau_min_g);
    if ~isnan(arw_empirical_deg_sqrt_h)
        fprintf('    → Use R_gyro based on ARW: ~%.2e rad²/s² per axis\n', ...
                (deg2rad(arw_empirical_deg_sqrt_h) / sqrt(3600))^2);
    end
    fprintf('    Accel white noise dominates at high rates (< %.1f s).\n', tau_min_a);
    if ~isnan(vrw_empirical_ug_sqrt_Hz)
        fprintf('    → Use R_accel based on VRW: ~%.2e m²/s⁴ per axis\n\n', ...
                (vrw_empirical_m_s2_sqrt_Hz)^2);
    end
    
    fprintf(' 2. PROCESS NOISE (Matrix Q - Bias Random Walk):\n');
    fprintf('    Gyro bias drifts significantly after ~%.1f s.\n', 3*tau_min_g);
    if ~isnan(rrw_empirical_rad_s2)
        fprintf('    → Model bias drift with q_bias_gyro ~%.2e (rad/s)²/s\n', ...
                (rrw_empirical_rad_s2 * sqrt(IMU.dt))^2);
    end
    fprintf('    Accel bias drifts significantly after ~%.1f s.\n', 3*tau_min_a);
    if ~isnan(arw_a_empirical_m_s2_sqrt_s)
        fprintf('    → Model bias drift with q_bias_accel ~%.2e (m/s²)²/s\n\n', ...
                (arw_a_empirical_m_s2_sqrt_s * sqrt(IMU.dt))^2);
    end
    
    fprintf(' 3. OPTIMAL UPDATE RATE:\n');
    fprintf('    External sensor updates (STR, MAG) should occur every ~%.1f-%.1f s\n', ...
            tau_min_g, 2*tau_min_g);
    fprintf('    to correct bias drift before random walk dominates.\n\n');
    
    fprintf(' 4. INTEGRATION ERROR GROWTH (Open-Loop):\n');
    fprintf('    Gyro: Orientation error grows as √t after %.1f s\n', 3*tau_min_g);
    fprintf('    → Δθ ~ %.2f deg after %.0f s without correction\n', ...
            rad2deg(theta_err(1, end)), t(end));
    fprintf('    Accel: Velocity error grows as √t after %.1f s\n', 3*tau_min_a);
    vel_err_final = sum(accel_err(1,:)) * IMU.dt;
    fprintf('    → Δv ~ %.2f m/s after %.0f s without correction\n', ...
          vel_err_final, t(end));
    
fprintf('\n=== IMU Error Analysis Completed ===\n');

%% ========================================================================
% 5. RESULTS VISUALIZATION
% =========================================================================

% Generate comprehensive validation plots
saveFlag = 1;
plotIMUResults(t, omega_true_b, f_true_b, imu_meas, saveFlag);

fprintf('\n=== All tasks completed successfully ===\n');