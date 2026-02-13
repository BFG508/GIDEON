%==========================================================================
% mainIMU.m - Inertial Measurement Unit Simulation
%==========================================================================
%
% DESCRIPTION:
%   Simulation framework for generating high-fidelity synthetic IMU data
%   (groscope and accelerometer) for AOCS/GNC development. This script
%   models realistic error sources including deterministic errors (bias, 
%   scale factor, misalignment) and stochastic noise processes (ARW/VRW,
%   bias instability, RRW).
%
% PURPOSE:
%   Demonstrates high-fidelity IMU error modeling for spacecraft AOCS/GNC
%   applications using realistic noise processes and deterministic errors.
%   Validates error models through Allan variance analysis and provides
%   filter design guidelines for EFK/EKF implementations.
%
% ALGORITHM:
%   The measurement model follows the standard equation for each sensor type:
%     u_meas = T_det * u_true + b_static + b_dynamic + noise
% 
%   where:
%     T_det    : M_SF * M_NonOrth * DCM_Mounting (Deterministic Transform) 
%     b_static : Turn-on bias                    (Constant per run)
%     b_dynamic: Dynamic bias drift              (Random Walk / Gauss-Markov process)
%     noise    : High-frequency white noise      (ARW/VRW)
% 
% WORKFLOW:
%   1. Define IMU hardware parameters (gyro/accel specifications) and
%      generate mounting misalignment and deterministic error matrices
%   2. Create truth trajectory (angular rates and specific forces)
%   3. Simulate IMU measurements with realistic error corruption
%   4. Validate error models: RMS statistics, bias drift analysis, 
%      Allan deviation analysis for noise characterization
%   5. Generate diagnostic plots and filter design recommendations
%
% KEY FEATURES:
%   - Space-grade IMU specifications (MEMS/FOG)
%   - Realistic error models: ARW/VRW, bias instability, RRW
%   - Deterministic errors: scale factor, non-orthogonality, mounting
%   - Allan deviation analysis for noise parameter extraction
%   - EKF filter tuning guidelines (Q and R matrices)
%   - Integration error growth quantification
%   - Publication-quality figures (saved as FIG, PNG, SVG)
%
% OUTPUTS:
%   - imuMeas:            Structure with corrupted measurements
%   - Error statistics:   RMS rate/accel errors, bias drift analysis
%   - Allan deviation:    Noise parameters (ARW, VRW, RRW, bias instability)
%   - Filter guidelines:  Recommended Q/R matrices for EFK
%   - Integration errors: Orientation/velocity drift without correction
%   - Diagnostic figures: 7 plots in Figures/ directory
%
% CONFIGURATION:
%   Edit Section 1 (initializeIMU.m) for:
%   - IMU hardware specifications (gyro/accel noise parameters)
%   - Mounting misalignment errors
%   - Sampling rate
%   - Deterministic error limits (bias, scale factor, non-orthogonality)
%
%   Edit Section 2 for:
%   - Simulation duration
%   - Truth trajectory profile (angular rates, specific forces)
%   - Dynamics scenarios (fine-pointing, slew, tumbling)
%==========================================================================

clear;
close all;
clc;

%% ========================================================================
% 1. IMU HARDWARE PARAMETERS (typical values from Space-Grade MEMS/FOG datasheets)
% =========================================================================

IMU = initializeIMU();

%% ========================================================================
% 2. SIMULATION SETTINGS & GROUND TRUTH GENERATION
% =========================================================================

% Simulation time settings
T_total = 60 * 45;    % Total simulation time [s]
t = 0:IMU.dt:T_total; % Time vector [s]
N = numel(t);

fprintf('\n=== Generating Truth Trajectory ===\n');
fprintf(' Simulation time: %.1f s (%d samples)\n', T_total, N);

% --- True Angular Rate Profile ---
% Define a smooth, low-rate rotational profile representative of:
%  - A spacecraft in fine-pointing mode with small residual rates
%  - Useful to test IMU noise and bias behaviour without large dynamics
    omegaTrue_body = zeros(3, N);        % [rad/s]
    
    % Constant small roll rate (e.g. residual wheel imbalance)
    omegaTrue_body(1, :) = deg2rad(0.03) * ones(1, N);           % [rad/s]
    
    % Slow sinusoidal pitch motion (e.g. thermal flexing / structural mode)
    fPitch = 1/300;                                              % [Hz]
    omegaTrue_body(2, :) = deg2rad(0.02) * sin(2*pi*fPitch * t); % [rad/s]
    
    % Even slower yaw drift (e.g. long-term pointing drift)
    fYaw = 1/500;                                                % [Hz]
    omegaTrue_body(3, :) = deg2rad(0.01) * cos(2*pi*fYaw * t);   % [rad/s]
    
    fprintf(' Maximum truth angular rates: [%.3f, %.3f, %.3f] deg/s\n', ...
        max(rad2deg(omegaTrue_body(1,:))), ...
        max(rad2deg(omegaTrue_body(2,:))), ...
        max(rad2deg(omegaTrue_body(3,:))));

% --- True Specific Force Profile ---  
% Define a smooth profile representative of:
%  - Spacecraft in LVLH-like frame
%  - No translational dynamics modelled
    g0             = 9.80665;                   % [m/s^2]
    forceTrue_body = repmat([0; 0; -g0], 1, N); % [m/s^2]

%% ========================================================================
% 3. GENERATE SYNTHETIC IMU MEASUREMENTS
% =========================================================================

meas = generateIMUMeasurements(t, omegaTrue_body, forceTrue_body, IMU);

%% ========================================================================
% 4. VALIDATION AND ERROR ANALYSIS
% =========================================================================

fprintf('\n=== IMU Validation and Error Analysis ===\n');

% Gyro statistic errors
gyroErr  = meas.gyro.omegaBody - omegaTrue_body; % [rad/s]
thetaErr = cumsum(gyroErr, 2) * IMU.dt;          % [rad]

    % RMS rate error per axis
    rateRMS = rad2deg(sqrt(mean(gyroErr.^2, 2))); % [deg/s]
    
    fprintf('\n--- Gyro Rate RMS Error ---\n');
    fprintf(' RMS(ω_x): %.4f deg/s\n', rateRMS(1));
    fprintf(' RMS(ω_y): %.4f deg/s\n', rateRMS(2));
    fprintf(' RMS(ω_z): %.4f deg/s\n', rateRMS(3));
    
    % Orientation error growth (integrated small-angle error)
    thetaFinal = rad2deg(thetaErr(:, end));
    
    fprintf('\n--- Integrated Orientation Error (no correction) ---\n');
    fprintf(' After %.1f s:\n', t(end));
    fprintf('  Δθ_x: %.3f deg\n', thetaFinal(1));
    fprintf('  Δθ_y: %.3f deg\n', thetaFinal(2));
    fprintf('  Δθ_z: %.3f deg\n', thetaFinal(3));
    
    % Gyro dynamic bias statistics
    gyroMeanBiasDyn = mean(rad2deg(meas.gyro.biasDyn) * 3600, 2);
     gyroStdBiasDyn =  std(rad2deg(meas.gyro.biasDyn) * 3600, 0, 2);
    
    fprintf('\n--- Gyro Dynamic Bias Statistics ---\n');
    fprintf(' x) Mean bias: %.3f deg/h, Std: %.3f deg/h\n', ...
        gyroMeanBiasDyn(1), gyroStdBiasDyn(1));
    fprintf(' y) Mean bias: %.3f deg/h, Std: %.3f deg/h\n', ...
        gyroMeanBiasDyn(2), gyroStdBiasDyn(2));
    fprintf(' z) Mean bias: %.3f deg/h, Std: %.3f deg/h\n', ...
        gyroMeanBiasDyn(3), gyroStdBiasDyn(3));

% Accelerometer statistic errors
accelErr  = meas.accel.forceBody - forceTrue_body; % [m/s²]

    % RMS rate error per axis
    forceRMS = rad2deg(sqrt(mean(accelErr.^2, 2))); % [m/s²]
    
    fprintf('\n--- Accel Force RMS Error ---\n');
    fprintf(' RMS(f_x): %.4f m/s²\n', forceRMS(1));
    fprintf(' RMS(f_y): %.4f m/s²\n', forceRMS(2));
    fprintf(' RMS(f_z): %.4f m/s²\n', forceRMS(3));
    
    % Gyro dynamic bias statistics
    accelMeanBiasDyn = mean(meas.accel.biasDyn, 2);
     accelStdBiasDyn =  std(meas.accel.biasDyn, 0, 2);
    
    fprintf('\n--- Accel Dynamic Bias Statistics ---\n');
    fprintf(' x) Mean bias: %.3f m/s², Std: %.3f m/s²\n', ...
        accelMeanBiasDyn(1), accelStdBiasDyn(1));
    fprintf(' y) Mean bias: %.3f m/s², Std: %.3f m/s²\n', ...
        accelMeanBiasDyn(2), accelStdBiasDyn(2));
    fprintf(' z) Mean bias: %.3f m/s², Std: %.3f m/s²\n', ...
        accelMeanBiasDyn(3), accelStdBiasDyn(3));

% Allan Deviation Analysis (Noise characterization)
fprintf('\n--- Allan Deviation Analysis (X-axis representative) ---\n');

    % Compute Allan Deviation for gyro (X-axis)
    [gyroTau, gyroADEV] = computeAllanDeviation(gyroErr(1,:), IMU.dt);
    
    % Initialize empirical parameters
    ARWEmpirical = NaN;
    RRWEmpirical = NaN;
    
    % Find characteristic points in Allan Deviation curve
    [gyroADEVMin, gyroADEVMin_idx] = min(gyroADEV);
    gyroTauMin = gyroTau(gyroADEVMin_idx);
    
    % Estimate Angle Random Walk (ARW) from left slope (tau = 1 s)
    gyroTau_1_idx = find(gyroTau >= 1.0, 1, 'first');
    if ~isempty(gyroTau_1_idx) && gyroTau_1_idx > 1
        ARWEmpirical = rad2deg(gyroADEV(gyroTau_1_idx)) * sqrt(3600);
    end
    
    % Bias Instability is the minimum value
    biasInstEmpirical = rad2deg(gyroADEVMin) * 3600;
    
    % Rate Random Walk (RRW) from right slope (tau = 10 s)
    gyroTau_10_idx = find(gyroTau >= 10, 1, 'first');
    if ~isempty(gyroTau_10_idx)
        RRWEmpirical = gyroADEV(gyroTau_10_idx) / sqrt(gyroTau(gyroTau_10_idx));
        RRWEmpirical = rad2deg(RRWEmpirical) * 3600 * sqrt(3600);
    end
    
    fprintf('\n Gyro Noise Parameters (Empirical from ADEV):\n');
    if ~isnan(ARWEmpirical)
        fprintf('  Angle Random Walk (ARW):        %.4f deg/√h  (Spec: %.4f)\n', ...
                ARWEmpirical, IMU.gyro.ARW);
    end
    fprintf('  Bias Instability:               %.4f deg/h   (Spec: %.4f)\n', ...
            biasInstEmpirical, IMU.gyro.biasInstability);
    if ~isnan(RRWEmpirical)
        fprintf('  Rate Random Walk (RRW):         %.4f deg/h/√h (Spec: %.4f)\n', ...
                RRWEmpirical, IMU.gyro.RRW);
    end
    fprintf('  Optimal averaging time (τ_min): %.2f s\n', gyroTauMin);
    

    % Compute Allan Deviation for accelerometer (X-axis)
    [accelTau, accelADEV] = computeAllanDeviation(accelErr(1,:), IMU.dt);
    
    % Initialize empirical parameterss
    VRWEmpirical      = NaN;
    accelARWEmpirical = NaN;
    
    % Find characteristic points in Allan Deviation curve
    [accelADEVMin, accelADEVMin_idx] = min(accelADEV);
    accelTauMin = accelTau(accelADEVMin_idx);
    
    % Estimate Velocity Random Walk (VRW) from left slope (tau = 1 s)
    idx_1s_a = find(accelTau >= 1.0, 1, 'first');
    if ~isempty(idx_1s_a) && idx_1s_a > 1
        VRWEmpirical = accelADEV(idx_1s_a) / sqrt(2) / (1e-6 * g0);
    end
    
    % Bias Instability is the minimum value
    accelBiasInstEmpirical = accelADEVMin / (1e-6 * 9.80665);
    
    % Acceleration Random Walk  from right slope (tau = 10 s)
    accelTau_10_idx = find(accelTau >= 10, 1, 'first');
    if ~isempty(accelTau_10_idx)
        accelARWEmpirical = accelADEV(accelTau_10_idx) / sqrt(accelTau(accelTau_10_idx)) * sqrt(3600);
    end
    
    fprintf('\n Accel Noise Parameters (Empirical from ADEV):\n');
    if ~isnan(VRWEmpirical)
        fprintf('  Velocity Random Walk (VRW):     %.2f μg/√Hz  (Spec: %.2f)\n', ...
                VRWEmpirical, IMU.accel.VRW);
    end
    fprintf('  Bias Instability:               %.2f μg      (Spec: %.2f)\n', ...
            accelBiasInstEmpirical, IMU.accel.biasInstability);
    if ~isnan(accelARWEmpirical)
        fprintf('  Accel Random Walk (ARW):        %.4f m/s/√h (Spec: %.4f)\n', ...
                accelARWEmpirical, IMU.accel.ARW);
    end
    fprintf('  Optimal averaging time (τ_min): %.2f s\n', accelTauMin);
    
% Implications for EFK
fprintf('\n--- Implications for EFK Design ---\n');

fprintf(' 1. MEASUREMENT NOISE (Matrix R):\n');
fprintf('    Gyro white noise dominates at high rates (< %.1f s).\n', gyroTauMin);
if ~isnan(ARWEmpirical)
    fprintf('    → Use R_gyro based on ARW: ~%.2e rad²/s² per axis\n', ...
            (deg2rad(ARWEmpirical) / sqrt(3600))^2);
end
fprintf('    Accel white noise dominates at high rates (< %.1f s).\n', accelTauMin);
if ~isnan(VRWEmpirical)
    fprintf('    → Use R_accel based on VRW: ~%.2e m²/s⁴ per axis\n\n', ...
            (VRWEmpirical)^2);
end

fprintf(' 2. PROCESS NOISE (Matrix Q - Bias Random Walk):\n');
fprintf('    Gyro bias drifts significantly after ~%.1f s.\n', 3*gyroTauMin);
if ~isnan(RRWEmpirical)
    fprintf('    → Model bias drift with qBiasGyro ~%.2e (rad/s)²/s\n', ...
            (RRWEmpirical * sqrt(IMU.dt))^2);
end
fprintf('    Accel bias drifts significantly after ~%.1f s.\n', 3*accelTauMin);
if ~isnan(accelARWEmpirical)
    fprintf('    → Model bias drift with qBiasAccel ~%.2e (m/s²)²/s\n\n', ...
            (accelARWEmpirical * sqrt(IMU.dt))^2);
end

fprintf(' 3. OPTIMAL UPDATE RATE:\n');
fprintf('    External sensor updates (STR, MAG) should occur every ~%.1f-%.1f s\n', ...
        gyroTauMin, 2*gyroTauMin);
fprintf('    to correct bias drift before random walk dominates.\n\n');

fprintf(' 4. INTEGRATION ERROR GROWTH (Open-Loop):\n');
fprintf('    Gyro: Orientation error grows as √t after %.1f s\n', 3*gyroTauMin);
fprintf('    → Δθ ~ %.2f deg after %.0f s without correction\n', ...
        rad2deg(thetaErr(1, end)), t(end));
fprintf('    Accel: Velocity error grows as √t after %.1f s\n', 3*accelTauMin);
vel_err_final = sum(accelErr(1,:)) * IMU.dt;
fprintf('    → Δv ~ %.2f m/s after %.0f s without correction\n', ...
      vel_err_final, t(end));

%% ========================================================================
% 5. RESULTS VISUALIZATION
% =========================================================================

% Generate comprehensive validation plots
saveFlag = 1;
plotIMUResults(t, omegaTrue_body, forceTrue_body, meas, saveFlag);

fprintf('\n=== All tasks completed successfully ===\n');