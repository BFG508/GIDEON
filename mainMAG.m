%==========================================================================
% mainMAG.m - Magnetometer Simulation
%==========================================================================
%
% DESCRIPTION:
%   Simulation framework for generating high-fidelity synthetic 3-axis
%   magnetometer data for spacecraft AOCS/GNC development. This script
%   models environmental magnetic fields (IGRF) and realistic sensor error
%   sources including deterministic errors (hard iron, soft iron, 
%   non-orthogonality) and stochastic noise processes (white noise, bias 
%   instability).
%
% PURPOSE:
%   Demonstrates high-fidelity magnetometer error modeling for spacecraft
%   AOCS/GNC applications using realistic orbital dynamics, Earth magnetic
%   field models, and sensor noise. Validates error models through residual
%   analysis and provides ground truth for calibration algorithm testing.
%
% ALGORITHM:
%   The measurement model follows the standard equation:
%     B_meas = M_det * B_true + b_hard_iron + b_dynamic + noise
% 
%   where:
%     M_det       : M_SoftIron * M_NonOrth * DCM_Mounting (Deterministic Transform)
%     b_hard_iron : Static hard iron bias                 (Constant per run)
%     b_dynamic   : Dynamic bias drift                    (Gauss-Markov process)
%     noise       : High-frequency white noise            (Gaussian)
% 
% WORKFLOW:
%   1. Define MAG hardware parameters (sensor specifications) and generate
%      mounting misalignment and deterministic error matrices
%   2. Propagate orbital trajectory (J2 + Drag + SRP perturbations)
%   3. Compute Earth magnetic field using IGRF-13 model
%   4. Generate truth attitude trajectory and rotate field to body frame
%   5. Simulate MAG measurements with realistic error corruption
%   6. Validate error models: RMS errors, noise statistics, angular errors
%   7. Generate diagnostic plots for calibration validation
%
% KEY FEATURES:
%   - Space-grade magnetometer specifications (Fluxgate/AMR)
%   - Realistic error models: hard/soft iron, bias instability, white noise
%   - Deterministic errors: scale factor, non-orthogonality, mounting
%   - IGRF-13 geomagnetic field model with orbital propagation
%   - Attitude dynamics integration for body-frame field
%   - Calibration validation: 3D field clouds, projection plots
%   - Publication-quality figures (saved as FIG, PNG, SVG)
%
% OUTPUTS:
%   - meas:               Structure with corrupted measurements
%   - Error statistics:   RMS field errors, noise verification, angular errors
%   - Ground truth:       Orbital position, attitude, magnetic field vectors
%   - Calibration data:   3D point clouds for ellipsoid fitting algorithms
%   - Diagnostic figures: 4 plots in Figures/MAG/ directory
%
% CONFIGURATION:
%   Edit Section 1 (initializeMAG.m) for:
%   - MAG hardware specifications (noise density, bias instability)
%   - Mounting misalignment errors
%   - Sampling rate
%   - Deterministic error limits (hard iron, soft iron, non-orthogonality)
%
%   Edit Section 2 for:
%   - Simulation duration and epoch
%   - Orbital elements (SMA, ECC, INC, RAAN, AOP, TA)
%   - Spacecraft parameters (mass, drag/SRP areas)
%   - Attitude dynamics profile (angular rates)
%==========================================================================

clear;
close all;
clc;

%% ========================================================================
% 1. MAG HARDWARE PARAMETERS (typical values from Space-Grade Fluxgate/AMR datasheets)
% =========================================================================

MAG = initializeMAG();

%% ========================================================================
% 2. SIMULATION SETTINGS & GROUND TRUTH GENERATION
% =========================================================================

% Simulation time settings
epoch = datetime(2026,2,10,0,0,0,'TimeZone','UTC'); % Start date
T_total = 60*45;                                    % Total simulation time [s]
t = 0:MAG.dt:T_total;                               % Time vector [s]
N = numel(t);

fprintf('\n=== Generating Truth Trajectory ===\n');
fprintf(' Simulation time: %.1f s (%d samples)\n', T_total, N);

% --- True Orbital Position (ECI frame) ---
fprintf(' Propagating orbit (J2 + Drag + SRP) ... \n');

    % Spacecraft Parameters
    scParams.mass     = 100.0; % [kg]
    scParams.areaDrag = 1.0;   % [m²]
    scParams.areaSRP  = 1.2;   % [m²]
    scParams.Cd       = 2.0;   %  [-] Drag Coeff
    scParams.Cr       = 1.5;   %  [-]  SRP Coeff (1 = Absorption, 2 = Reflection)
    
    % Orbital Elements
    orbitalElems.SMA   = 6378e3 + 200e3 + 400e3 * rand();
    orbitalElems.ECC   =                  0.01 * rand();
    orbitalElems.INC   =                   180 * rand();
    orbitalElems.RAAN  =                   360 * rand();
    orbitalElems.AOP   =                   360 * rand();
    orbitalElems.TA    =                   360 * rand();
    orbitalElems.epoch = epoch;
    
    [rECI, ~] = simulateOrbit(t, orbitalElems, scParams);

% --- ECI to ECEF and Geodetic Coordinates (LLA) ---
    % Rotate ECI positions to ECEF
    DCM_ECI2ECEF = dcmeci2ecef('IAU-2000/2006', epoch + seconds(t'));
    rECEF = zeros(3,N);
    for k = 1:N
        rECEF(:,k) = DCM_ECI2ECEF(:,:,k) * rECI(:,k);
    end
    
    % Convert ECEF to Geodetic (LLA) using WGS84
    LLA = ecef2lla(rECEF', 'WGS84');
    lat = LLA(:,1)'; % [deg] 1xN
    lon = LLA(:,2)'; % [deg] 1xN
    alt = LLA(:,3)'; %   [m] 1xN

% --- Magnetic Field Truth ---
fprintf(' Using IGRF model ... \n');

    B_ECI = zeros(3,N);
    decimalYr = 2025 * ones(1,N);
    
    XYZ_NED                     = igrfmagm(alt, lat, lon, decimalYr, 13);
    [B_ECEFx, B_ECEFy, B_ECEFz] = ned2ecefv(XYZ_NED(:,1), XYZ_NED(:,2), XYZ_NED(:,3), ...
                                            lat(:), lon(:));
    B_ECEF                      = [B_ECEFx'; B_ECEFy'; B_ECEFz'];
    
    DCM_ECEF2ECI = permute(DCM_ECI2ECEF, [2, 1, 3]);
    for k = 1:N
        B_ECI(:,k) = DCM_ECEF2ECI(:,:,k) * B_ECEF(:,k);
    end

% --- Truth Attitude & Body Frame Field ---
    angleTrue = rand * 2*pi;                    % Rotation angle [rad]
     axisTrue = randn(3,1);
     axisTrue = axisTrue / norm(axisTrue);      % Unit rotation axis
    
    qTrue      = zeros(4,N);
    qTrue(:,1) = [cos(angleTrue/2); ...         % Scalar part
                  sin(angleTrue/2) * axisTrue]; % Vector part
    
    omegaTrue_body      = zeros(3,N); % [rad/s] body rates
    omegaTrue_body(1,:) = deg2rad(0.15);
    omegaTrue_body(2,:) = deg2rad(0.10) * sin(2*pi*(1/300)*t);
    omegaTrue_body(3,:) = deg2rad(0.08) * cos(2*pi*(1/500)*t);
    
    for k = 2:N
        qDot       = 1/2 * skewMatrix(omegaTrue_body(:,k-1)) * qTrue(:,k-1);
        
        qTrue(:,k) = qTrue(:,k-1) + qDot * MAG.dt;
        qTrue(:,k) = qTrue(:,k) / norm(qTrue(:,k));
    end
    
    BTrue_body = zeros(3,N);
    for k = 1:N
        DCM_ECI2B  = quat2dcm(qTrue(:,k));
        BTrue_body(:,k) = DCM_ECI2B * B_ECI(:,k);
    end

fprintf('\n=== MAG Ground Truth Complete (Mean ||B||: %.0f nT) ===\n', mean(vecnorm(BTrue_body)));

%% ========================================================================
% 3. GENERATE SYNTHETIC MAGNETOMETER MEASUREMENTS
% =========================================================================

meas = generateMAGMeasurements(t, BTrue_body, MAG);

%% ========================================================================
% 4. VALIDATION AND ERROR ANALYSIS
% =========================================================================

fprintf('\n=== MAG Validation and Error Analysis ===\n');

% 1. Reference Truth in Sensor Frame (Ideal Alignment)
%    Corresponds to Step 1 in generateMAGMeasurements: B_MAG = DCM_mount * B_Body
BTrue_MAG = MAG.DCM_mounting * BTrue_body;

% 2. RMS Error Analysis (Analog Signal vs. Ideal Truth)
%    Comparison using the unquantized analog signal (BClean) to see intrinsic sensor errors
%    (includes calibration errors, bias, and noise)
BErr  = meas.BClean - BTrue_MAG;
rms_B = sqrt(mean(BErr.^2, 2));

fprintf('\n--- RMS Error vs. Ideal MAG Truth (Analog) ---\n');
fprintf('  RMS(B_x):  %.2f nT\n', rms_B(1));
fprintf('  RMS(B_y):  %.2f nT\n', rms_B(2));
fprintf('  RMS(B_z):  %.2f nT\n', rms_B(3));

% 3. Noise Verification
%    Reconstruct noise by subtracting known deterministic and bias components from the analog signal
%    w_hat = BClean - (BDet + biasTotal)
noiseEst = meas.BClean - (meas.BDet + meas.biasTotal); 
stdNoise = std(noiseEst, 0, 2);

fprintf('\n--- Empirical White Noise Estimation ---\n');
fprintf('  Measured Std:    [%.3f, %.3f, %.3f] nT\n', stdNoise(1), stdNoise(2), stdNoise(3));
fprintf('  Theoretical RMS: %.3f nT\n', MAG.sigmaWhite);

% 4. Angular Alignment Error
%    Compare the direction of the final digital measurement vs. the ideal sensor truth
uTrue = BTrue_MAG ./ vecnorm(BTrue_MAG, 2, 1);
uMeas = meas.B    ./ vecnorm(   meas.B, 2, 1); % Using final quantized output

cosAng = max(min(sum(uTrue .* uMeas, 1), 1), -1); % Clamp to [-1, 1] for numerical safety
angErr = rad2deg(acos(cosAng));

fprintf('\n--- Angular Alignment Error ---\n');
fprintf('  Mean:   %.4f deg\n', mean(angErr));
fprintf('  Std:    %.4f deg\n', std(angErr));
fprintf('  Max:    %.4f deg\n', max(angErr));

%% ========================================================================
% 5. RESULTS VISUALIZATION
% =========================================================================

% Generate comprehensive validation plots
saveFlag = 1;
plotMAGResults(t, BTrue_body, meas, MAG, saveFlag);

fprintf('\n=== All tasks completed successfully ===\n');