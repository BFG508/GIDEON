%==========================================================================
% mainMAG.m - Magnetometer Simulation
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
%==========================================================================

clear;
close all;
clc;

%% ========================================================================
% 1. MAG HARDWARE PARAMETERS (typical values from Space-Grade Fluxgate/AMR datasheets)
% =========================================================================

% --- General Configuration ---
MAG.rate = 10;           % Sampling frequency [Hz]
MAG.dt   = 1 / MAG.rate; % Sampling time step [s]

% Dynamic Range & Resolution       
MAG.range      = 100000; % [nT] +/- full scale
MAG.resolution = 0.5;    % [nT] ADC quantization step

% Mounting misalignment errors (typical mechanical tolerances)
rollErr  =  0.002; % Mounting error around X-axis [deg]
pitchErr = -0.005; % Mounting error around Y-axis [deg]
yawErr   =  0.003; % Mounting error around Z-axis [deg]

% Construct misalignment DCM (small perturbation from ideal alignment)
% This matrix is shared between gyro and accel (common mechanical mounting).
MAG.DCM_mounting = angle2dcm(deg2rad(yawErr), ...
                             deg2rad(pitchErr), ...
                             deg2rad(rollErr), ...
                             'ZYX');

% --- Stochastic Noise Parameters ---
    MAG.noiseDensity = 0.1;               % White noise spectral density [nT/sqrt(Hz)]
    
    % Anti-alias / measurement bandwidth model for noise integration
    MAG.antiAlias.type = "onePole";       % "idealNyquist" | "onePole"
    MAG.antiAlias.fc   = 0.25 * MAG.rate; % [Hz]
    
    % Bias drift as bounded 1st-order Gauss-Markov process
    MAG.biasInstability = 5.0;            % Steady-state    [nT] (1σ, per axis)
    MAG.biasTau         = 3600;           % Correlation time [s]

% --- Deterministic Errors ---
    MAG.hardIronLim      = 500.0;         % Static hard-iron magnitude limit [nT] (per axis, uniform)
    MAG.SF               = 1000;          % Diagonal scale errors [ppm] (1σ)
    MAG.softIronCoupling = 500;           % Off-diagonal symmetric coupling [ppm] (1σ)
    MAG.nonortho         = 0.1;           % Small non-orthogonality angles [deg] (1σ)

% --- Derived Matrices & Pre-Allocation ---
fprintf('\n=== Initializing MAG Model ===\n');

    % Soft iron / Scale factor (symmetric matrix near identity)
    SFErr  =               (MAG.SF * 1e-6) .* randn(3,1);
    SICErr = (MAG.softIronCoupling * 1e-6) .* randn(3,1); % xy, xz, yz couplings
    MAG.M_SoftIron = eye(3) + [ SFErr(1), SICErr(1), SICErr(2);
                               SICErr(1),  SFErr(2), SICErr(3);
                               SICErr(2), SICErr(3),  SFErr(3)];
     
    fprintf(' Soft-Iron / Scale-Factor Matrix: [%+.6f  %+.6f  %+.6f]\n', MAG.M_SoftIron(1,1), MAG.M_SoftIron(1,2), MAG.M_SoftIron(1,3));
    fprintf('                                  [%+.6f  %+.6f  %+.6f]\n', MAG.M_SoftIron(2,1), MAG.M_SoftIron(2,2), MAG.M_SoftIron(2,3));
    fprintf('                                  [%+.6f  %+.6f  %+.6f]\n', MAG.M_SoftIron(3,1), MAG.M_SoftIron(3,2), MAG.M_SoftIron(3,3));
    
    % Non-orthogonality (small-angle upper triangular correction)
    nonorthErr = deg2rad(MAG.nonortho) .* randn(3,1);
    MAG.M_Nonorth = [             1, -nonorthErr(3),  nonorthErr(2);
                     -nonorthErr(3),              1, -nonorthErr(1);
                      nonorthErr(2), -nonorthErr(1),              1];
    
    % Deterministic MAG matrix acting in MAG frame (after mounting)
    MAG.M_Deterministic = MAG.M_SoftIron * MAG.M_Nonorth * MAG.DCM_mounting;

% Static hard-iron bias (MAG axes)
MAG.biasHardIron = MAG.hardIronLim * (2*rand(3,1) - 1);

fprintf(' Hard Iron Bias (MAG):            [%.2f, %.2f, %.2f] nT\n', MAG.biasHardIron(1), MAG.biasHardIron(2), MAG.biasHardIron(3));

% White noise RMS per sample from density and assumed bandwidth
if string(MAG.antiAlias.type) == "idealNyquist"
    BW = MAG.antiAlias.fc/2;
elseif string(MAG.antiAlias.type) == "onePole"
    BW = 1.57 * MAG.antiAlias.fc; % ENBW for 1-pole low-pass
else
    BW = MAG.antiAlias.fc/2;
end
MAG.sigmaWhite = MAG.noiseDensity * sqrt(BW);

fprintf(' White Noise RMS:                 %.4f nT (per axis)\n', MAG.sigmaWhite);

fprintf('=== MAG Initialization Complete ===\n');

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