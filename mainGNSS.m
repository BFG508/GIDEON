%==========================================================================
% mainGNSS.m - Global Navigation Satellite System (GNSS) Simulation
%==========================================================================
%
% DESCRIPTION:
%   Simulation framework for generating high-fidelity synthetic GNSS PVT 
%   (Position, Velocity, Time) data for spacecraft Navigation filters. 
%   This script models realistic error sources typically found in space-grade 
%   GNSS receivers, including the kinematic Lever Arm effect, wideband white 
%   noise, and temporally correlated errors (Gauss-Markov) caused by 
%   ionospheric delays, ephemeris errors, and clock wander.
%
% PURPOSE:
%   Demonstrates high-fidelity translational sensor modeling. Validates the
%   combination of high-frequency and low-frequency noise characteristics,
%   and quantifies the impact of attitude dynamics on antenna position 
%   (Lever Arm coupling). Provides ground truth for PV-EKF algorithm testing.
%
% ALGORITHM:
%   The measurement model simulates the PVT output:
%     r_meas = r_CG + DCM_B2I * L_arm + b_pos_GM + noise_pos_white
%     v_meas = v_CG + DCM_B2I * (omega x L_arm) + b_vel_GM + noise_vel_white
%
% WORKFLOW:
%   1. Define GNSS hardware parameters via initializeGNSS.m
%   2. Initialize spacecraft and orbital parameters via helper functions
%   3. Generate centralized ground truth (Orbit & Attitude)
%   4. Simulate GNSS measurements with Lever Arm and realistic error corruption
%   5. Validate error models: RMS analysis, Gauss-Markov verification, and 
%      Lever Arm displacement analysis.
%   6. Generate diagnostic plots for filter tuning
%
% CONFIGURATION:
%   - Sensor parameters:      Edit initializeGNSS.m
%   - Orbital elements:       Edit initializeOrbit.m
%   - Spacecraft & Attitude:  Edit initializeSpacecraft.m
%==========================================================================

clear;
close all;
clc;

%% ========================================================================
% 1. GNSS HARDWARE PARAMETERS (Space-Grade Receiver specs)
% =========================================================================

GNSS = initializeGNSS();

%% ========================================================================
% 2. SIMULATION SETTINGS & GROUND TRUTH GENERATION
% =========================================================================

% Simulation time settings
epoch   = datetime(2026,1,1,0,0,0,'TimeZone','UTC'); % Start date
T_total = 60 * 60;                                   % Total simulation time [s]
t       = 0:GNSS.dt:T_total;                         % Time vector [s]
N       = numel(t);

% Initialize orbit and S/C models
orbitalElems          = initializeOrbit(epoch);
[scParams, attParams] = initializeSpacecraft();

% Generate centralized ground truth
truth = generateGroundTruth(t, epoch, orbitalElems, scParams, attParams);

% Extract true translation and attitude variables for GNSS coupling
rECI_true  = truth.rECI;
vECI_true  = truth.vECI;
qTrue      = truth.qTrue;
omegaTrue  = truth.omegaTrue;

fprintf('\n=== GNSS Ground Truth Complete ===\n');

%% ========================================================================
% 3. GENERATE SYNTHETIC GNSS MEASUREMENTS
% =========================================================================

meas = generateGNSSMeasurements(t, rECI_true, vECI_true, qTrue, omegaTrue, GNSS);

%% ========================================================================
% 4. VALIDATION AND ERROR ANALYSIS
% =========================================================================

fprintf('\n=== GNSS Validation and Error Analysis ===\n');

% --- 4.1 Position Error Analysis ---
% Calculate total error compared to the TRUE ANTENNA position (Clean)
posErr = meas.rECI - meas.rClean;
rmsPos = sqrt(mean(posErr.^2, 2));

% Expected theoretical total standard deviation per axis: sqrt(sigma_W^2 + sigma_GM^2)
expectedPosRMS = sqrt(GNSS.sigmaPosWhite^2 + GNSS.sigmaPosGM^2);

fprintf('\n--- Position Measurement Error (vs Antenna Truth) ---\n');
fprintf('  RMS(X): %.2f m\n', rmsPos(1));
fprintf('  RMS(Y): %.2f m\n', rmsPos(2));
fprintf('  RMS(Z): %.2f m\n', rmsPos(3));
fprintf('  Theoretical Expected 1-sigma: %.2f m\n', expectedPosRMS);

% --- 4.2 Velocity Error Analysis ---
velErr = meas.vECI - meas.vClean;
rmsVel = sqrt(mean(velErr.^2, 2));

expectedVelRMS = sqrt(GNSS.sigmaVelWhite^2 + GNSS.sigmaVelGM^2);

fprintf('\n--- Velocity Measurement Error (vs Antenna Truth) ---\n');
fprintf('  RMS(Vx): %.4f m/s\n', rmsVel(1));
fprintf('  RMS(Vy): %.4f m/s\n', rmsVel(2));
fprintf('  RMS(Vz): %.4f m/s\n', rmsVel(3));
fprintf('  Theoretical Expected 1-sigma: %.4f m/s\n', expectedVelRMS);

% --- 4.3 Gauss-Markov Correlated Bias Verification ---
stdPosGM = std(meas.posBiasDyn, 0, 2);
stdVelGM = std(meas.velBiasDyn, 0, 2);

fprintf('\n--- Gauss-Markov Correlated Bias (Ephemeris/Iono/Clock) ---\n');
fprintf('  Position GM Std: [%.2f, %.2f, %.2f] m    (Target: %.2f m)\n', ...
        stdPosGM(1), stdPosGM(2), stdPosGM(3), GNSS.sigmaPosGM);
fprintf('  Velocity GM Std: [%.4f, %.4f, %.4f] m/s (Target: %.4f m/s)\n', ...
        stdVelGM(1), stdVelGM(2), stdVelGM(3), GNSS.sigmaVelGM);

% --- 4.4 Lever Arm Effect Analysis ---
% Difference between S/C Center of Mass and actual Antenna phase center
leverArmDisplacement = meas.rClean - rECI_true;
maxDisplacement      = max(vecnorm(leverArmDisplacement));
meanLeverVel         = mean(vecnorm(meas.vClean - vECI_true));

fprintf('\n--- Kinematic Lever Arm Effect ---\n');
fprintf('  Static lever arm (body):      [%.2f, %.2f, %.2f] m\n', GNSS.leverArm);
fprintf('  Max dynamic 3D offset:        %.3f m (Theoretical = %.3f m)\n', ...
        maxDisplacement, norm(GNSS.leverArm));
fprintf('  Mean velocity disturbance:    %.4f m/s induced by attitude rates\n', meanLeverVel);

% --- 4.5 Satellite Visibility ---
fprintf('\n--- Constellation Tracking ---\n');
fprintf('  Satellites tracked: Mean = %.1f, Min = %d, Max = %d\n', ...
        mean(meas.nSats), min(meas.nSats), max(meas.nSats));

%% ========================================================================
% 5. RESULTS VISUALIZATION
% =========================================================================

% Generate comprehensive validation plots
saveFlag = 1;
plotGNSSResults(t, truth, meas, GNSS, saveFlag);

fprintf('\n=== All tasks completed successfully ===\n');