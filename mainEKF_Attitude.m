%==========================================================================
% mainEKF_Attitude.m - Multiplicative Extended Kalman Filter (MEKF) for Attitude
%==========================================================================
%
% DESCRIPTION:
%   Complete simulation and validation framework for spacecraft attitude
%   estimation using a Multiplicative Extended Kalman Filter (MEKF). This 
%   script fuses high-rate inertial data (IMU/Gyroscope) with low-rate 
%   vector observations from Star Trackers (STR) and Magnetometers (MAG)
%   to estimate the true orientation and dynamic gyro bias.
%
% PURPOSE:
%   Demonstrates state-of-the-art sensor fusion for spacecraft AOCS/GNC.
%   Validates filter convergence, covariance bounding (3-sigma consistency),
%   and robustness against realistic sensor noise (ARW, RRW, scale factors) 
%   and multi-rate measurement architectures.
%
% ALGORITHM:
%   The MEKF estimates a 6-state error vector: x = [δθ; Δβ]
%     1. Prediction: Integrates gyro kinematics using the Hamilton passive 
%        quaternion convention. Propagates the covariance matrix (P) using
%        the state transition matrix (Φ) and tuned process noise (Q).
%     2. Update: When STR or MAG measurements arrive, computes the Kalman 
%        Gain (K) using linearized observation matrices (H) and measurement 
%        noise (R). Updates the error state.
%     3. Reset: Multiplicatively injects the angle error (δθ) into the 
%        nominal quaternion, additively updates the bias, and resets the 
%        error state to zero.
%
% WORKFLOW:
%   1. Initialize hardware parameters (IMU, MAG, 2x STR, Star Catalog)
%   2. Generate ground truth trajectory (LEO orbit, Nadir-pointing attitude)
%   3. Simulate multi-rate noisy sensor measurements
%   4. Initialize MEKF (Uses TRIAD algorithm for the first frame)
%   5. Execute EKF loop (Predict @ 120Hz, Update MAG @ 10Hz, STR @ 4Hz)
%   6. Validate results: Attitude error, bias drift, NEES, covariance bounds
%   7. Generate diagnostic plots (publication-ready figures)
%
% KEY FEATURES:
%   - Multi-rate asynchronous sensor fusion (120Hz / 10Hz / 4Hz)
%   - Space-grade error models with split covariance tuning
%   - Realistic LEO orbital dynamics and magnetic field (IGRF) propagation
%   - Hamilton passive quaternion algebra
%   - TRIAD-based autonomous initialization
%   - Comprehensive statistical validation (Average NEES, RMS bounds)
%
% OUTPUTS:
%   - qEst:               Estimated attitude quaternion history [4xN]
%   - biasEst:            Estimated gyroscope bias history [3xN]
%   - PHist:              Covariance matrix history [6x6xN]
%   - Diagnostic figures: Attitude error, bias error, uncertainty evolution, 
%                         NEES consistency, and angular rates.
%
% CONFIGURATION:
%   Edit Section 1/2 for: Orbital elements, simulation duration, pointing mode.
%   Edit initializeEKF_Att.m for: Tuning parameters (Q and R matrices).
%==========================================================================

close all;
clear;
clc;

%% ========================================================================
% 1. HARDWARE PARAMETERS & SENSOR INITIALIZATION
% =========================================================================

fprintf('\n=== EKF Initialization ===\n');

% IMU (Gyroscope + Accelerometer)
IMU = initializeIMU();

% Magnetometer
MAG = initializeMAG();

% Star Tracker (2 units)
nSTR = 2;
STR  = initializeSTR(nSTR);

% Load stellar catalog for STR simulation
catalog = loadHipparcosStellarCatalog(STR(1).magnitudeLimit);

%% ========================================================================
% 2. TRUTH TRAJECTORY GENERATION
% =========================================================================
fprintf('\n=== Truth Trajectory Generation ===\n');

% Simulation settings
epoch   = datetime(2026,1,1,0,0,0,'TimeZone','UTC'); % Start date
T_total = 60 * 60;                                   % Total simulation time   [s]
t       = 0:IMU.dt:T_total;                          % Time vector at IMU rate [s]
N       = numel(t);

fprintf(' Simulation duration: %.1f min (%d samples at %.0f Hz)\n', ...
        T_total/60, N, IMU.rate);

% Initialize Orbit and Spacecraft Models
orbitalElems          = initializeOrbit(epoch);
[scParams, attParams] = initializeSpacecraft();

% Generate Centralized Ground Truth
groundTruth = generateGroundTruth(t, epoch, orbitalElems, scParams, attParams);

% Extract truth variables for convenience
qTrue       = groundTruth.qTrue;     % Attitude quaternion (ECI to Body),          4xN
omegaTrue   = groundTruth.omegaTrue; % Angular velocity (Body frame)      [rad/s], 3xN
rECI        = groundTruth.rECI;      % Position in ECI                    [m],     3xN
vECI        = groundTruth.vECI;      % Velocity in ECI                    [m/s],   3xN
B_ECI       = groundTruth.B_ECI;     % Magnetic field in ECI              [nT],    3xN

% Rotate magnetic field to body frame for MAG
BTrue_body = zeros(3,N);
for k = 1:N
    DCM_ECI2B       = quat2dcm(qTrue(:,k));
    BTrue_body(:,k) = DCM_ECI2B * B_ECI(:,k);
end

% Compute specific force for accelerometer (1g assumption for static-like tuning)
g0             = 9.80665;                   % [m/s²]
forceTrue_body = repmat([0; 0; -g0], 1, N); % [m/s²]

fprintf('=== Truth Generation Complete ===\n');

%% ========================================================================
% 3. SENSOR MEASUREMENTS GENERATION
% =========================================================================

fprintf('\n=== Generating Sensor Measurements ===\n');

% IMU Measurements (Gyro + Accel at 120 Hz)
imuMeas = generateIMUMeasurements(t, omegaTrue, forceTrue_body, IMU);

% Magnetometer Measurements (at 10 Hz)
magMeas = generateMAGMeasurements(t, BTrue_body, MAG);

% Star Tracker Measurements (at 4 Hz, only when stars are visible)
strMeas = cell(1,N);
for k = 1:N
    if mod(k-1, round(IMU.rate/STR(1).rate)) == 0
        DCMTrue_k   = quat2dcm(qTrue(:,k));
        omegaBody_k = omegaTrue(:,k);
        strMeas{k}  = generateSTRMeasurements(STR, nSTR, catalog, DCMTrue_k, omegaBody_k);
    else
        strMeas{k} = []; % No measurement at this time step
    end
end

fprintf('=== Sensor Data Generation Complete ===\n');

%% ========================================================================
% 4. INITIALIZE MEKF
% =========================================================================

EKF = initializeEKF_Att(IMU, MAG, STR);

% Initialize with first TRIAD solution (if STR data available)
if ~isempty(strMeas{1})
    fprintf('\n Initializing attitude with TRIAD solution ... ');
    [qInit, ~, ~] = solveTRIADAttitude(strMeas{1}, nSTR);
    EKF.qNom      = qInit;
    fprintf('Done.\n');
else
    warning('No STR data at t = 0, using identity quaternion.');
end

%% ========================================================================
% 5. EKF PROPAGATION AND UPDATE LOOP
% =========================================================================

fprintf('\n=== Running EKF ===\n');

% Pre-allocate history
qEst         = zeros(4,N);
biasEst      = zeros(3,N);
PHist        = zeros(6,6,N);

qEst(:,1)    = EKF.qNom;
biasEst(:,1) = EKF.x(4:6);
PHist(:,:,1) = EKF.P;

for k = 2:N
    % 1) Prediction Step (Gyro propagation at IMU rate)
    EKF = predictEKF_Att(EKF, imuMeas.gyro.omegaBody(:,k), IMU.dt);
    
    % 2) Correction Step (Sensor updates at their respective rates)
        % Magnetometer update (10 Hz)
        if mod(k - 1, round(IMU.rate/MAG.rate)) == 0
            EKF = updateMAG(EKF, magMeas.B(:,k), B_ECI(:,k));
        end

        % Star Tracker update (4 Hz, when available)
        if ~isempty(strMeas{k})
            EKF = updateSTR(EKF, strMeas{k}, nSTR);
        end
    
    % 3) Store history
    qEst(:,k)    = EKF.qNom;
    biasEst(:,k) = EKF.x(4:6);
    PHist(:,:,k) = EKF.P;
    
    % 4) Progress indicator
    if mod(k, round(N/10)) == 0
        fprintf(' Progress: %.0f%%\n', 100*k/N);
    end
end

fprintf('=== EKF Propagation Complete ===\n');

%% ========================================================================
% 6. ERROR ANALYSIS AND VISUALIZATION
% =========================================================================

% Generate comprehensive validation plots
saveFlag = 1;
plotEKF_AttResults(t, groundTruth, qEst, biasEst, PHist, imuMeas, saveFlag);

fprintf('\n=== All tasks completed successfully ===\n');