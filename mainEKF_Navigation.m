%==========================================================================
% mainEKF_Navigation.m - Position & Velocity Extended Kalman Filter (PV-EKF)
%==========================================================================
%
% DESCRIPTION:
%   Complete simulation and validation framework for spacecraft translational
%   navigation using an Extended Kalman Filter (PV-EKF). This script fuses 
%   high-rate inertial data (IMU/Accelerometer) with low-rate Position, 
%   Velocity, and Time (PVT) observations from a GNSS receiver to estimate 
%   the spacecraft's true orbital state and dynamic accelerometer bias.
%
% PURPOSE:
%   Demonstrates state-of-the-art translational sensor fusion for S/C GNC.
%   Validates orbital propagation inside the filter, covariance bounding, 
%   and the filter's ability to reject GNSS correlated noise (Gauss-Markov)
%   and lever-arm kinematic disturbances.
%
% ALGORITHM:
%   The PV-EKF estimates a 9-state vector: x = [r_ECI; v_ECI; b_accel]
%     1. Prediction: Integrates orbital dynamics (central gravity + J2) and
%        adds the rotated IMU specific force to propagate the state. 
%        Propagates the covariance matrix (P) using the state transition 
%        matrix (Φ) and tuned process noise (Q).
%     2. Update: When GNSS measurements arrive (typically 1 Hz), computes 
%        the Kalman Gain (K) to update position and velocity, while 
%        simultaneously estimating the accelerometer bias.
%
% WORKFLOW:
%   1. Initialize hardware parameters (IMU and GNSS)
%   2. Initialize spacecraft and orbital parameters
%   3. Generate centralized ground truth (Orbit, Attitude, Environment)
%   4. Simulate multi-rate noisy sensor measurements (including Lever Arm)
%   5. Initialize PV-EKF (Uses first GNSS fix for initial state)
%   6. Execute EKF loop (Predict @ 120Hz, Update GNSS @ 1Hz)
%   7. Validate results: Position/Velocity errors, bias drift, NEES
%   8. Generate diagnostic plots
%
% CONFIGURATION:
%   - Sensor parameters:      Edit initializeIMU.m, initializeGNSS.m
%   - Orbital elements:       Edit initializeOrbit.m
%   - Spacecraft & Attitude:  Edit initializeSpacecraft.m
%   - Filter Tuning (Q/R):    Edit initializeEKF_Nav.m
%==========================================================================

close all;
clear;
clc;

%% ========================================================================
% 1. HARDWARE PARAMETERS & SENSOR INITIALIZATION
% =========================================================================
fprintf('\n=== Navigation EKF Initialization ===\n');

% IMU (Gyroscope + Accelerometer)
IMU = initializeIMU();

% GNSS Receiver
GNSS = initializeGNSS();

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

% Specific Force (Accelerometer Truth)
% In free-fall orbit, gravity is NOT measured by the IMU. 
% For a passive LEO spacecraft, specific force is strictly drag & SRP (~1e-6 m/s²).
% We set it to zero to cleanly simulate a free-flying unpowered spacecraft.
forceTrue_body = zeros(3,N); 

fprintf('=== Truth Generation Complete ===\n');

%% ========================================================================
% 3. SENSOR MEASUREMENTS GENERATION
% =========================================================================
fprintf('\n=== Generating Sensor Measurements ===\n');

% IMU Measurements (Gyro + Accel at 120 Hz)
imuMeas = generateIMUMeasurements(t, omegaTrue, forceTrue_body, IMU);

% GNSS Measurements (Position & Velocity at 1 Hz)
gnssMeas = generateGNSSMeasurements(t, rECI, vECI, qTrue, omegaTrue, GNSS);

fprintf('=== Sensor Data Generation Complete ===\n');

%% ========================================================================
% 4. INITIALIZE PV-EKF
% =========================================================================
EKF = initializeEKF_Nav(IMU, GNSS);

% Initialize filter states with the very first GNSS measurement
% Assuming first measurement is available at t=0
EKF.x(1:3) = gnssMeas.rECI(:,1);
EKF.x(4:6) = gnssMeas.vECI(:,1);
EKF.x(7:9) = zeros(3,1); % Initial guess for accel bias is 0

fprintf('\n Initialized PV-EKF state with first GNSS reading.\n');

%% ========================================================================
% 5. EKF PROPAGATION AND UPDATE LOOP
% =========================================================================
fprintf('\n=== Running PV-EKF ===\n');

% Pre-allocate history
rEst         = zeros(3,N);
vEst         = zeros(3,N);
biasEst      = zeros(3,N);
PHist        = zeros(9,9,N);

rEst(:,1)    = EKF.x(1:3);
vEst(:,1)    = EKF.x(4:6);
biasEst(:,1) = EKF.x(7:9);
PHist(:,:,1) = EKF.P;

% GNSS update step interval (e.g., 120 samples if IMU=120Hz and GNSS=1Hz)
gnssStep = round(IMU.rate * GNSS.dt);

for k = 2:N
    % 1) Prediction Step (Integrate Orbit + IMU Accel at high rate)
    % Note: In a fully coupled S/C, qTrue would be replaced by qEst from the Attitude EKF.
    % Here we use qTrue to cleanly isolate the Navigation filter's performance.
    EKF = predictEKF_Nav(EKF, imuMeas.accel.forceBody(:,k), qTrue(:,k), IMU.dt);
    
    % 2) Correction Step (GNSS update at low rate)
    if mod(k - 1, gnssStep) == 0
        % Update using GNSS measurements. We pass true attitude/rates to 
        % allow the filter to compensate for the Lever Arm effect internally.
        EKF = updateGNSS(EKF, gnssMeas.rECI(:,k), gnssMeas.vECI(:,k), ...
                         GNSS.leverArm, qTrue(:,k), omegaTrue(:,k));
    end
    
    % 3) Store history
    rEst(:,k)    = EKF.x(1:3);
    vEst(:,k)    = EKF.x(4:6);
    biasEst(:,k) = EKF.x(7:9);
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
plotEKF_NavResults(t, groundTruth, rEst, vEst, biasEst, PHist, saveFlag);

fprintf('\n=== All tasks completed successfully ===\n');