%==========================================================================
% mainEKF.m - Multiplicative Extended Kalman Filter (MEKF) for Attitude
%==========================================================================
clear; close all; clc;

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
epoch   = datetime(2026,2,13,0,0,0,'TimeZone','UTC'); % Start date
T_total = 60 * 45;                                    % Total simulation time [s]
t       = 0:IMU.dt:T_total;                           % Time vector at IMU rate [s]
N       = numel(t);

fprintf(' Simulation duration: %.1f min (%d samples at %.0f Hz)\n', ...
        T_total/60, N, IMU.rate);

% Spacecraft Parameters
scParams.mass     = 100.0; % [kg]
scParams.areaDrag = 1.0;   % [m²]
scParams.areaSRP  = 1.2;   % [m²]
scParams.Cd       = 2.0;   %  [-] Drag Coeff
scParams.Cr       = 1.5;   %  [-]  SRP Coeff (1 = Absorption, 2 = Reflection)

% Orbital elements (randomized LEO)
orbitalElems.SMA   = 6378e3 + 200e3 + 400e3 * rand(); % [m]
orbitalElems.ECC   =                   0.01 * rand(); % [-]
orbitalElems.INC   =                    180 * rand(); % [deg]
orbitalElems.RAAN  =                    360 * rand(); % [deg]
orbitalElems.AOP   =                    360 * rand(); % [deg]
orbitalElems.TA    =                    360 * rand(); % [deg]
orbitalElems.epoch = epoch;

% Attitude parameters (Nadir Pointing with realistic perturbations)
attParams.I      = diag(5 + 45*[rand(), ...
                                rand(), ...
                                rand()]);                % Inertia tensor [kg·m²]
attParams.dipole =  [-0.1 +  0.2*rand(); 
                     -0.1 +  0.2*rand(); 
                     0.05 + 0.25*rand()];                % Magnetic dipole [A·m²]
attParams.Kp     =   0.1 +   1.9*rand();                 % Proportional gain [N·m]
attParams.Ki     = 0.001 + 0.099*rand();                 % Integral gain [N·m/s]
attParams.Kd     =   0.5 +   4.5*rand();                 % Derivative gain [N·m·s]
attParams.mode   = 'LVLH';                               % Nadir pointing mode

% Generate complete truth trajectory
groundTruth = generateGroundTruth(t, epoch, orbitalElems, scParams, attParams);

% Extract truth variables for convenience
qTrue       = groundTruth.qTrue;     % Attitude quaternion (ECI to Body),          4xN
omegaTrue   = groundTruth.omegaTrue; % Angular velocity (body frame)      [rad/s], 3xN
rECI        = groundTruth.rECI;      % Position in ECI                    [m],     3xN
vECI        = groundTruth.vECI;      % Velocity in ECI                    [m/s],   3xN
B_ECI       = groundTruth.B_ECI;     % Magnetic field in ECI              [nT],    3xN

% Rotate magnetic field to body frame for MAG
BTrue_body = zeros(3,N);
for k = 1:N
    DCM_ECI2B       = quat2dcm(qTrue(:,k));
    BTrue_body(:,k) = DCM_ECI2B * B_ECI(:,k);
end

% Compute specific force for accelerometer
g0             = 9.80665;                   % [m/s²]
forceTrue_body = repmat([0; 0; -g0], 1, N); % [m/s²]

fprintf('=== Truth Generation Complete ===\n');

%% ========================================================================
% 3. SENSOR MEASUREMENTS GENERATION
% =========================================================================

fprintf('\n=== Generating Sensor Measurements ===\n');

% IMU Measurements (Gyro + Accel at 120 Hz)
fprintf(' Generating IMU measurements ... ');
imuMeas = generateIMUMeasurements(t, omegaTrue, forceTrue_body, IMU);
fprintf('Done.\n');

% Magnetometer Measurements (at 10 Hz)
fprintf(' Generating MAG measurements ... ');
magMeas = generateMAGMeasurements(t, BTrue_body, MAG);
fprintf('Done.\n');

% Star Tracker Measurements (at 4 Hz, only when stars are visible)
fprintf(' Generating STR measurements ... ');
strMeas = cell(1,N);
for k = 1:N
    % STR operates at its own rate
    if mod(k-1, round(IMU.rate/STR(1).rate)) == 0
        DCMTrue_k   = quat2dcm(qTrue(:,k));
        omegaBody_k = omegaTrue(:,k);
        strMeas{k}  = generateSTRMeasurements(STR, nSTR, catalog, DCMTrue_k, omegaBody_k);
    else
        strMeas{k} = []; % No measurement at this time step
    end
end
fprintf('Done.\n');

fprintf('=== Sensor Data Generation Complete ===\n');

%% ========================================================================
% 4. INITIALIZE MEKF
% =========================================================================

EKF = initializeEKF_Att(IMU, MAG, STR);

% Initialize with first TRIAD solution (if STR data available)
if ~isempty(strMeas{1})
    fprintf('\n Initializing attitude with QUEST solution ... ');
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