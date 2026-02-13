function EKF = initializeEKF_Att(IMU, MAG, STR)
%==========================================================================
% initializeMEKF - Initialize Multiplicative Extended Kalman Filter (MEKF)
%                  for spacecraft attitude estimation
%
% Inputs:
%   IMU - IMU parameter structure       (from initializeIMU.m)
%   MAG - MAG parameter structure       (from initializeMAG.m)
%   STR - STR parameter array structure (from initializeSTR.m)
%
% Outputs:
%   EKF - EKF state structure with fields:
%         .x      - Error state vector [6x1]: [δθ; biasGyro]
%         .P      - Error covariance matrix [6x6]
%         .qNom   - Nominal attitude quaternion [4x1] (scalar-first)
%         .Q      - Process noise covariance [6x6]
%         .R_MAG  - Magnetometer measurement noise [3x3]
%         .R_STR  - Star Tracker measurement noise [3x3]
%
% MEKF State Definition:
%   x = [δθ (3x1); biasGyro (3x1)]
%   where:
%     δθ      : Attitude error (small angle, Gibbs vector or MRP)
%     biasGyro: Gyroscope bias error [rad/s]
%
%   The true attitude is maintained separately as a quaternion:
%     q_true = δq(δθ) ⊗ qNom
%==========================================================================

    fprintf('\n=== Initializing MEKF ===\n');
    
    %% ===================================================================
    % 1. ERROR STATE INITIALIZATION
    % ====================================================================

    % Initial error state (zero error assumption)
    EKF.x = zeros(6,1);  % [δθ; biasGyro]
    
    %% ===================================================================
    % 2. INITIAL COVARIANCE MATRIX
    % ====================================================================

    % Attitude uncertainty: assume ±10° initial pointing error (1σ)
    sigmaAtt0 = deg2rad(10); % [rad]
    PAtt      = (sigmaAtt0^2) * eye(3);
    
    % Gyro bias uncertainty: use turn-on bias repeatability spec
    sigmaBias0 = deg2rad(IMU.gyro.biasStaticLim / 3600);  % [rad/s]
    PBias      = (sigmaBias0^2) * eye(3);
    
    EKF.P = blkdiag(PAtt, PBias);
    
    fprintf(' Initial covariance:\n');
    fprintf('   - Attitude uncertainty:  %.2f deg (1σ)\n', ...
            rad2deg(sqrt(PAtt(1,1))));
    fprintf('   - Gyro bias uncertainty: %.4f deg/h (1σ)\n', ...
            rad2deg(sqrt(PBias(1,1))) * 3600);
    
    %% ===================================================================
    % 3. NOMINAL ATTITUDE
    % ====================================================================

    EKF.qNom = [1; 0; 0; 0]; % Identity quaternion (ECI = Body initially)
    
    %% ===================================================================
    % 4. PROCESS NOISE COVARIANCE (from Allan Deviation analysis)
    % ====================================================================

    EKF.Q = computeIMUNoise(IMU) * 10^2;
    
    fprintf(' Process noise (Q):\n');
    fprintf('   - ARW noise (gyro):       %.2e (rad/s)²/s\n', EKF.Q(1,1));
    fprintf('   - RRW noise (bias drift): %.2e (rad/s)²/s\n', EKF.Q(4,4));
    
    %% ===================================================================
    % 5. MEASUREMENT NOISE COVARIANCES
    % ====================================================================

    % Magnetometer noise
    EKF.R_MAG = computeMAGNoise(MAG) + (MAG.hardIronLim)^2 * eye(3);
    fprintf(' Magnetometer noise (R):  %.2f nT (1σ per axis)\n', ...
            sqrt(EKF.R_MAG(1,1)));
    
    % Star Tracker noise
    EKF.R_STR = computeSTRNoise(STR);
    fprintf(' Star Tracker noise (R):  %.2f arcsec (1σ cross-boresight)\n', ...
            rad2deg(sqrt(EKF.R_STR(1,1))) * 3600);
    
    fprintf('=== MEKF Initialization Complete ===\n');

end