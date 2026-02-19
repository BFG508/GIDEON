function EKF = initializeEKF_Nav(IMU, GNSS)
%==========================================================================
% initializeEKF_Nav - Initialize Position & Velocity Extended Kalman Filter
%                     (PV-EKF) for spacecraft translational navigation.
%
% Inputs:
%   IMU  - IMU parameter structure       (from initializeIMU.m)
%   GNSS - GNSS parameter structure      (from initializeGNSS.m)
%
% Outputs:
%   EKF - EKF state structure with fields:
%         .x      - Full state vector [9x1]: [rECI; vECI; biasAccel]
%         .P      - Error covariance matrix [9x9]
%         .Q      - Discrete process noise covariance [9x9]
%         .R_GNSS - GNSS measurement noise covariance [6x6]
%
% PV-EKF State Definition:
%   x = [r (3x1); v (3x1); bA (3x1)]
%   where:
%     r : Position vector in ECI frame [m]
%     v : Velocity vector in ECI frame [m/s]
%     bA: Accelerometer dynamic bias in body frame [m/s²]
%==========================================================================

    fprintf('\n=== Initializing PV-EKF (Navigation) ===\n');

    % Extract sampling time
    dt = IMU.dt;
    
    %% ===================================================================
    % 1. INITIAL FULL STATE VECTOR (x)
    % ====================================================================
    % The filter states will be properly initialized in the main loop using 
    % the first valid GNSS measurement. Here we pre-allocate the 9x1 vector.
    EKF.x = zeros(9, 1);
    
    %% ===================================================================
    % 2. INITIAL ERROR COVARIANCE MATRIX (P0)
    % ====================================================================
    % Defines the initial uncertainty of our state estimate.
    % - Position: ~100 m initial uncertainty
    % - Velocity: ~10 m/s initial uncertainty
    % - Accel Bias: ~0.05 m/s² initial uncertainty
    
    sigmaR0 = 100.0; % [m]
    sigmaV0 = 10.0;  % [m/s]
    sigmaB0 = 0.05;  % [m/s²]
    
    P0Pos  = (sigmaR0^2) * eye(3);
    P0Vel  = (sigmaV0^2) * eye(3);
    P0Bias = (sigmaB0^2) * eye(3);
    
    EKF.P = blkdiag(P0Pos, P0Vel, P0Bias);
    
    %% ===================================================================
    % 3. PROCESS NOISE COVARIANCE (Q)
    % ====================================================================
    % Computes the continuous-time Power Spectral Density (PSD) of the 
    % process noise and discretizes it for the EKF propagation step.
    
    % --- 3.1 Base Noise from IMU Datasheet ---
    [qVRW, qARW] = computeIMUNoise_Nav(IMU);
    
    % --- 3.2 Process Noise Tuning Factors ---
    % Since the S/C operates in LEO, the orbital dynamics (gravity, drag) 
    % are well modeled, but density variations (NRLMSISE-00) and attitude 
    % coupling inject unmodeled perturbations. We inflate Q slightly to 
    % keep the filter "awake" and prevent covariance collapse.
    
    tuningFactorVRW  = 1e0; % Inflate velocity uncertainty growth
    tuningFactorBias = 1e0; % Inflate bias random walk to track dynamics
    
    QV = (qVRW * tuningFactorVRW)  * eye(3);
    QB = (qARW * tuningFactorBias) * eye(3);
    
    % --- 3.3 Discrete-time Q Matrix Approximation ---
    % For a standard PV-EKF, process noise enters through velocity (accel 
    % integration) and bias states. Position noise is coupled via velocity.
    % Simplified diagonal discrete approximation: Q_d ≈ Q_c * dt
    Q11 = 1/3 * QV * dt^3; % Position noise variance
    Q12 = 1/2 * QV * dt^2; % Position-Velocity cross-correlation
    Q21 = Q12;
    Q22 = QV * dt;         % Velocity noise variance
    Q33 = QB * dt;         % Bias Random Walk variance
    
    O3 = zeros(3,3);
    
    EKF.Q = [Q11, Q12,  O3;
             Q21, Q22,  O3;
              O3,  O3, Q33];
    
    fprintf(' Process noise (Q_discrete) - Tuned:\n');
    fprintf('   - Position integrated: %.2e m²\n', EKF.Q(1,1));
    fprintf('   - VRW (Velocity)     : %.2e (m/s)²\n', EKF.Q(4,4));
    fprintf('   - ARW (Bias RW)      : %.2e (m/s²)²\n', EKF.Q(7,7));
    
    %% ===================================================================
    % 4. MEASUREMENT NOISE COVARIANCE (R)
    % ====================================================================
    % For GNSS, we have both high-frequency thermal noise (White) and 
    % slowly varying correlated errors (Gauss-Markov).
    % Since our 9-state EKF does not explicitly estimate the Gauss-Markov 
    % states (to save computational load), we MUST inflate the R matrix by 
    % adding the Gauss-Markov variance. This prevents the EKF from becoming 
    % overconfident during correlated ephemeris/ionosphere drifts.
    
    % Position measurement variance (White + GM)
    varPos = GNSS.sigmaPosWhite^2 + GNSS.sigmaPosGM^2;
    
    % Velocity measurement variance (White + GM)
    varVel = GNSS.sigmaVelWhite^2 + GNSS.sigmaVelGM^2;
    
    RPos = varPos * eye(3);
    RVel = varVel * eye(3);
    
    EKF.R_GNSS = blkdiag(RPos, RVel);
    
    fprintf('\n Measurement noise (R) - GNSS (White + GM Inflated):\n');
    fprintf('   - Position: %.2f m² (Sigma: %.2f m)\n', EKF.R_GNSS(1,1), sqrt(EKF.R_GNSS(1,1)));
    fprintf('   - Velocity: %.4f (m/s)² (Sigma: %.3f m/s)\n', EKF.R_GNSS(4,4), sqrt(EKF.R_GNSS(4,4)));
    
    fprintf('=== PV-EKF Initialization Complete ===\n');

end