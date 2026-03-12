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
%         .x      - Full state vector [15x1]: [rECI; vECI; biasAccel; rGM; vGM]
%         .P      - Error covariance matrix [15x15]
%         .Q      - Discrete process noise covariance [15x15]
%         .R_GNSS - GNSS measurement noise covariance [6x6]
%
% PV-EKF State Definition:
%   x = [r (3x1); v (3x1); bA (3x1); rGM (3x1); vGM (3x1)]
%   where:
%     r  : Position vector in ECI frame [m]
%     v  : Velocity vector in ECI frame [m/s]
%     bA : Accelerometer dynamic bias in body frame [m/s²]
%     rGM: GNSS correlated position error (Gauss-Markov) [m]
%     vGM: GNSS correlated velocity error (Gauss-Markov) [m/s]
%==========================================================================

    fprintf('\n=== Initializing PV-EKF (Navigation) ===\n');

    % Extract sampling time
    dt = IMU.dt;
    
    %% ===================================================================
    % 1. INITIAL FULL STATE VECTOR (x)
    % ====================================================================
    % The filter states will be properly initialized in the main loop using 
    % the first valid GNSS measurement. Here we pre-allocate the 15x1 vector.
    EKF.x = zeros(15, 1);
    
    %% ===================================================================
    % 2. INITIAL ERROR COVARIANCE MATRIX (P0)
    % ====================================================================
    % Defines the initial uncertainty of our state estimate.
    % - Position: ~100 m initial uncertainty
    % - Velocity: ~10 m/s initial uncertainty
    % - Accel Bias: ~0.05 m/s² initial uncertainty
    % - GNSS GM: Initialized to steady-state variance of the Markov process
    
    sigmaR0 = 100.0; % [m]
    sigmaV0 = 10.0;  % [m/s]
    sigmaB0 = 0.05;  % [m/s²]
    
    P0Pos   = (sigmaR0^2) * eye(3);
    P0Vel   = (sigmaV0^2) * eye(3);
    P0Bias  = (sigmaB0^2) * eye(3);
    P0PosGM = (GNSS.sigmaPosGM^2) * eye(3);
    P0VelGM = (GNSS.sigmaVelGM^2) * eye(3);
    
    EKF.P = blkdiag(P0Pos, P0Vel, P0Bias, P0PosGM, P0VelGM);
    
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
    
    tuningFactorVRW  = 1e0;  % Inflate velocity uncertainty growth
    tuningFactorBias = 2e1;  % Inflate bias random walk to track dynamics
    tuningFactorGM   = 2e-1; % Tuning for the Gauss-Markov states
    
    QV = (qVRW * tuningFactorVRW)  * eye(3);
    QB = (qARW * tuningFactorBias) * eye(3);

    % Discrete process noise for the Gauss-Markov states
    % Q_discrete = 2 * sigma^2 / tau * dt
    Q_GM_Pos = tuningFactorGM * (2 * GNSS.sigmaPosGM^2 / GNSS.tauPos) * dt * eye(3);
    Q_GM_Vel = tuningFactorGM * (2 * GNSS.sigmaVelGM^2 / GNSS.tauVel) * dt * eye(3);
    
    % --- 3.3 Discrete-time Q Matrix Approximation ---
    % For a standard PV-EKF, process noise enters through velocity (accel 
    % integration) and bias states. Position noise is coupled via velocity.
    % Simplified diagonal discrete approximation: Q_d ≈ Q_c * dt
    Q11 = 1/3 * QV * dt^3;            % Position noise variance (pure kinematics)
    Q12 = 1/2 * QV * dt^2;            % Position-Velocity cross-correlation
    Q21 = Q12;
    Q22 = QV * dt;                    % Velocity noise variance (pure kinematics)
    Q33 = QB * dt;                    % Bias Random Walk variance

    O3 = zeros(3,3);
    
    EKF.Q = [Q11, Q12,  O3,       O3,       O3;
             Q21, Q22,  O3,       O3,       O3;
              O3,  O3, Q33,       O3,       O3;
              O3,  O3,  O3, Q_GM_Pos,       O3;
              O3,  O3,  O3,       O3, Q_GM_Vel];
    
    fprintf(' Process noise (Q_discrete) - Tuned:\n');
    fprintf('   - Position integrated: %.2e m²\n', EKF.Q(1,1));
    fprintf('   - VRW (Velocity)     : %.2e (m/s)²\n', EKF.Q(4,4));
    fprintf('   - ARW (Bias RW)      : %.2e (m/s²)²\n', EKF.Q(7,7));
    fprintf('   - GM Position        : %.2e m²\n', EKF.Q(10,10));
    fprintf('   - GM Velocity        : %.2e (m/s)²\n', EKF.Q(13,13));
    
    %% ===================================================================
    % 4. MEASUREMENT NOISE COVARIANCE (R)
    % ====================================================================
    % For GNSS, the high-frequency thermal noise (White) is considered.
    % The slowly varying correlated errors (Gauss-Markov) are explicitly 
    % estimated as states within the EKF, so we don't inflate R.
    
    % Position measurement variance (White only)
    varPos = GNSS.sigmaPosWhite^2;
    
    % Velocity measurement variance (White only)
    varVel = GNSS.sigmaVelWhite^2;
    
    RPos = varPos * eye(3);
    RVel = varVel * eye(3);
    
    EKF.R_GNSS = blkdiag(RPos, RVel);
    
    fprintf('\n Measurement noise (R) - GNSS (White ONLY):\n');
    fprintf('   - Position: %.2f m² (σ: %.2f m)\n', EKF.R_GNSS(1,1), sqrt(EKF.R_GNSS(1,1)));
    fprintf('   - Velocity: %.4f (m/s)² (σ: %.3f m/s)\n', EKF.R_GNSS(4,4), sqrt(EKF.R_GNSS(4,4)));
    
    fprintf('=== PV-EKF Initialization Complete ===\n');

end