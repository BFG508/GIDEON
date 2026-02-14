function EKF = initializeEKF_Nav(IMU, GNSS)
%==========================================================================
% initializeEKF_Nav - Initialize Position & Velocity Extended Kalman Filter
%                     (PV-EKF) for spacecraft translational navigation.
%
% Inputs:
%   IMU          - IMU parameter structure        (from initializeIMU.m)
%   GNSS         - GNSS parameter structure       (from initializeGNSS.m)
%   orbitalElems - Orbital elements structure     (from initializeOrbit.m)
%
% Outputs:
%   EKF - EKF state structure with fields:
%         .x      - Full state vector [9x1]: [r_ECI; v_ECI; biasAccel]
%         .P      - Error covariance matrix [9x9]
%         .Q      - Discrete process noise covariance [9x9]
%         .R_GNSS - GNSS measurement noise covariance [6x6]
%
% PV-EKF State Definition:
%   x = [r (3x1); v (3x1); b_a (3x1)]
%   where:
%     r  : Position vector in ECI frame [m]
%     v  : Velocity vector in ECI frame [m/s]
%     b_a: Accelerometer dynamic bias in body frame [m/s²]
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
    
    sigma_r0 = 100.0; % [m]
    sigma_v0 = 10.0;  % [m/s]
    sigma_b0 = 0.05;  % [m/s²]
    
    P0_pos  = (sigma_r0^2) * eye(3);
    P0_vel  = (sigma_v0^2) * eye(3);
    P0_bias = (sigma_b0^2) * eye(3);
    
    EKF.P = blkdiag(P0_pos, P0_vel, P0_bias);
    
    %% ===================================================================
    % 3. PROCESS NOISE COVARIANCE (Q)
    % ====================================================================
    % Computes the continuous-time Power Spectral Density (PSD) of the 
    % process noise and discretizes it for the EKF propagation step.
    
    g0 = 9.80665; % Standard gravity [m/s²]
    
    % --- 3.1 Base Noise from IMU Datasheet ---
    % Velocity Random Walk (VRW): Wideband noise on acceleration measuring
    accelVRW_metric = (IMU.accel.VRW * 1e-6) * g0; % [m/s²/sqrt(Hz)]
    q_vrw = accelVRW_metric^2;                     % PSD: [(m/s²)²/Hz]
    
    % Acceleration Random Walk (ARW): Random walk driving the bias drift
    accelARW_metric = IMU.accel.ARW / sqrt(3600);  % [m/s²/sqrt(s)]
    q_arw = accelARW_metric^2;                     % PSD: [(m/s²)²/Hz]
    
    % --- 3.2 Process Noise Tuning Factors ---
    % Since the S/C operates in LEO, the orbital dynamics (gravity, drag) 
    % are well modeled, but density variations (NRLMSISE-00) and attitude 
    % coupling inject unmodeled perturbations. We inflate Q slightly to 
    % keep the filter "awake" and prevent covariance collapse.
    
    tuningFactorVRW  = 1e0; % Inflate velocity uncertainty growth
    tuningFactorBias = 1e0; % Inflate bias random walk to track dynamics
    
    Q_v = (q_vrw * tuningFactorVRW)  * eye(3);
    Q_b = (q_arw * tuningFactorBias) * eye(3);
    
    % --- 3.3 Discrete-time Q Matrix Approximation ---
    % For a standard PV-EKF, process noise enters through velocity (accel 
    % integration) and bias states. Position noise is coupled via velocity.
    % Simplified diagonal discrete approximation: Q_d ≈ Q_c * dt
    
    Q_pos_discrete = (1/3) * Q_v * dt^3; % Position integration noise
    Q_vel_discrete = Q_v * dt;           % Velocity integration noise
    Q_bias_discrete = Q_b * dt;          % Bias random walk noise
    
    EKF.Q = blkdiag(Q_pos_discrete, Q_vel_discrete, Q_bias_discrete);
    
    fprintf(' Process noise (Q_discrete) - Tuned:\n');
    fprintf('   - Position integr: %.2e m²\n', EKF.Q(1,1));
    fprintf('   - VRW (velocity):  %.2e (m/s)²\n', EKF.Q(4,4));
    fprintf('   - ARW (bias RW):   %.2e (m/s²)²\n', EKF.Q(7,7));
    
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
    var_pos = GNSS.sigmaPosWhite^2 + GNSS.sigmaPosGM^2;
    
    % Velocity measurement variance (White + GM)
    var_vel = GNSS.sigmaVelWhite^2 + GNSS.sigmaVelGM^2;
    
    R_pos = var_pos * eye(3);
    R_vel = var_vel * eye(3);
    
    EKF.R_GNSS = blkdiag(R_pos, R_vel);
    
    fprintf('\n Measurement noise (R) - GNSS (White + GM Inflated):\n');
    fprintf('   - Position: %.2f m² (Sigma: %.2f m)\n', EKF.R_GNSS(1,1), sqrt(EKF.R_GNSS(1,1)));
    fprintf('   - Velocity: %.4f (m/s)² (Sigma: %.3f m/s)\n', EKF.R_GNSS(4,4), sqrt(EKF.R_GNSS(4,4)));
    
    fprintf('=== PV-EKF Initialization Complete ===\n');

end